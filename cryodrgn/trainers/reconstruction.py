"""Base classes for reconstruction model learning engines."""

import os
import pickle
import contextlib
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Union
import yaml
from datetime import datetime as dt
import logging

import numpy as np
import torch
from torch import nn
import cryodrgn.utils
from cryodrgn import __version__, ctf
from cryodrgn.dataset import ImageDataset, TiltSeriesData, make_dataloader
from cryodrgn.lattice import Lattice
from cryodrgn.models.cryodrgnai import MyDataParallel
from cryodrgn.trainers._base import BaseConfigurations, BaseTrainer

try:
    import apex.amp as amp  # type: ignore  # PYR01
except ImportError:
    pass


@dataclass
class ReconstructionModelConfigurations(BaseConfigurations):
    """The abstract base class for parameter sets used by volume learning models.

    Arguments
    ---------
    > inherited from `BaseConfigurations`:
        verbose     An integer specifiying the verbosity level for this engine, with
                    the default value of 0 generally specifying no/minimum verbosity.
        seed        A non-negative integer used to fix the stochasticity of the random
                    number generators used by this engine for reproducibility.
                    The default is to not fix stochasticity and thus use a different
                    random seed upon each run of the engine.
        test_installation   Only perform a smoke test that this module has been
                            installed correctly and exit immediately without running
                            anything if this boolean value is set to `True`.
                            Default is not to run this test.

    model       A label for the reconstruction algorithm to be used — must be either
                `hps` for cryoDRGN v3 models or `amort` for cryoDRGN-AI models.
    zdim        The dimensionality of the latent space of conformations.
                Thus zdim=0 for homogeneous models and zdim>0 for hetergeneous models.
    num_epochs  The total number of epochs to use when training the model, not including
                pretraining epoch(s).

    dataset     Label for the particle dataset to be used as input for the model.
                If used, remaining input parameters can be omitted.
    particles   Path to the stack of particle images to use as input for the model.
                Must be a (.mrcs/.txt/.star/.cs file).
    ctf         Path to the file storing contrast transfer function parameters used to
                process the input particle images.
    poses       Path to the input particle poses data (.pkl).
    datadir     Path prefix to particle stack if loading relative paths from
                a .star or .cs file.
    ind         Path to a numpy array saved as a .pkl used to filter input particles.

    pose_estimation     Whether to perform ab-initio reconstruction ("abinit"),
                        reconstruction using fixed poses ("fixed"), or reconstruction
                        with SGD refinement of poses ("refine").
                        Default is to use fixed poses if poses file is given
                        and ab-initio otherwise.

    load:       Load model from given weights.<epoch>.pkl output file saved from
                previous run of the engine.
                Can also be given as "latest", in which the latest saved epoch in the
                given output directory will be used.
    lazy        Whether to use lazy loading of data into memory in smaller batches.
                Necessary if input dataset is too large to fit into memory.

    batch_size      The number of input images to use at a time when updating
                    the learning algorithm.

    multigpu        Whether to use all available GPUs available on this machine.
                    The default is to use only one GPU.

    log_interval    Print a log message every `N` number of training images.
    checkpoint      Save model results to file every `N` training epochs.

    pe_type     Label for the type of positional encoding to use.
    pe_dim      Number of frequencies to use in the positional encoding (default: D/2).
    volume_domain   Representation to use in the volume
                    decoder ("hartley" or "fourier").
    """

    # A parameter belongs to this configuration set if and only if it has a type and a
    # default value defined here, note that children classes inherit these parameters
    model: str = "autodec"
    zdim: int = 8
    num_epochs: int = 30
    dataset: str = None
    particles: str = None
    ctf: str = None
    poses: str = None
    datadir: str = None
    ind: str = None
    labels: str = None
    # Whether to start with given poses to some degree, or do ab initio pose search
    pose_estimation: str = None
    no_trans: bool = False
    # loading checkpoints from previous experiment runs
    load: str = None
    load_poses: str = None
    # Using a lazy data loader to minimize memory usage, shuffling training minibatches
    # and controlling their size, applying parallel computation and data processing
    batch_size: int = 8
    batch_size_known_poses: int = None
    batch_size_sgd: int = None
    batch_size_hps: int = None
    lazy: bool = False
    shuffle: bool = False
    shuffler_size: int = 0
    amp: bool = True
    multigpu: bool = False
    max_threads: int = 16
    num_workers: int = 0
    # How often to print log messages and save trained model data
    log_interval: int = 1000
    checkpoint: int = 1
    # If using a cryo-ET dataset, description of dataset tilting parameters
    tilt: bool = False
    n_tilts: int = 11
    tilt_deg: int = 45
    angle_per_tilt: int = 3
    dose_per_tilt: float = 2.97
    # masking
    window: bool = True
    window_r: float = 0.85
    # data normalization
    invert_data: bool = True
    data_norm: float = None
    use_real: bool = False
    # how long to do pretraining for, and what type
    pretrain: int = 10000
    reset_optim_after_pretrain: bool = False
    # Other parameters used in learning
    weight_decay: float = 0.0
    learning_rate: float = 1e-4
    pose_learning_rate: float = 1e-4
    l_extent: float = 0.5
    l_start: int = 12
    l_end: int = 32
    t_extent: float = 20.0
    t_ngrid: int = 7
    t_xshift: int = 0
    t_yshift: int = 0
    hidden_layers: int = 3
    hidden_dim: int = None
    pe_type: str = "gaussian"
    pe_dim: int = None
    volume_domain: str = None
    activation: str = "relu"
    feat_sigma: float = 0.5
    base_healpy: int = 2
    volume_optim_type: str = "adam"
    pose_sgd_emb_type: str = "quat"
    verbose_time: bool = False
    tdim: int = None
    tlayers: int = None
    t_emb_dim: int = None

    def __post_init__(self) -> None:
        """Parsing given configuration parameter values and checking their validity."""
        super().__post_init__()

        if self.model not in {"autodec", "autoenc"}:
            raise ValueError(
                f"Given model `{self.model}` not in currently supported model types:\n"
                f"`autodec` ... cryoDRGN v4 (CryoDRGN-AI)\n"
                f"`autoenc` ... cryoDRGN v3,v2,v1 \n"
            )
        if not isinstance(self.zdim, int) or self.zdim < 0:
            raise ValueError(
                f"Given latent space dimension {self.zdim=} "
                f"is not zero (homogeneous reconstruction) "
                f"or a positive integer (heterogeneous reconstruction)!"
            )
        if not isinstance(self.num_epochs, int) or self.num_epochs <= 0:
            raise ValueError(
                f"Given number of training epochs {self.zdim=} "
                f"is not a positive integer!"
            )

        # handling different ways of specifying the input data, starting with a
        # file containing the data files
        paths_file = os.environ.get("CRYODRGN_DATASETS")
        data_paths = cryodrgn.utils.load_yaml(paths_file) if paths_file else None
        if self.dataset is not None:
            if os.path.exists(self.dataset):
                paths = cryodrgn.utils.load_yaml(self.dataset)

                # resolve paths relative to the dataset file if they look relative
                for k in list(paths):
                    if paths[k] and not os.path.isabs(paths[k]):
                        paths[k] = os.path.abspath(os.path.join(self.dataset, paths[k]))

            elif data_paths and self.dataset not in data_paths:
                raise ValueError(
                    f"Given dataset {self.dataset} is not a "
                    "label in the list of known datasets!"
                )

            elif data_paths is None:
                raise ValueError(
                    "To specify datasets using a label, first specify"
                    "a .yaml catalogue of datasets using the "
                    "environment variable $CRYODRGN_DATASETS!"
                )

            # you can also give the dataset as a label in the global dataset list
            else:
                paths = data_paths[self.dataset]

            for k, v in paths.items():
                setattr(self, k, v)

        if isinstance(self.ind, str):
            if not self.ind.isnumeric():
                self.ind = os.path.abspath(self.ind)

        if self.pose_estimation is None:
            self.pose_estimation = "fixed" if self.poses else "abinit"
        if self.pose_estimation == "refine":
            if not self.poses:
                raise ValueError(
                    "Specify an input file (poses=) if refining ground truth poses!"
                )
            if self.volume_domain != "hartley":
                raise ValueError("Need to use --domain hartley if doing pose SGD")
        elif self.pose_estimation == "fixed":
            if not self.poses:
                raise ValueError(
                    "Specify an input file (poses=) if using ground truth poses!"
                )

        if self.batch_size_known_poses is None:
            self.batch_size_known_poses = self.batch_size
        if self.batch_size_sgd is None:
            self.batch_size_sgd = self.batch_size
        if self.batch_size_hps is None:
            self.batch_size_hps = self.batch_size

    def _file_dict(self) -> dict[str, Union[str, dict[str, str]]]:
        """Organizing the parameter values for use in human-readable formats."""
        configs = super()._file_dict()

        dataset_args = dict(
            particles=self.particles,
            poses=self.poses,
            ctf=self.ctf,
            ind=self.ind,
            dataset=self.dataset,
            labels=self.labels,
            datadir=self.datadir,
            data_norm=self.data_norm,
            invert_data=self.invert_data,
            use_real=self.use_real,
            window=self.window,
            window_r=self.window_r,
            tilt=self.tilt,
            load=self.load,
            load_poses=self.load_poses,
        )
        model_args = dict(
            model=self.model,
            zdim=self.zdim,
            hidden_layers=self.hidden_layers,
            hidden_dim=self.hidden_dim,
            pose_estimation=self.pose_estimation,
            no_trans=self.no_trans,
            pe_type=self.pe_type,
            pe_dim=self.pe_dim,
            feat_sigma=self.feat_sigma,
            volume_domain=self.volume_domain,
            activation=self.activation,
            weight_decay=self.weight_decay,
            learning_rate=self.learning_rate,
            pose_learning_rate=self.pose_learning_rate,
            tdim=self.tdim,
            tlayers=self.tlayers,
            t_emb_dim=self.t_emb_dim,
            volume_optim_type=self.volume_optim_type,
            pose_sgd_emb_type=self.pose_sgd_emb_type,
            l_extent=self.l_extent,
            l_start=self.l_start,
            l_end=self.l_end,
            t_extent=self.t_extent,
            t_ngrid=self.t_ngrid,
            t_xshift=self.t_xshift,
            t_yshift=self.t_yshift,
            base_healpy=self.base_healpy,
            tilt=self.tilt,
            n_tilts=self.n_tilts,
            tilt_deg=self.tilt_deg,
            angle_per_tilt=self.angle_per_tilt,
            dose_per_tilt=self.dose_per_tilt,
        )
        train_args = dict(
            num_epochs=self.num_epochs,
            checkpoint=self.checkpoint,
            verbose_time=self.verbose_time,
            log_interval=self.log_interval,
            multigpu=self.multigpu,
            max_threads=self.max_threads,
            batch_size=self.batch_size,
            batch_size_known_poses=self.batch_size_known_poses,
            batch_size_sgd=self.batch_size_sgd,
            batch_size_hps=self.batch_size_hps,
            amp=self.amp,
            lazy=self.lazy,
            shuffle=self.shuffle,
            shuffler_size=self.shuffler_size,
            pretrain=self.pretrain,
            reset_optim_after_pretrain=False,
            num_workers=self.num_workers,
        )

        return dict(
            dataset_args=dataset_args,
            model_args=model_args,
            train_args=train_args,
            **configs,
        )


class ReconstructionModelTrainer(BaseTrainer, ABC):
    """Abstract base class for volume reconstruction model training engines.

    Arguments
    ---------
    configs (dict)      The raw configuration parameters for this engine.
                        Will be parsed by the engine's configuration class.

    Attributes
    ----------
    label               String used to refer to this engine for e.g. logging messages.
    configs (ReconstructionModelConfigurations)
                        The parsed parameter configs for this engine.

    outdir (str)        The path where output produced by the engine will be saved.
    logger (Logger)     Logging utility used to create info and warning messages.
    use_cuda (bool)     Whether we are using CUDA GPUs (or otherwise CPUs).
    data (ImageDataset)     The processed particle stack data used as model input.

    reconstruction_model (nn.Module)    The PyTorch learning algorithm used to
                                        reconstruct volumes from the input image stack.
    """

    configs: ReconstructionModelConfigurations
    config_cls = ReconstructionModelConfigurations
    label = "CryoDRGN Reconstruction Engine"
    activations = {"relu": nn.ReLU, "leaky_relu": nn.LeakyReLU}

    # options for optimizers to use
    optim_types = {"adam": torch.optim.Adam, "lbfgs": torch.optim.LBFGS}

    @abstractmethod
    def make_reconstruction_model(self, weights=None) -> nn.Module:
        pass

    def __init__(self, configs: dict[str, Any], outdir: str) -> None:
        super().__init__(configs, outdir)

        # set the device
        torch.manual_seed(self.configs.seed)
        self.use_cuda = torch.cuda.is_available()
        self.device = torch.device("cuda" if self.use_cuda else "cpu")
        self.activation = self.activations[self.configs.activation]

        if self.use_cuda:
            self.logger.info(f"Using GPU {self.device}")
            self.n_prcs = torch.cuda.device_count()
            self.logger.info(f"Number of available gpus: {self.n_prcs}")
        else:
            self.logger.warning(f"No GPUs detected, using {self.device} instead!")
            self.n_prcs = 1

        # load index filter
        if self.configs.ind is not None:
            if self.configs.ind.isnumeric():
                self.logger.info(f"Keeping the first {self.configs.ind} particles")
                self.ind = np.arange(int(self.configs.ind))
            else:
                self.logger.info(f"Filtering image dataset with {configs['ind']}")
                self.ind = cryodrgn.utils.load_pkl(self.configs.ind)

                if self.configs.tilt:
                    particle_ind = pickle.load(open(self.configs.ind, "rb"))
                    pt, tp = TiltSeriesData.parse_particle_tilt(self.configs.particles)
                    self.ind = TiltSeriesData.particles_to_tilts(pt, particle_ind)
        else:
            self.ind = None

        if self.configs.tilt:
            self.tilts = torch.tensor(
                cryodrgn.utils.xrot(self.configs.tilt_deg).astype(np.float32),
                device=self.device,
            )
        else:
            self.tilts = None

        # load dataset
        self.logger.info(f"Loading dataset from {self.configs.particles}")
        if self.configs.pose_estimation in {"fixed", "refine"}:
            use_poses = self.configs.poses
        else:
            use_poses = None

        if not self.configs.tilt:
            self.data = ImageDataset(
                mrcfile=self.configs.particles,
                lazy=self.configs.lazy,
                norm=self.configs.data_norm,
                invert_data=self.configs.invert_data,
                keepreal=self.configs.use_real,
                ind=self.ind,
                window=self.configs.window,
                datadir=self.configs.datadir,
                window_r=self.configs.window_r,
                max_threads=self.configs.max_threads,
                poses_gt_pkl=use_poses,
                device=self.device,
            )
            self.particle_count = self.data.N

        else:
            self.data = TiltSeriesData(
                tiltstar=self.configs.particles,
                ind=self.ind,
                ntilts=self.configs.n_tilts,
                angle_per_tilt=self.configs.angle_per_tilt,
                window_r=self.configs.window_r,
                datadir=self.configs.datadir,
                max_threads=self.configs.max_threads,
                dose_per_tilt=self.configs.dose_per_tilt,
                device=self.device,
                poses_gt_pkl=use_poses,
            )
            self.particle_count = self.data.Np

        self.image_count = self.data.N
        self.resolution = self.data.D

        if self.configs.ctf:
            if self.configs.use_real:
                raise NotImplementedError(
                    "Not implemented with real-space encoder."
                    "Use phase-flipped images instead"
                )

            self.logger.info(f"Loading ctf params from {self.configs.ctf}")
            ctf_params = ctf.load_ctf_for_training(
                self.resolution - 1, self.configs.ctf
            )

            if self.ind is not None:
                ctf_params = ctf_params[self.ind, ...]
            assert ctf_params.shape == (self.image_count, 8), ctf_params.shape

            if self.configs.tilt:
                ctf_params = np.concatenate(
                    (ctf_params, self.data.ctfscalefactor.reshape(-1, 1)),
                    axis=1,  # type: ignore
                )
                self.data.voltage = float(ctf_params[0, 4])

            self.ctf_params = torch.tensor(ctf_params, device=self.device)
            self.apix = float(self.ctf_params[0, 0])
        else:
            self.ctf_params = None
            self.apix = None

        self.logger.info("Building lattice...")
        self.lattice = Lattice(
            self.resolution, extent=self.configs.l_extent, device=self.device
        )

        self.logger.info("Initializing volume model...")
        self.reconstruction_model = self.make_reconstruction_model()
        self.reconstruction_model.to(self.device)
        self.logger.info(self.reconstruction_model)

        # parallelize
        if self.configs.multigpu and torch.cuda.device_count() > 1:
            self.logger.info(f"Using {torch.cuda.device_count()} GPUs!")
            self.configs.batch_size *= self.n_prcs
            self.logger.info(f"Increasing batch size to {self.configs.batch_size}")
        elif self.configs.multigpu:
            self.logger.warning(
                f"WARNING: --multigpu selected, but {torch.cuda.device_count()} "
                "GPUs detected"
            )

        if self.configs.multigpu:
            self.configs.batch_size_known_poses = (
                self.configs.batch_size_known_poses * self.n_prcs
            )
            self.configs.batch_size_hps = self.configs.batch_size_hps * self.n_prcs
            self.configs.batch_size_sgd = self.configs.batch_size_sgd * self.n_prcs

        cpu_count = os.cpu_count() or 1
        if self.configs.num_workers > cpu_count:
            self.logger.warning(f"Reducing workers to {cpu_count} cpus")
            self.num_workers = cpu_count
        else:
            self.num_workers = self.configs.num_workers

        self.data_iterator = make_dataloader(
            self.data,
            batch_size=self.configs.batch_size,
            shuffler_size=self.configs.shuffler_size,
            shuffle=self.configs.shuffle,
            num_workers=self.num_workers,
            seed=self.configs.seed,
        )

        self.optimizer_types = {"hypervolume": self.configs.volume_optim_type}
        self.optimizers = {
            "hypervolume": self.optim_types[self.optimizer_types["hypervolume"]](
                self.reconstruction_model.parameters(),
                lr=self.configs.learning_rate,
                weight_decay=self.configs.weight_decay,
            )
        }

        # Mixed precision training
        self.scaler = None
        if self.configs.amp:
            if self.configs.batch_size % 8 != 0:
                self.logger.warning(
                    f"Batch size {self.configs.batch_size} not divisible by 8 "
                    f"and thus not optimal for AMP training!"
                )
            if (self.data.D - 1) % 8 != 0:
                self.logger.warning(
                    f"Image size {self.data.D - 1} not divisible by 8 "
                    f"and thus not optimal for AMP training!"
                )
            if self.configs.zdim % 8 != 0:
                self.logger.warning(
                    f"Z dimension {self.configs.zdim} is not a multiple of 8 "
                    "-- AMP training speedup is not optimized!"
                )

            try:  # Mixed precision with apex.amp
                (
                    self.reconstruction_model,
                    self.optimizers["hypervolume"],
                ) = amp.initialize(
                    self.reconstruction_model,
                    self.optimizers["hypervolume"],
                    opt_level="O1",
                )
            except:  # noqa: E722
                # Mixed precision with pytorch (v1.6+)
                self.scaler = torch.cuda.amp.grad_scaler.GradScaler()

        if self.scaler is not None:
            try:
                self.amp_mode = torch.amp.autocast("cuda")
            except AttributeError:
                self.amp_mode = torch.cuda.amp.autocast_mode.autocast()
        else:
            self.amp_mode = contextlib.nullcontext()

        # counters used across training iterations
        self.current_epoch = None
        self.accum_losses = None
        self.total_batch_count = None
        self.epoch_batch_count = None
        self.epoch_losses = None
        self.total_images_seen = None
        self.epoch_images_seen = None
        self.conf_search_particles = None
        self.beta = None
        self.epoch_start_time = None
        self.base_pose = None
        self.base_poses = list()
        self.predicted_rots = None
        self.predicted_trans = None
        self.predicted_conf = None

        # Initializing model from a saved checkpoint if given by user
        self.checkpoint = None
        if self.configs.load:
            if self.configs.load == "latest" or self.configs.load is True:
                load_path = cryodrgn.utils.find_latest_output(self.outdir)
            else:
                load_path = self.configs.load

            self.logger.info(f"Loading checkpoint from `{load_path}` ...")
            self.checkpoint = torch.load(
                load_path, weights_only=False, map_location=self.device
            )
            state_dict = self.checkpoint["model_state_dict"]

            # With ab-initio reconstruction we also need to load the poses saved by
            # the model; these are looked for in the same directory as the checkpoint
            if self.configs.pose_estimation == "abinit":
                if self.configs.load_poses is None or self.configs.load_poses is True:
                    epoch = load_path.split(".")[-2]
                    load_poses_path = os.path.join(self.outdir, f"pose.{epoch}.pkl")
                    if not os.path.exists(load_poses_path):
                        load_poses_path = os.path.joqin(
                            os.path.dirname(load_path), f"pose.{epoch}.pkl"
                        )
                else:
                    load_poses_path = self.configs.load_poses

                assert os.path.exists(load_poses_path)
                self.logger.info(f"Loading poses from {load_poses_path}...")
                rot, trans = cryodrgn.utils.load_pkl(load_poses_path)
                if np.any(trans > 1):
                    raise RuntimeError(
                        "Old pose format detected."
                        "Translations must be in units of fraction of box."
                    )

                # Convert translations to pixel units to feed back to the model
                self.predicted_rots = rot
                self.predicted_trans = trans * self.model_resolution

            if "base_shifts" in state_dict:
                state_dict.pop("base_shifts")

            self.logger.info(
                self.reconstruction_model.load_state_dict(state_dict, strict=False)
            )
            if "output_mask_radius" in self.checkpoint:
                self.reconstruction_model.output_mask.update_radius(
                    self.checkpoint["output_mask_radius"]
                )
            self.start_epoch = self.checkpoint["epoch"] + 1
        else:
            self.start_epoch = 1

    def train(self) -> None:
        """Train the model to reconstruct volumes from the input particle stack."""

        paths_file = os.environ.get("CRYODRGN_DATASETS")
        if self.configs.particles is None or not os.path.isfile(self.configs.particles):
            raise ValueError(
                "Must specify either a dataset label stored in "
                f"{paths_file} or the path to an existing particles file to train!"
            )

        os.makedirs(self.outdir, exist_ok=True)
        self.save_configs()
        if self.configs.verbose:
            self.logger.setLevel(logging.DEBUG)

        self.logger.addHandler(
            logging.FileHandler(os.path.join(self.outdir, "run.log"))
        )
        self.logger.info(f"cryoDRGN {__version__}")
        self.logger.info(f"Model Configurations:\n{str(self.configs)}")

        if self.configs.pose_estimation != "fixed" and not self.configs.load:
            self.predicted_rots = np.empty((self.image_count, 3, 3))
            self.predicted_trans = (
                np.empty((self.image_count, 2)) if not self.configs.no_trans else None
            )
        self.predicted_conf = (
            np.empty((self.particle_count, self.configs.zdim))
            if self.configs.zdim > 0
            else None
        )

        train_start_time = dt.now()
        self.accum_losses = dict()
        if self.do_pretrain:
            self.current_epoch = 0
            self.pretrain()

        self.logger.info("--- Training Starts Now ---")
        self.current_epoch = self.start_epoch
        self.total_batch_count = 0
        self.total_images_seen = 0
        self.conf_search_particles = 0

        while self.current_epoch <= self.configs.num_epochs:
            self.epoch_start_time = dt.now()
            self.epoch_batch_count = 0
            self.epoch_images_seen = 0
            for k in self.accum_losses:
                self.accum_losses[k] = 0

            will_make_checkpoint = self.will_make_checkpoint
            if will_make_checkpoint:
                self.logger.info("Will make a full summary at the end of this epoch")

            self.begin_epoch()

            for batch in self.data_iterator:
                self.total_batch_count += 1
                self.epoch_batch_count += 1
                len_y = len(batch["indices"])
                self.total_images_seen += len_y
                self.epoch_images_seen += len_y
                self.conf_search_particles += len_y

                self.reconstruction_model.train()
                self.optimizers["hypervolume"].zero_grad()
                with self.amp_mode:
                    (
                        losses,
                        tilt_ind,
                        ind,
                        rot,
                        trans,
                        z_mu,
                        z_logvar,
                    ) = self.train_batch(batch)

                if isinstance(rot, torch.Tensor):
                    rot = rot.detach().cpu().numpy()
                if isinstance(trans, torch.Tensor):
                    trans = trans.detach().cpu().numpy()

                ind_tilt = tilt_ind if tilt_ind is not None else ind
                if rot is not None and self.predicted_rots is not None:
                    self.predicted_rots[ind_tilt] = rot.reshape(-1, 3, 3)
                if trans is not None and self.predicted_trans is not None:
                    self.predicted_trans[ind_tilt] = trans.reshape(-1, 2)
                if self.base_pose is not None and self.base_poses is not None:
                    self.base_poses.append((ind.cpu().numpy(), self.base_pose))

                # logging
                for loss_k, loss_val in losses.items():
                    loss_v = (
                        loss_val if isinstance(loss_val, float) else loss_val.item()
                    )
                    if loss_k in self.accum_losses:
                        self.accum_losses[loss_k] += loss_v * len(ind)
                    else:
                        self.accum_losses[loss_k] = loss_v * len(ind)

                # scalar summary
                if self.epoch_images_seen % self.configs.log_interval < len_y:
                    self.print_batch_summary(losses)

                if self.epoch_images_seen > self.data.N:
                    break

            self.end_epoch()
            self.print_epoch_summary()

            if self.configs.verbose_time:
                torch.cuda.synchronize()

            # image and pose summary
            if will_make_checkpoint:
                self.save_epoch_data()

            self.current_epoch += 1

        t_total = dt.now() - train_start_time
        self.logger.info(
            f"Finished in {t_total} ({t_total / self.configs.num_epochs} per epoch)"
        )

    @property
    def will_make_checkpoint(self) -> bool:
        """Is this engine in a training epoch where model results will be saved?"""

        make_chk = self.current_epoch == self.configs.num_epochs
        make_chk |= (
            self.configs.checkpoint
            and self.current_epoch % self.configs.checkpoint == 0
        )
        make_chk |= self.in_pose_search_step

        return make_chk

    @property
    def model_module(self) -> nn.Module:
        if isinstance(self.reconstruction_model, (nn.DataParallel, MyDataParallel)):
            model_mdl = self.reconstruction_model.module
        else:
            model_mdl = self.reconstruction_model

        return model_mdl

    @property
    def model_resolution(self) -> int:
        return (
            self.model_module.lattice.D
            if self.configs.zdim > 0
            else self.model_module.D
        )

    def pretrain(self) -> None:
        """Iterate the model before main learning epochs for better initializations."""

        self.logger.info("Will make a full summary at the end of this epoch")
        self.logger.info(f"Will pretrain on {self.configs.pretrain} particles")
        self.epoch_start_time = dt.now()
        self.epoch_images_seen = 0
        self.accum_losses = {"total": 0.0}

        pretrain_dataloader = make_dataloader(
            self.data,
            batch_size=self.configs.batch_size_known_poses,
            num_workers=self.configs.num_workers,
            shuffler_size=self.configs.shuffler_size,
            shuffle=self.configs.shuffle,
            seed=self.configs.seed,
        )
        epoch_label = f"Pretrain Epoch [{self.current_epoch}/{self.configs.num_epochs}]"

        for batch in pretrain_dataloader:
            len_y = len(batch["indices"])
            self.epoch_images_seen += len_y
            self.pretrain_batch(batch)

            if self.configs.verbose_time:
                torch.cuda.synchronize()

            if self.epoch_images_seen % self.configs.log_interval < len_y:
                self.logger.info(
                    f"{epoch_label}; {self.epoch_images_seen} images seen; "
                    f"avg. loss={self.average_losses['total']:.4g}"
                )

            if self.epoch_images_seen >= self.configs.pretrain:
                break

        self.print_epoch_summary()
        self.save_epoch_data()

        # Reset the model after pretraining if asked for
        if self.configs.reset_optim_after_pretrain:
            self.logger.info(">> Resetting optimizer after pretrain")
            self.optimizers["hypervolume"] = torch.optim.Adam(
                self.reconstruction_model.parameters(),
                lr=self.configs.learning_rate,
                weight_decay=self.configs.weight_decay,
            )

    def get_configs(self) -> dict[str, Any]:
        """Retrieves all given and inferred configurations for downstream use.

        Note that we need this in addition to the class attribute `configs` which is a
        `ReconstructionModelConfigurations` object in order to 1) define a way of saving
        these configurations to file in a more structured, hierarchical, human-readable
        format where the configurations are sorted and partitioned according to the
        facet of the engine they concern and 2) provide a method of retrieving
        configuration values that need to be computed from the particle stack (such as
        `data.norm`) without having to reload the stack from file.
        """
        configs = self.configs.file_dict
        lattice_args = dict(
            D=self.lattice.D,
            l_extent=self.lattice.extent,
            ignore_DC=self.lattice.ignore_DC,
        )

        # Properties inferred from the input dataset only once it has been loaded
        dataset_args = dict(
            norm=self.data.norm,
            apix=self.apix,
            n_particles=self.particle_count,
            n_images=self.image_count,
        )
        configs["dataset_args"].update(dataset_args)
        configs["lattice_args"] = lattice_args

        return configs

    def save_configs(self) -> None:
        """Saves all given and inferred configurations to file."""

        configs = self.get_configs()
        if "version" not in configs:
            configs["version"] = cryodrgn.__version__
        if "time" not in configs:
            configs["time"] = dt.now()

        with open(os.path.join(self.outdir, "config.yaml"), "w") as f:
            yaml.dump(configs, f, default_flow_style=False, sort_keys=False)

    def begin_epoch(self) -> None:
        """Actions to perform at the beginning of each training epoch."""
        pass

    def end_epoch(self) -> None:
        """Actions to perform at the end of each training epoch."""
        pass

    @abstractmethod
    def train_batch(self, batch: dict[str, torch.Tensor]) -> tuple:
        pass

    def pretrain_batch(self, batch: dict[str, torch.Tensor]) -> None:
        """Process a single batch of data during a pretraining epoch of the model.

        This method is called during the pretraining phase of model training, which
        occurs before the main training epochs. It uses the same training logic as the
        main training phase by delegating to the train_batch method, but is called with
        a different dataloader configuration optimized for pretraining.

        Parameters
        ----------
        batch : dict[str, torch.Tensor]
            A dictionary containing the batch data tensors needed for training. The
            exact contents depend on the specific model implementation.
        """
        self.train_batch(batch)

    @abstractmethod
    def save_epoch_data(self) -> None:
        """Save current model state and results at the conclusion of a training epoch.

        This method is called at the end of each training epoch for which checkpointing
        is enabled. It saves the current state of the model to files in the output
        directory, including:
            - Model weights and parameters (weights.*.pkl)
            - Predicted poses and translations for ab-initio models (poses.*.pkl)
            - Latent space embeddings for heterogeneous models (z.*.pkl)
            - A reconstructed volume, chosen from within the latent space if
              necessary (reconstruct.*.mrc)

        The format and contents of these saved files should be consistent across
        reconstruction model implementations.

        """
        pass

    def print_batch_summary(self, losses: dict[str, float]) -> None:
        """Create a summary at the end of a training batch and print it to the log."""
        raise NotImplementedError

    def print_epoch_summary(self) -> None:
        """Create a summary at the end of a training epoch and print it to the log."""
        self.configs: ReconstructionModelConfigurations
        epoch_lbl = f"[{self.current_epoch}/{self.configs.num_epochs}]"

        if self.configs.zdim > 0:
            loss_str = ", ".join(
                [
                    f"{loss_k} = {loss_val:.4g}"
                    for loss_k, loss_val in self.average_losses.items()
                    if loss_k != "total"
                ]
            )
        else:
            loss_str = ""

        time_str, mcrscd_str = str(dt.now() - self.epoch_start_time).split(".")
        time_str = ".".join([time_str, mcrscd_str[:3]])
        log_msg = (
            f"===> Training Epoch {epoch_lbl} === Finished in {time_str}\n"
            f"\t\t\t\t\t\t             |> "
            f"Avg. Losses: Total = {self.average_losses['total']:.4g}\n"
        )
        if loss_str:
            log_msg += f"\t\t\t\t\t\t             |>              {loss_str}"

        self.logger.info(log_msg)

    @property
    def in_pose_search_step(self) -> bool:
        """Whether we are in a pose search epoch of the training stage."""
        return False

    @property
    def in_pretraining(self) -> bool:
        """Whether we are in a pretraining epoch of the training stage."""
        return self.current_epoch is not None and self.current_epoch == 0

    @property
    def epoch_type(self) -> str:
        """A human-readable label for the type of the current training epoch."""
        if self.in_pretraining:
            etype = "pre"
        elif self.in_pose_search_step:
            etype = "hps"
        else:
            etype = "fix"

        return etype

    @property
    def epoch_lbl(self) -> str:
        """A human-readable label for the current training epoch."""
        return ".".join([self.epoch_type, str(self.current_epoch)])

    @property
    def average_losses(self) -> dict[str, float]:
        """Calculate the running mean losses for the current training epoch."""

        return {
            loss_k: (
                loss_val / self.epoch_images_seen
                if loss_k in {"total", "gen", "kld"}
                else loss_val
            )
            for loss_k, loss_val in self.accum_losses.items()
        }
