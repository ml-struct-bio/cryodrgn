"""Reconstructing volume(s) from picked cryoEM and cryoET particles."""

import argparse
import os
import pickle
import yaml
import logging
from datetime import datetime as dt
import numpy as np
import time
from typing_extensions import Any
from types import SimpleNamespace
import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader

from cryodrgn.mrcfile import write_mrc
from cryodrgn.lattice import Lattice
from cryodrgn import utils, dataset, ctf
from cryodrgn.losses import kl_divergence_conf, l1_regularizer, l2_frequency_bias
from cryodrgn.models_ai import DrgnAI, MyDataParallel
from cryodrgn.masking import CircularMask, FrequencyMarchingMask
from cryodrgn.commands_drgnai.analyze import ModelAnalyzer


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "particles",
        type=os.path.abspath,
        help="Input particles (.mrcs, .star, .cs, or .txt)",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=os.path.abspath,
        required=True,
        help="Working directory containing out/ for outputs",
    )

    # Basic training controls
    parser.add_argument(
        "--load", type=str, help="Load a previous checkpoint weights.pkl"
    )
    parser.add_argument("--seed", type=int, default=np.random.randint(0, 100000))
    parser.add_argument("-v", "--verbose", action="store_true")

    # Dataset loading
    parser.add_argument("--ctf", type=os.path.abspath, help="CTF parameters (.pkl)")
    parser.add_argument(
        "--datadir", type=os.path.abspath, help="Path prefix for star/cs images"
    )
    parser.add_argument("--ind", type=os.path.abspath, help="Subset indices pkl to use")
    parser.add_argument("--relion31", action="store_true")
    parser.add_argument("--uninvert-data", dest="invert_data", action="store_false")
    parser.add_argument("--lazy", action="store_true")
    parser.add_argument("--max-threads", type=int, default=16)

    # Logging
    parser.add_argument("--log-interval", type=int, default=10000)
    parser.add_argument("--log-heavy-interval", type=int, default=5)
    parser.add_argument("--verbose-time", action="store_true")

    # Data loading and parallelism
    parser.add_argument("--no-shuffle", dest="shuffle", action="store_false")
    parser.add_argument("--num-workers", type=int, default=2)
    parser.add_argument("--fast-dataloading", action="store_true")
    parser.add_argument("--shuffler-size", type=int, default=32768)
    parser.add_argument("--multigpu", action="store_true")
    parser.add_argument(
        "--no-amp",
        action="store_false",
        dest="amp",
        help="Disable torch.amp mixed precision",
    )

    # Batch sizes
    parser.add_argument("--batch-size-hps", type=int, default=8)
    parser.add_argument("--batch-size-sgd", type=int, default=32)
    parser.add_argument("--batch-size-known-poses", type=int, default=128)

    # Optimizers
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--lr-pose-table", type=float, default=1e-3)
    parser.add_argument("--lr-conf-table", type=float, default=1e-2)
    parser.add_argument("--lr-conf-encoder", type=float, default=1e-4)
    parser.add_argument("--wd", type=float, default=0.0)
    parser.add_argument(
        "--hypervolume-optimizer-type", choices=("adam",), default="adam"
    )
    parser.add_argument(
        "--pose-table-optimizer-type", choices=("adam", "lbfgs"), default="adam"
    )
    parser.add_argument(
        "--conf-table-optimizer-type", choices=("adam", "lbfgs"), default="adam"
    )
    parser.add_argument(
        "--conf-encoder-optimizer-type", choices=("adam",), default="adam"
    )

    # Scheduling
    parser.add_argument("--n-imgs-pose-search", type=int, default=500000)
    parser.add_argument("--epochs-sgd", type=int, default=100)
    parser.add_argument("--pose-only-phase", type=int, default=0)

    # Masking
    parser.add_argument(
        "--output-mask", choices=("circ", "frequency_marching"), default="circ"
    )
    parser.add_argument("--add-one-frequency-every", type=int, default=100000)
    parser.add_argument("--n-frequencies-per-epoch", type=int, default=10)
    parser.add_argument("--max-freq", type=int)
    parser.add_argument("--window-radius-gt-real", type=float, default=0.85)

    # Losses
    parser.add_argument("--beta-conf", type=float, default=0.0)
    parser.add_argument("--trans-l1-regularizer", type=float, default=0.0)
    parser.add_argument("--l2-smoothness-regularizer", type=float, default=0.0)

    # Z / heterogeneity
    parser.add_argument("--variational-het", action="store_true")
    parser.add_argument("--zdim", type=int, default=4)
    parser.add_argument("--std-z-init", type=float, default=0.1)
    parser.add_argument("--use-conf-encoder", action="store_true")
    parser.add_argument("--depth-cnn", type=int, default=5)
    parser.add_argument("--channels-cnn", type=int, default=32)
    parser.add_argument("--kernel-size-cnn", type=int, default=3)
    parser.add_argument("--resolution-encoder", type=int)

    # Hypervolume
    parser.add_argument("--explicit-volume", action="store_true")
    parser.add_argument("--hypervolume-layers", type=int, default=3)
    parser.add_argument("--hypervolume-dim", type=int, default=256)
    parser.add_argument("--pe-type", choices=("gaussian",), default="gaussian")
    parser.add_argument("--pe-dim", type=int, default=64)
    parser.add_argument("--feat-sigma", type=float, default=0.5)
    parser.add_argument("--hypervolume-domain", choices=("hartley",), default="hartley")
    parser.add_argument("--pe-type-conf", choices=(None, "geom"), default=None)
    parser.add_argument(
        "--initial-conf", type=os.path.abspath, help="Initial conformations (.pkl)"
    )

    # Pretrain
    parser.add_argument("--n-imgs-pretrain", type=int, default=10000)

    # Pose search
    parser.add_argument("--l-start", type=int, default=12)
    parser.add_argument("--l-end", type=int, default=32)
    parser.add_argument("--n-iter", type=int, default=4)
    parser.add_argument("--t-extent", type=float, default=20.0)
    parser.add_argument("--t-n-grid", type=int, default=7)
    parser.add_argument("--t-x-shift", type=float, default=0.0)
    parser.add_argument("--t-y-shift", type=float, default=0.0)
    parser.add_argument("--no-trans-search-at-pose-search", action="store_true")
    parser.add_argument("--n-kept-poses", type=int, default=8)
    parser.add_argument("--base-healpy", type=int, default=2)
    parser.add_argument("--no-trans", action="store_true")

    parser.add_argument(
        "--norm",
        type=float,
        nargs=2,
        default=None,
        help="Data normalization as shift, 1/scale (default: mean, std of dataset)",
    )

    parser.add_argument(
        "--no-analysis",
        dest="do_analysis",
        action="store_false",
        help="Do not run analysis on the final training epoch",
    )


class ModelTrainer:
    """An engine for training the DRGN-AI reconstruction model on particle data.

    The two key methods of this engine class are the `__init__()` method, in which
    model parameters and data structures are initialized, and `train()`, in which the
    model is trained in batches over the particle input data.

    Attributes
    ----------
    configs (TrainingConfigurations)
        Values of all user-set parameters controlling the behaviour of the model.

    outdir (str):   Folder `out/` within the experiment working directory where
                    model results will be saved.

    n_particles_dataset (int):  The number of picked particles in the data.
    pretraining (bool):     Whether we are in the pretraining stage.
    epoch (int):    Which training epoch the model is presently in.

    logger (logging.Logger):    Utility for printing and writing information
                                about the model as it is running.
    """

    # options for optimizers to use
    optim_types = {"adam": torch.optim.Adam, "lbfgs": torch.optim.LBFGS}

    # placeholders for runtimes
    run_phases = [
        "dataloading",
        "to_gpu",
        "ctf",
        "encoder",
        "decoder",
        "decoder_coords",
        "decoder_query",
        "loss",
        "backward",
        "to_cpu",
    ]

    def make_dataloader(self, batch_size: int) -> DataLoader:
        return dataset.make_dataloader(
            self.data,
            batch_size=batch_size,
            num_workers=self.configs.num_workers,
            shuffle=self.configs.shuffle,
            seed=self.configs.seed,
        )

    def __init__(self, outdir: str, config_vals: dict[str, Any]) -> None:
        """Initialize model parameters and variables.

        Arguments
        ---------
        outdir:         Location on file where model results will be saved.
        config_vals:    Parsed model parameter values provided by the user.
        """

        self.logger = logging.getLogger(__name__)
        self.outdir = outdir
        self.configs = SimpleNamespace(**config_vals)

        # Create the output folder for model results and log file for model training
        os.makedirs(self.outdir, exist_ok=True)
        self.logger.addHandler(
            logging.FileHandler(os.path.join(self.outdir, "run.log"))
        )

        # Parallelize training across GPUs if --multigpu config option is selected
        gpu_count = torch.cuda.device_count()
        if self.configs.multigpu and gpu_count > 1:
            self.n_prcs = int(gpu_count)
            self.logger.info(f"Using {gpu_count} GPUs!")
            if self.configs.batch_size_known_poses is not None:
                new_batch_size = self.configs.batch_size_known_poses * self.n_prcs
                self.logger.info(
                    f"Increasing batch size for known poses to {new_batch_size}"
                )
                self.configs.batch_size_known_poses = new_batch_size
            if self.configs.batch_size_hps is not None:
                new_batch_size = self.configs.batch_size_hps * self.n_prcs
                self.logger.info(f"Increasing batch size for HPS to {new_batch_size}")
                self.configs.batch_size_hps = new_batch_size
            if self.configs.batch_size_sgd is not None:
                new_batch_size = self.configs.batch_size_sgd * self.n_prcs
                self.logger.info(f"Increasing batch size for SGD to {new_batch_size}")
                self.configs.batch_size_sgd = new_batch_size
        elif self.configs.multigpu:
            self.n_prcs = 1
            self.logger.warning(
                f"--multigpu selected, but only {gpu_count} GPUs detected!"
            )
        elif gpu_count > 1:
            self.n_prcs = 1
            self.logger.warning(
                f"Using one GPU in spite of {gpu_count} available GPUs "
                f"because --multigpu is not being used!"
            )
        else:
            self.n_prcs = 1

        np.random.seed(self.configs.seed)
        torch.manual_seed(self.configs.seed)

        # Set the compute device, using the first available GPU
        self.use_cuda = torch.cuda.is_available()
        self.device = torch.device("cuda:0" if self.use_cuda else "cpu")
        self.logger.info(f"Use cuda {self.use_cuda}")

        # Load the index used to filter particles, if given
        if self.configs.ind is not None:
            if isinstance(self.configs.ind, int):
                self.logger.info(f"Keeping {self.configs.ind} particles")
                self.index = np.arange(self.configs.ind)

            elif isinstance(self.configs.ind, str):
                if not os.path.exists(self.configs.ind):
                    raise ValueError(
                        "Given subset index file "
                        f"`{self.configs.ind}` does not exist!"
                    )

                self.logger.info(f"Filtering dataset with {self.configs.ind}")
                self.index = pickle.load(open(self.configs.ind, "rb"))

        else:
            self.index = None

        self.logger.info("Creating dataset")
        if self.configs.norm is None:
            norm_mean, norm_std = None, None
        elif isinstance(self.configs.norm, (list, tuple)):
            norm_mean, norm_std = self.configs.norm
        else:
            norm_mean, norm_std = self.configs.norm, 1.0

        if norm_mean is not None and norm_std is not None:
            data_norm = (norm_mean, norm_std)
        elif norm_mean is not None:
            data_norm = (norm_mean, 1.0)
        elif norm_std is not None:
            data_norm = (0.0, norm_std)
        else:
            data_norm = None
        if data_norm is not None:
            self.logger.info(
                f"Manually overriding data normalization: (mean, std) = {data_norm}"
            )

        self.data = dataset.ImageDataset(
            self.configs.particles,
            norm=data_norm,
            keepreal=True,
            invert_data=self.configs.invert_data,
            ind=self.index,
            window_r=self.configs.window_radius_gt_real,
            max_threads=self.configs.max_threads,
            lazy=self.configs.lazy,
            datadir=self.configs.datadir,
        )
        self.n_particles_dataset = self.data.N
        self.n_tilts_dataset = self.data.N
        self.resolution = self.data.D

        # Load contrast transfer function parameters, if given
        if self.configs.ctf is not None:
            self.logger.info(f"Loading ctf params from {self.configs.ctf}")
            ctf_params = ctf.load_ctf_for_training(
                self.resolution - 1, self.configs.ctf
            )

            if self.index is not None:
                self.logger.info("Filtering dataset")
                ctf_params = ctf_params[self.index]

            assert ctf_params.shape == (self.n_tilts_dataset, 8)
            self.ctf_params = torch.tensor(ctf_params)
            self.ctf_params = self.ctf_params.to(self.device)

        else:
            self.ctf_params = None

        self.logger.info("Building lattice")
        self.lattice = Lattice(self.resolution, extent=0.5, device=self.device)

        # Set up the output mask
        if self.configs.output_mask == "circ":
            radius = (
                self.lattice.D // 2
                if self.configs.max_freq is None
                else self.configs.max_freq
            )
            self.output_mask = CircularMask(self.lattice, radius)

        elif self.configs.output_mask == "frequency_marching":
            self.output_mask = FrequencyMarchingMask(
                self.lattice,
                self.lattice.D // 2,
                radius=self.configs.l_start_fm,
                add_one_every=self.configs.add_one_frequency_every,
            )

        else:
            raise NotImplementedError

        # Set up the pose search parameters
        self.epochs_pose_search = max(
            2, self.configs.n_imgs_pose_search // self.n_particles_dataset + 1
        )
        ps_params = {
            "l_min": self.configs.l_start,
            "l_max": self.configs.l_end,
            "t_extent": self.configs.t_extent,
            "t_n_grid": self.configs.t_n_grid,
            "niter": self.configs.n_iter,
            "nkeptposes": self.configs.n_kept_poses,
            "base_healpy": self.configs.base_healpy,
            "t_xshift": self.configs.t_x_shift,
            "t_yshift": self.configs.t_y_shift,
            "no_trans_search_at_pose_search": self.configs.no_trans_search_at_pose_search,
            "tilting_func": None,
        }

        # CNN
        cnn_params = {
            "conf": self.configs.use_conf_encoder,
            "depth_cnn": self.configs.depth_cnn,
            "channels_cnn": self.configs.channels_cnn,
            "kernel_size_cnn": self.configs.kernel_size_cnn,
        }

        # Conformational encoder
        if self.configs.zdim > 0:
            self.logger.info(
                "Heterogeneous reconstruction with " f"zdim = {self.configs.zdim}"
            )
        else:
            self.logger.info("Homogeneous reconstruction")

        conf_regressor_params = {
            "z_dim": self.configs.zdim,
            "std_z_init": self.configs.std_z_init,
            "variational": self.configs.variational_het,
        }

        # Hypervolume
        hyper_volume_params = {
            "explicit_volume": self.configs.explicit_volume,
            "n_layers": self.configs.hypervolume_layers,
            "hidden_dim": self.configs.hypervolume_dim,
            "pe_type": self.configs.pe_type,
            "pe_dim": self.configs.pe_dim,
            "feat_sigma": self.configs.feat_sigma,
            "domain": self.configs.hypervolume_domain,
            "extent": self.lattice.extent,
            "pe_type_conf": self.configs.pe_type_conf,
        }

        will_use_point_estimates = self.configs.epochs_sgd >= 1
        self.logger.info("Initializing model...")

        self.model = DrgnAI(
            self.lattice,
            self.output_mask,
            self.n_particles_dataset,
            self.n_tilts_dataset,
            cnn_params,
            conf_regressor_params,
            hyper_volume_params,
            resolution_encoder=self.configs.resolution_encoder,
            no_trans=self.configs.no_trans,
            use_gt_poses=False,
            use_gt_trans=False,
            will_use_point_estimates=will_use_point_estimates,
            ps_params=ps_params,
            verbose_time=self.configs.verbose_time,
            pretrain_with_gt_poses=False,
        )

        # Initialization from a checkpoint saved to file from a previous training run
        if self.configs.load:
            self.logger.info(f"Loading checkpoint from {self.configs.load}")
            checkpoint = torch.load(self.configs.load)
            state_dict = checkpoint["model_state_dict"]

            if "base_shifts" in state_dict:
                state_dict.pop("base_shifts")

            self.logger.info(self.model.load_state_dict(state_dict, strict=False))
            self.start_epoch = checkpoint["epoch"] + 1

            if "output_mask_radius" in checkpoint:
                self.output_mask.update_radius(checkpoint["output_mask_radius"])

        else:
            self.start_epoch = 0 if self.configs.n_imgs_pretrain > 0 else 1

        # Move to GPU and parallelize the model if necessary
        self.logger.info(self.model)
        parameter_count = sum(
            p.numel() for p in self.model.parameters() if p.requires_grad
        )
        self.logger.info(f"{parameter_count} parameters in model")

        # TODO: Replace with DistributedDataParallel
        if self.n_prcs > 1:
            self.model = MyDataParallel(self.model)

        self.logger.info("Model initialized. Moving to GPU...")
        self.model.to(self.device)
        self.model.output_mask.binary_mask = self.model.output_mask.binary_mask.cpu()

        self.optimizers = dict()
        self.optimizer_types = dict()

        # Hypervolume
        hyper_volume_params = [{"params": list(self.model.hypervolume.parameters())}]

        self.optimizers["hypervolume"] = self.optim_types[
            self.configs.hypervolume_optimizer_type
        ](hyper_volume_params, lr=self.configs.lr)
        self.optimizer_types["hypervolume"] = self.configs.hypervolume_optimizer_type

        # Pose table
        if self.configs.epochs_sgd > 0:
            pose_table_params = [{"params": list(self.model.pose_table.parameters())}]
            self.optimizers["pose_table"] = self.optim_types[
                self.configs.pose_table_optimizer_type
            ](pose_table_params, lr=self.configs.lr_pose_table)
            self.optimizer_types["pose_table"] = self.configs.pose_table_optimizer_type

        # Z-latent-space conformations
        if self.configs.zdim > 0:
            if self.configs.use_conf_encoder:
                conf_encoder_params = [
                    {
                        "params": (
                            list(self.model.conf_cnn.parameters())
                            + list(self.model.conf_regressor.parameters())
                        )
                    }
                ]

                self.optimizers["conf_encoder"] = self.optim_types[
                    self.configs.conf_encoder_optimizer_type
                ](
                    conf_encoder_params,
                    lr=self.configs.lr_conf_encoder,
                    weight_decay=self.configs.wd,
                )
                self.optimizer_types[
                    "conf_encoder"
                ] = self.configs.conf_encoder_optimizer_type

            else:
                conf_table_params = [
                    {"params": list(self.model.conf_table.parameters())}
                ]

                self.optimizers["conf_table"] = self.optim_types[
                    self.configs.conf_table_optimizer_type
                ](conf_table_params, lr=self.configs.lr_conf_table)

                self.optimizer_types[
                    "conf_table"
                ] = self.configs.conf_table_optimizer_type

        self.optimized_modules = []

        # Complete initialization from a previous checkpoint
        if self.configs.load:
            checkpoint = torch.load(self.configs.load)

            for key in self.optimizers:
                self.optimizers[key].load_state_dict(
                    checkpoint["optimizers_state_dict"][key]
                )

        # Data loaders used to iterated over training batches of input images
        self.data_generator_pose_search = self.make_dataloader(
            batch_size=self.configs.batch_size_hps
        )
        self.data_generator = self.make_dataloader(
            batch_size=self.configs.batch_size_known_poses
        )
        self.data_generator_latent_optimization = self.make_dataloader(
            batch_size=self.configs.batch_size_sgd
        )

        # Save configurations within the output directory for future reference
        cfg_path = os.path.join(self.outdir, "config.yaml")
        data_norm_mean = float(self.data.norm[0])
        data_norm_std = float(self.data.norm[1])
        payload = {
            "training": dict(vars(self.configs)),
            "data_norm_mean": data_norm_mean,
            "data_norm_std": data_norm_std,
        }
        with open(cfg_path, "w") as f:
            yaml.safe_dump(payload, f)

        epsilon = 1e-8
        # Booleans used to track the current state of the training process
        self.log_latents = False
        self.pose_only = True
        self.pretraining = False
        self.is_in_pose_search_step = False
        self.use_point_estimates = False
        self.first_switch_to_point_estimates = True
        self.first_switch_to_point_estimates_conf = True

        if self.configs.load is not None:
            if self.start_epoch > self.epochs_pose_search:
                self.first_switch_to_point_estimates = False
            self.first_switch_to_point_estimates_conf = False

        self.use_kl_divergence = (
            not self.configs.zdim == 0
            and self.configs.variational_het
            and self.configs.beta_conf >= epsilon
        )
        self.use_trans_l1_regularizer = (
            self.configs.trans_l1_regularizer >= epsilon
            and not self.configs.use_gt_trans
            and not self.configs.no_trans
        )
        self.use_l2_smoothness_regularizer = (
            self.configs.l2_smoothness_regularizer >= epsilon
        )

        if self.configs.load:
            self.num_epochs = self.start_epoch + self.configs.epochs_sgd
            self.num_epochs += max(self.epochs_pose_search - self.start_epoch, 0)
        else:
            self.num_epochs = self.epochs_pose_search + self.configs.epochs_sgd

        self.n_particles_pretrain = (
            self.configs.n_imgs_pretrain
            if self.configs.n_imgs_pretrain >= 0
            else self.n_particles_dataset
        )

        # Placeholders for predicted latent variables, last input/output batch, losses
        self.in_dict_last = None
        self.y_pred_last = None

        self.predicted_rots = np.empty((self.n_tilts_dataset, 3, 3))
        self.predicted_trans = (
            np.empty((self.n_tilts_dataset, 2)) if not self.configs.no_trans else None
        )
        self.predicted_conf = (
            np.empty((self.n_particles_dataset, self.configs.zdim))
            if self.configs.zdim > 0
            else None
        )

        self.predicted_logvar = (
            np.empty((self.n_particles_dataset, self.configs.zdim))
            if self.configs.zdim > 0 and self.configs.variational_het
            else None
        )

        self.mask_particles_seen_at_last_epoch = np.zeros(self.n_particles_dataset)
        self.mask_tilts_seen_at_last_epoch = np.zeros(self.n_tilts_dataset)

        # Counters used to track the progress of the training process
        self.epoch = None
        self.run_times = {phase: [] for phase in self.run_phases}
        self.current_epoch_particles_count = 0
        self.total_batch_count = 0
        self.total_particles_count = 0
        self.batch_idx = 0
        self.cur_loss = None
        self.norm_mean, self.norm_std = self.data.norm

        # Activating Automatic Mixed Precision (AMP) model training through `torch.amp`
        self.scaler = None
        if self.configs.amp:
            if self.configs.pose_table_optimizer_type == "lbfgs":
                raise ValueError("AMP is not compatible with the lbfgs optimizer!")
            if (self.data.D - 1) % 8 != 0:
                self.logger.warning(
                    f"torch.amp mixed precision training is not optimized; "
                    f"image box size {self.data.D-1} is not a multiple of 8!"
                )
            if self.configs.zdim > 0 and self.configs.zdim % 8 != 0:
                self.logger.warning(
                    f"torch.amp mixed precision training is not optimized; "
                    f"{self.configs.zdim=} is not a multiple of 8!"
                )
            if self.configs.batch_size_hps % 8 != 0:
                self.logger.warning(
                    f"torch.amp mixed precision training is not optimized; "
                    f"{self.configs.batch_size_hps=} is not a multiple of 8!"
                )
            if self.configs.batch_size_sgd % 8 != 0:
                self.logger.warning(
                    f"torch.amp mixed precision training is not optimized; "
                    f"{self.configs.batch_size_sgd=} is not a multiple of 8!"
                )
            if self.configs.batch_size_known_poses % 8 != 0:
                self.logger.warning(
                    f"torch.amp mixed precision training is not optimized; "
                    f"{self.configs.batch_size_known_poses=} is not a multiple of 8!"
                )
            if self.configs.hypervolume_dim % 8 != 0:
                self.logger.warning(
                    f"torch.amp mixed precision training is not optimized; "
                    f"{self.configs.hypervolume_dim=} is not a multiple of 8!"
                )

            self.logger.info("Using Automatic Mixed Precision training via torch.amp")
            self.scaler = torch.amp.GradScaler()

    def train(self):
        self.logger.info("--- Training Starts Now ---")
        t_0 = dt.now()

        self.predicted_rots = (
            np.eye(3).reshape(1, 3, 3).repeat(self.n_tilts_dataset, axis=0)
        )
        self.predicted_trans = (
            np.zeros((self.n_tilts_dataset, 2)) if not self.configs.no_trans else None
        )
        self.predicted_conf = (
            np.zeros((self.n_particles_dataset, self.configs.zdim))
            if self.configs.zdim > 0
            else None
        )

        self.total_batch_count = 0
        self.total_particles_count = 0
        for epoch in range(self.start_epoch, self.num_epochs + 1):
            te = dt.now()

            self.epoch = epoch
            self.mask_particles_seen_at_last_epoch = np.zeros(self.n_particles_dataset)
            self.mask_tilts_seen_at_last_epoch = np.zeros(self.n_tilts_dataset)
            self.current_epoch_particles_count = 0
            self.optimized_modules = ["hypervolume"]

            self.pose_only = (
                self.total_particles_count < self.configs.pose_only_phase
                or self.configs.zdim == 0
                or epoch == 0
            )
            self.pretraining = self.epoch == 0
            self.is_in_pose_search_step = 0 < epoch <= self.epochs_pose_search
            self.use_point_estimates = epoch > self.epochs_pose_search

            n_max_particles = self.n_particles_dataset
            data_generator = self.data_generator

            # Pre-training
            if self.pretraining:
                n_max_particles = self.n_particles_pretrain
                self.logger.info(f"Will pretrain on {n_max_particles} particles")

            # HPS
            elif self.is_in_pose_search_step:
                n_max_particles = self.n_particles_dataset
                self.logger.info(f"Will use pose search on {n_max_particles} particles")
                data_generator = self.data_generator_pose_search

            # SGD
            elif self.use_point_estimates:
                if self.first_switch_to_point_estimates:
                    self.first_switch_to_point_estimates = False
                    self.logger.info("Switched to autodecoding poses")
                    self.logger.info(
                        "Initializing pose table from " "hierarchical pose search"
                    )
                    self.model.pose_table.initialize(
                        self.predicted_rots, self.predicted_trans
                    )
                    self.model.to(self.device)

                self.logger.info(
                    "Will use latent optimization on "
                    f"{self.n_particles_dataset} particles"
                )
                data_generator = self.data_generator_latent_optimization
                self.optimized_modules.append("pose_table")

            # GT poses
            else:
                raise RuntimeError("GT poses are not supported in this mode")

            # Z-latent-space conformations
            if not self.pose_only:
                if self.configs.use_conf_encoder:
                    self.optimized_modules.append("conf_encoder")

                else:
                    if self.first_switch_to_point_estimates_conf:
                        self.first_switch_to_point_estimates_conf = False

                        if self.configs.initial_conf is not None:
                            self.logger.info(
                                "Initializing conformation table " "from given z's"
                            )
                            self.model.conf_table.initialize(
                                utils.load_pkl(self.configs.initial_conf)
                            )

                        self.model.to(self.device)

                    self.optimized_modules.append("conf_table")

            will_make_summary = (
                (
                    self.configs.log_heavy_interval
                    and (epoch - 1) % self.configs.log_heavy_interval == 0
                )
                or self.is_in_pose_search_step
                or self.pretraining
                or epoch == self.num_epochs
            )
            self.log_latents = will_make_summary

            if will_make_summary:
                self.logger.info("Will make a full summary at the end of this epoch")

            for key in self.run_times.keys():
                self.run_times[key] = []

            end_time = time.time()
            self.cur_loss = 0

            # Inner loop
            for batch_idx, in_dict in enumerate(data_generator):
                self.batch_idx = batch_idx

                # with torch.autograd.detect_anomaly():
                self.train_step(in_dict, end_time=end_time)
                if self.configs.verbose_time:
                    torch.cuda.synchronize()

                end_time = time.time()

                if self.current_epoch_particles_count > n_max_particles:
                    break

            total_loss = self.cur_loss / self.current_epoch_particles_count
            self.logger.info(
                f"# =====> {self.epoch_type()} Epoch: {self.epoch} "
                f"finished in {dt.now() - te}; "
                f"total loss = {format(total_loss, '.6f')}"
            )

            # Image and pose summary at the end of each epoch
            if will_make_summary:
                self.save_latents()
                self.save_volume()
                self.save_model()

            # Update output mask -- epoch-based scaling
            if hasattr(self.output_mask, "update_epoch") and self.use_point_estimates:
                self.output_mask.update_epoch(self.configs.n_frequencies_per_epoch)

        t_total = dt.now() - t_0
        self.logger.info(
            f"Finished in {t_total} ({t_total / self.num_epochs} per epoch)"
        )

    def get_ctfs_at(self, index):
        batch_size = len(index)
        ctf_params_local = (
            self.ctf_params[index] if self.ctf_params is not None else None
        )

        if ctf_params_local is not None:
            freqs = self.lattice.freqs2d.unsqueeze(0).expand(
                batch_size, *self.lattice.freqs2d.shape
            ) / ctf_params_local[:, 0].view(batch_size, 1, 1)

            ctf_local = ctf.compute_ctf(
                freqs, *torch.split(ctf_params_local[:, 1:], 1, 1)
            ).view(batch_size, self.resolution, self.resolution)
        else:
            ctf_local = None

        return ctf_local

    def train_step(self, in_dict, end_time):
        if self.configs.verbose_time:
            torch.cuda.synchronize()
            self.run_times["dataloading"].append(time.time() - end_time)

        # Update output mask -- image-based scaling
        if hasattr(self.output_mask, "update") and self.is_in_pose_search_step:
            self.output_mask.update(self.total_particles_count)

        if self.is_in_pose_search_step:
            self.model.ps_params["l_min"] = self.configs.l_start

            if self.configs.output_mask == "circ":
                self.model.ps_params["l_max"] = self.configs.l_end
            else:
                self.model.ps_params["l_max"] = min(
                    self.output_mask.current_radius, self.configs.l_end
                )

        y_gt = in_dict["y"]
        ind = in_dict["index"]
        in_dict["tilt_index"] = in_dict["index"]
        ind_tilt = in_dict["tilt_index"]
        self.total_batch_count += 1
        batch_size = len(y_gt)
        self.total_particles_count += batch_size
        self.current_epoch_particles_count += batch_size

        # Move to GPU
        if self.configs.verbose_time:
            torch.cuda.synchronize()
        start_time_gpu = time.time()

        for key in in_dict.keys():
            if in_dict[key] is not None:
                in_dict[key] = in_dict[key].to(self.device)

        if self.configs.verbose_time:
            torch.cuda.synchronize()
            self.run_times["to_gpu"].append(time.time() - start_time_gpu)

        # Zero grad
        for key in self.optimized_modules:
            self.optimizers[key].zero_grad()

        # Forward pass
        if self.scaler is not None:
            with torch.amp.autocast("cuda"):
                latent_variables_dict, y_pred, y_gt_processed = self.forward_pass(
                    in_dict
                )

                if self.n_prcs > 1:
                    self.model.module.is_in_pose_search_step = False
                else:
                    self.model.is_in_pose_search_step = False
                # Loss
                if self.configs.verbose_time:
                    torch.cuda.synchronize()
                start_time_loss = time.time()
                total_loss, all_losses = self.loss(
                    y_pred, y_gt_processed, latent_variables_dict
                )
                if self.configs.verbose_time:
                    torch.cuda.synchronize()
                    self.run_times["loss"].append(time.time() - start_time_loss)
        else:
            latent_variables_dict, y_pred, y_gt_processed = self.forward_pass(in_dict)

            if self.n_prcs > 1:
                self.model.module.is_in_pose_search_step = False
            else:
                self.model.is_in_pose_search_step = False

            # Loss
            if self.configs.verbose_time:
                torch.cuda.synchronize()

            start_time_loss = time.time()
            total_loss, all_losses = self.loss(
                y_pred, y_gt_processed, latent_variables_dict
            )

            if self.configs.verbose_time:
                torch.cuda.synchronize()
                self.run_times["loss"].append(time.time() - start_time_loss)

        # Backward pass
        if self.configs.verbose_time:
            torch.cuda.synchronize()
        start_time_backward = time.time()
        if self.scaler is not None:
            self.scaler.scale(total_loss).backward()
        else:
            total_loss.backward()
        self.cur_loss += total_loss.item() * len(ind)

        for key in self.optimized_modules:
            if self.optimizer_types[key] == "adam":
                if self.scaler is not None:
                    self.scaler.step(self.optimizers[key])
                else:
                    self.optimizers[key].step()

            elif self.optimizer_types[key] == "lbfgs":

                def closure():
                    self.optimizers[key].zero_grad()
                    (
                        _latent_variables_dict,
                        _y_pred,
                        _y_gt_processed,
                    ) = self.forward_pass(in_dict)
                    _loss, _ = self.loss(
                        _y_pred, _y_gt_processed, _latent_variables_dict
                    )
                    _loss.backward()
                    return _loss.item()

                self.optimizers[key].step(closure)

            else:
                raise NotImplementedError

        if self.scaler is not None:
            self.scaler.update()

        if self.configs.verbose_time:
            torch.cuda.synchronize()

            self.run_times["backward"].append(time.time() - start_time_backward)

        # Detach from GPU
        if self.log_latents:
            self.in_dict_last = in_dict
            self.y_pred_last = y_pred

            if self.configs.verbose_time:
                torch.cuda.synchronize()

            start_time_cpu = time.time()
            rot_pred, trans_pred, conf_pred, logvar_pred = self.detach_latent_variables(
                latent_variables_dict
            )

            if self.configs.verbose_time:
                torch.cuda.synchronize()
                self.run_times["to_cpu"].append(time.time() - start_time_cpu)

            # Log
            if self.use_cuda:
                ind = ind.cpu()
                ind_tilt = ind_tilt.cpu()

            self.mask_particles_seen_at_last_epoch[ind] = 1
            self.mask_tilts_seen_at_last_epoch[ind_tilt] = 1
            self.predicted_rots[ind_tilt] = rot_pred.reshape(-1, 3, 3)

            if not self.configs.no_trans:
                self.predicted_trans[ind_tilt] = trans_pred.reshape(-1, 2)

            if self.configs.zdim > 0:
                self.predicted_conf[ind] = conf_pred

                if self.configs.variational_het:
                    self.predicted_logvar[ind] = logvar_pred

        else:
            self.run_times["to_cpu"].append(0.0)

        # Scalar summary
        if self.total_particles_count % self.configs.log_interval < batch_size:
            self.make_light_summary(all_losses)

    def detach_latent_variables(self, latent_variables_dict):
        rot_pred = latent_variables_dict["R"].detach().cpu().numpy()
        trans_pred = (
            latent_variables_dict["t"].detach().cpu().numpy()
            if not self.configs.no_trans
            else None
        )

        conf_pred = (
            latent_variables_dict["z"].detach().cpu().numpy()
            if self.configs.zdim > 0 and "z" in latent_variables_dict
            else None
        )

        logvar_pred = (
            latent_variables_dict["z_logvar"].detach().cpu().numpy()
            if self.configs.zdim > 0 and "z_logvar" in latent_variables_dict
            else None
        )

        return rot_pred, trans_pred, conf_pred, logvar_pred

    def forward_pass(self, in_dict):
        if self.configs.verbose_time:
            torch.cuda.synchronize()

        start_time_ctf = time.time()
        ctf_local = self.get_ctfs_at(in_dict["tilt_index"])

        if self.configs.verbose_time:
            torch.cuda.synchronize()
            self.run_times["ctf"].append(time.time() - start_time_ctf)

        # Forward pass
        if "hypervolume" in self.optimized_modules:
            self.model.hypervolume.train()
        else:
            self.model.hypervolume.eval()

        if hasattr(self.model, "conf_cnn"):
            if hasattr(self.model, "conf_regressor"):
                if "conf_encoder" in self.optimized_modules:
                    self.model.conf_cnn.train()
                    self.model.conf_regressor.train()
                else:
                    self.model.conf_cnn.eval()
                    self.model.conf_regressor.eval()

        if hasattr(self.model, "pose_table"):
            if "pose_table" in self.optimized_modules:
                self.model.pose_table.train()
            else:
                self.model.pose_table.eval()

        if hasattr(self.model, "conf_table"):
            if "conf_table" in self.optimized_modules:
                self.model.conf_table.train()
            else:
                self.model.conf_table.eval()

        in_dict["ctf"] = ctf_local
        if self.n_prcs > 1:
            self.model.module.pose_only = self.pose_only
            self.model.module.use_point_estimates = self.use_point_estimates
            self.model.module.pretrain = self.pretraining
            self.model.module.is_in_pose_search_step = self.is_in_pose_search_step
            self.model.module.use_point_estimates_conf = (
                not self.configs.use_conf_encoder
            )

        else:
            self.model.pose_only = self.pose_only
            self.model.use_point_estimates = self.use_point_estimates
            self.model.pretrain = self.pretraining
            self.model.is_in_pose_search_step = self.is_in_pose_search_step
            self.model.use_point_estimates_conf = not self.configs.use_conf_encoder

        out_dict = self.model(in_dict)

        self.run_times["encoder"].append(
            torch.mean(out_dict["time_encoder"].cpu())
            if self.configs.verbose_time
            else 0.0
        )
        self.run_times["decoder"].append(
            torch.mean(out_dict["time_decoder"].cpu())
            if self.configs.verbose_time
            else 0.0
        )
        self.run_times["decoder_coords"].append(
            torch.mean(out_dict["time_decoder_coords"].cpu())
            if self.configs.verbose_time
            else 0.0
        )
        self.run_times["decoder_query"].append(
            torch.mean(out_dict["time_decoder_query"].cpu())
            if self.configs.verbose_time
            else 0.0
        )

        latent_variables_dict = out_dict
        y_pred = out_dict["y_pred"]
        y_gt_processed = out_dict["y_gt_processed"]

        return latent_variables_dict, y_pred, y_gt_processed

    def loss(self, y_pred, y_gt, latent_variables_dict):
        """
        y_pred: [batch_size(, n_tilts), n_pts]
        y_gt: [batch_size(, n_tilts), n_pts]
        """
        all_losses = {}

        # Data loss
        data_loss = F.mse_loss(y_pred, y_gt)
        all_losses["Data Loss"] = data_loss.item()
        total_loss = data_loss

        # KL divergence
        if self.use_kl_divergence:
            kld_conf = kl_divergence_conf(latent_variables_dict)
            total_loss += self.configs.beta_conf * kld_conf / self.resolution**2
            all_losses["KL Div. Conf."] = kld_conf.item()

        # L1 regularization for translations
        if self.use_trans_l1_regularizer and self.use_point_estimates:
            trans_l1_loss = l1_regularizer(latent_variables_dict["t"])
            total_loss += self.configs.trans_l1_regularizer * trans_l1_loss
            all_losses["L1 Reg. Trans."] = trans_l1_loss.item()

        # L2 smoothness prior
        if self.use_l2_smoothness_regularizer:
            smoothness_loss = l2_frequency_bias(
                y_pred,
                self.lattice.freqs2d,
                self.output_mask.binary_mask,
                self.resolution,
            )
            total_loss += self.configs.l2_smoothness_regularizer * smoothness_loss
            all_losses["L2 Smoothness Loss"] = smoothness_loss.item()

        return total_loss, all_losses

    def make_light_summary(self, all_losses: dict[str, float]) -> None:
        """Creates a log describing progress within batches of a training epoch."""
        self.logger.info(
            f"# [Train Epoch: {self.epoch}/{self.num_epochs}] "
            f"[{self.current_epoch_particles_count}"
            f"/{self.n_particles_dataset} particles]"
        )

        if hasattr(self.output_mask, "current_radius"):
            all_losses["Mask Radius"] = self.output_mask.current_radius

        if self.model.trans_search_factor is not None:
            all_losses["Trans. Search Factor"] = self.model.trans_search_factor

        if self.configs.verbose_time:
            for key in self.run_times.keys():
                self.logger.info(
                    f"{key} time: {np.mean(np.array(self.run_times[key]))}"
                )

    def save_latents(self) -> None:
        """Write model's latent variables to file."""
        out_pose = os.path.join(self.outdir, f"pose.{self.epoch}.pkl")

        if self.configs.no_trans:
            with open(out_pose, "wb") as f:
                pickle.dump(self.predicted_rots, f)
        else:
            with open(out_pose, "wb") as f:
                pickle.dump((self.predicted_rots, self.predicted_trans), f)

        if self.configs.zdim > 0:
            out_conf = os.path.join(self.outdir, f"z.{self.epoch}.pkl")
            with open(out_conf, "wb") as f:
                pickle.dump(self.predicted_conf, f)

    def save_volume(self) -> None:
        """Write reconstructed volume to file."""
        out_mrc = os.path.join(self.outdir, f"reconstruct.{self.epoch}.mrc")

        self.model.hypervolume.eval()
        if hasattr(self.model, "conf_cnn"):
            if hasattr(self.model, "conf_regressor"):
                self.model.conf_cnn.eval()
                self.model.conf_regressor.eval()

        if hasattr(self.model, "pose_table"):
            self.model.pose_table.eval()
        if hasattr(self.model, "conf_table"):
            self.model.conf_table.eval()

        # For heterogeneous models reconstruct the volume at the latent coordinates
        # of the image whose embedding is closest to the mean of all embeddings
        if self.configs.zdim > 0:
            mean_z = np.mean(self.predicted_conf, axis=0)
            distances = np.linalg.norm(self.predicted_conf - mean_z, axis=1)
            closest_idx = np.argmin(distances)
            zval = self.predicted_conf[closest_idx].reshape(-1)
        else:
            zval = None

        vol = self.model.eval_volume(self.data.norm, zval=zval)
        write_mrc(out_mrc, vol.cpu().numpy().astype(np.float32))

    def save_model(self) -> None:
        """Write current PyTorch model state to file."""
        out_weights = os.path.join(self.outdir, f"weights.{self.epoch}.pkl")

        optimizers_state_dict = {}
        for key in self.optimizers.keys():
            optimizers_state_dict[key] = self.optimizers[key].state_dict()

        saved_objects = {
            "epoch": self.epoch,
            "model_state_dict": (
                self.model.module.state_dict()
                if self.n_prcs > 1
                else self.model.state_dict()
            ),
            "hypervolume_state_dict": (
                self.model.module.hypervolume.state_dict()
                if self.n_prcs > 1
                else self.model.hypervolume.state_dict()
            ),
            "hypervolume_params": self.model.hypervolume.get_building_params(),
            "optimizers_state_dict": optimizers_state_dict,
        }

        if hasattr(self.output_mask, "current_radius"):
            saved_objects["output_mask_radius"] = self.output_mask.current_radius

        torch.save(saved_objects, out_weights)

    def epoch_type(self) -> str:
        """Returns a label for the type of epoch currently being run."""
        if self.pretraining:
            return "Pretrain"
        elif self.is_in_pose_search_step:
            return "HPS"
        else:
            return "SGD"


def main(args: argparse.Namespace) -> None:
    # Build configs dict from args similar to TrainingConfigurations
    cfg = dict(
        particles=args.particles,
        ctf=args.ctf,
        datadir=args.datadir,
        ind=args.ind,
        relion31=args.relion31,
        invert_data=args.invert_data,
        load=args.load,
        lazy=args.lazy,
        max_threads=args.max_threads,
        log_interval=args.log_interval,
        log_heavy_interval=args.log_heavy_interval,
        verbose_time=args.verbose_time,
        shuffle=args.shuffle,
        num_workers=args.num_workers,
        shuffler_size=args.shuffler_size,
        multigpu=args.multigpu,
        amp=args.amp,
        batch_size_known_poses=args.batch_size_known_poses,
        batch_size_hps=args.batch_size_hps,
        batch_size_sgd=args.batch_size_sgd,
        hypervolume_optimizer_type=args.hypervolume_optimizer_type,
        pose_table_optimizer_type=args.pose_table_optimizer_type,
        conf_table_optimizer_type=args.conf_table_optimizer_type,
        conf_encoder_optimizer_type=args.conf_encoder_optimizer_type,
        lr=args.lr,
        lr_pose_table=args.lr_pose_table,
        lr_conf_table=args.lr_conf_table,
        lr_conf_encoder=args.lr_conf_encoder,
        wd=args.wd,
        n_imgs_pose_search=args.n_imgs_pose_search,
        epochs_sgd=args.epochs_sgd,
        pose_only_phase=args.pose_only_phase,
        output_mask=args.output_mask,
        add_one_frequency_every=args.add_one_frequency_every,
        n_frequencies_per_epoch=args.n_frequencies_per_epoch,
        max_freq=args.max_freq,
        window_radius_gt_real=args.window_radius_gt_real,
        beta_conf=args.beta_conf,
        trans_l1_regularizer=args.trans_l1_regularizer,
        l2_smoothness_regularizer=args.l2_smoothness_regularizer,
        variational_het=args.variational_het,
        zdim=args.zdim,
        std_z_init=args.std_z_init,
        use_conf_encoder=args.use_conf_encoder,
        depth_cnn=args.depth_cnn,
        channels_cnn=args.channels_cnn,
        kernel_size_cnn=args.kernel_size_cnn,
        resolution_encoder=args.resolution_encoder,
        explicit_volume=args.explicit_volume,
        hypervolume_layers=args.hypervolume_layers,
        hypervolume_dim=args.hypervolume_dim,
        pe_type=args.pe_type,
        pe_dim=args.pe_dim,
        feat_sigma=args.feat_sigma,
        hypervolume_domain=args.hypervolume_domain,
        pe_type_conf=args.pe_type_conf,
        n_imgs_pretrain=args.n_imgs_pretrain,
        l_start=args.l_start,
        l_end=args.l_end,
        n_iter=args.n_iter,
        t_extent=args.t_extent,
        t_n_grid=args.t_n_grid,
        t_x_shift=args.t_x_shift,
        t_y_shift=args.t_y_shift,
        no_trans_search_at_pose_search=args.no_trans_search_at_pose_search,
        n_kept_poses=args.n_kept_poses,
        base_healpy=args.base_healpy,
        no_trans=args.no_trans,
        seed=args.seed,
        norm=args.norm,
        initial_conf=args.initial_conf,
    )

    # Support --load latest (resolve from outdir or outdir/out)
    if cfg["load"] is not None:
        if cfg["load"].strip().lower() == "latest":
            weights_pkl, pose_pkl = utils.get_latest_checkpoint(cfg["outdir"])
            cfg["load"] = weights_pkl
            # Optionally expose pose file to initialize pose table later
            cfg["load_poses"] = pose_pkl if os.path.exists(pose_pkl) else None
        elif not os.path.exists(cfg["load"]):
            raise ValueError(
                f"Invalid load argument which must be a path to "
                f"a .pkl file or `latest`: {args.load}"
            )

    trainer = ModelTrainer(args.outdir, cfg)
    trainer.train()

    if args.do_analysis:
        anlz_cfgs = {
            "workdir": args.outdir,
            "epoch": trainer.epoch,
            "invert": cfg["invert_data"],
            "device": trainer.device,
            "skip_vol": False,
            "skip_umap": False,
            "pc": 2,
            "n_per_pc": 10,
            "ksample": 20,
            "apix": None,
            "flip": False,
            "downsample": None,
            "vol_start_index": 1,
        }
        cfg_file = os.path.join(args.outdir, "config.yaml")
        analyzer = ModelAnalyzer(args.outdir, anlz_cfgs, utils.load_yaml(cfg_file))
        analyzer.analyze()
