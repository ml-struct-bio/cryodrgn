"""Training engine for the cryoDRGN v4 (CryoDRGN-AI) reconstruction model.

This module contains the model training engine and the corresponding configuration
definitions for the aCryoDRGN-AI approach to particle reconstruction originally
introduced by Alex Levy in the drgnai package.

"""
import os
import pickle
import numpy as np
from dataclasses import dataclass
from typing import Any
import time

import torch
from torch import nn
import torch.nn.functional as F
from torch.utils.tensorboard import SummaryWriter

import cryodrgn.utils
import cryodrgn.ctf
from cryodrgn.mrcfile import write_mrc
from cryodrgn.dataset import make_dataloader
from cryodrgn.trainers import summary
from cryodrgn.models.losses import kl_divergence_conf, l1_regularizer, l2_frequency_bias
from cryodrgn.models.cryodrgnai import DRGNai, MyDataParallel
from cryodrgn.masking import CircularMask, FrequencyMarchingMask
from cryodrgn.trainers.reconstruction import (
    ReconstructionModelTrainer,
    ReconstructionModelConfigurations,
)


@dataclass
class SGDPoseSearchConfigurations(ReconstructionModelConfigurations):
    """The configurations used by the cryoDRGN v3 model training engine.

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

    > inherited from `ReconstructionModelConfigurations`:
        model       A label for the reconstruction algorithm to be used — must be either
                    `hps` for cryoDRGN v3 models or `amort` for cryoDRGN-AI models.
        zdim       The dimensionality of the latent space of conformations.
                    Thus zdim=0 for homogeneous models
                    and zdim>0 for hetergeneous models.
        num_epochs  The total number of epochs to use when training the model, not
                    including pretraining epoch(s).

        dataset     Label for the particle dataset to be used as input for the model.
                    If used, remaining input parameters can be omitted.
        particles   Path to the stack of particle images to use as input for the model.
                    Must be a (.mrcs/.txt/.star/.cs file).
        ctf         Path to the file storing contrast transfer function parameters
                    used to process the input particle images.
        poses       Path to the input particle poses data (.pkl).
        datadir     Path prefix to particle stack if loading relative paths from
                    a .star or .cs file.
        ind         Path to a numpy array saved as a .pkl used to filter
                    input particles.

        pose_estimation     Whether to perform ab-initio reconstruction ("abinit"),
                            reconstruction using fixed poses ("fixed"), or
                            reconstruction with SGD refinement of poses ("refine").
                            Default is to use fixed poses if poses file is given
                            and ab-initio otherwise.

        load:       Load model from given weights.<epoch>.pkl output file saved from
                    previous run of the engine.
                    Can also be given as "latest", in which the latest saved epoch
                    in the given output directory will be used.
        lazy        Whether to use lazy loading of data into memory in smaller batches.
                    Necessary if input dataset is too large to fit into memory.

        batch_size      The number of input images to use at a time when updating
                        the learning algorithm.

        multigpu        Whether to use all available GPUs available on this machine.
                        The default is to use only one GPU.

        log_interval    Print a log message every `N` number of training images.
        checkpoint      Save model results to file every `N` training epochs.

        pe_type     Label for the type of positional encoding to use.
        pe_dim      Number of frequencies to use in the positional encoding
                    (default: 64).
        volume_domain   Representation to use in the volume
                        decoder ("hartley" or "fourier").

    n_imgs_pose_search      The number of images the model needs to see in order to
                            complete the pose search portion of training.
    """

    # A parameter belongs to this configuration set if and only if it has a type and a
    # default value defined here, note that children classes inherit these parameters
    model: str = "cryodrgn-ai"

    # scheduling
    n_imgs_pose_search: int = 500000
    epochs_sgd: int = None
    pose_only_phase: int = 0
    use_gt_trans: bool = False
    invert_data: bool = False
    subtomo_averaging: bool = False
    # optimizers
    pose_table_optim_type: str = "adam"
    conf_table_optim_type: str = "adam"
    conf_encoder_optim_type: str = "adam"
    lr_conf_table: float = 1.0e-2
    lr_conf_encoder: float = 1.0e-4
    # masking
    output_mask: str = "circ"
    add_one_frequency_every: int = 100000
    n_frequencies_per_epoch: int = 10
    max_freq: int = None
    l_start_fm: int = 1
    # loss
    beta_conf: float = 0.0
    trans_l1_regularizer: float = 0.0
    l2_smoothness_regularizer: float = 0.0
    # conformations
    variational_het: bool = False
    std_z_init: float = 0.1
    use_conf_encoder: bool = False
    depth_cnn: int = 5
    channels_cnn: int = 32
    kernel_size_cnn: int = 3
    resolution_encoder: str = None
    initial_conf: str = None
    pe_type_conf: str = None
    # hypervolume
    volume_domain: str = "hartley"
    explicit_volume: bool = False
    # pre-training
    pretrain_with_gt_poses: bool = False
    # pose search
    n_iter: int = 4
    n_tilts_pose_search: int = 11
    average_over_tilts: bool = False
    no_trans_search_at_pose_search: bool = False
    n_kept_poses: int = 8
    # others
    palette_type: str = None

    # quick configs
    conf_estimation: str = None

    def __post_init__(self) -> None:
        super().__post_init__()

        if self.model != "cryodrgn-ai":
            raise ValueError(
                f"Mismatched model {self.model=}!=`cryodrgn-ai` "
                f"for {self.__class__.__name__}!"
            )
        if self.pose_estimation is not None:
            if self.pose_estimation == "refine":
                self.pose_learning_rate = 1.0e-4
        if self.conf_estimation is not None:
            if self.conf_estimation == "encoder":
                self.use_conf_encoder = True

        if self.explicit_volume and self.zdim >= 1:
            raise ValueError(
                "Explicit volumes do not support heterogeneous reconstruction."
            )
        if self.volume_optim_type not in {"adam"}:
            raise ValueError(
                f"Invalid value `{self.volume_optim_type=}` "
                f"for hypervolume optimizer type!"
            )
        if self.pose_table_optim_type not in {"adam", "lbfgs"}:
            raise ValueError(
                f"Invalid value `{self.pose_table_optim_type=}` "
                f"for pose table optimizer type!"
            )
        if self.conf_table_optim_type not in {"adam", "lbfgs"}:
            raise ValueError(
                f"Invalid value `{self.conf_table_optim_type=}` "
                f"for conformation table optimizer type!"
            )
        if self.conf_encoder_optim_type not in {"adam"}:
            raise ValueError(
                f"Invalid value `{self.conf_encoder_optim_type}` "
                "for conformation encoder optimizer type!"
            )

        if self.output_mask not in {"circ", "frequency_marching"}:
            raise ValueError(f"Invalid value {self.output_mask} for output_mask!")

        if self.pe_type not in {"gaussian"}:
            raise ValueError(f"Invalid value {self.pe_type} for pe_type!")

        if self.pe_type_conf not in {None, "geom"}:
            raise ValueError(f"Invalid value {self.pe_type_conf} for pe_type_conf!")

        if self.volume_domain not in {"hartley"}:
            raise ValueError(
                f"Invalid value {self.volume_domain} for hypervolume_domain."
            )

        if self.n_imgs_pose_search < 0:
            raise ValueError("n_imgs_pose_search must be greater than 0!")

        if self.use_conf_encoder and self.initial_conf:
            raise ValueError(
                "Conformations cannot be initialized when also using an encoder!"
            )

        if self.pose_estimation == "refine":
            self.n_imgs_pose_search = 0
            if self.poses is None:
                raise ValueError("Initial poses must be specified to be refined!")

        if self.subtomo_averaging:
            # TODO: Implement conformation encoder for subtomogram averaging.
            if self.use_conf_encoder:
                raise ValueError(
                    "Conformation encoder is not implemented "
                    "for subtomogram averaging!"
                )

            # TODO: Implement translation search for subtomogram averaging.
            if not (
                self.pose_estimation == "fixed"
                and (self.use_gt_trans or self.t_extent == 0.0)
            ):
                raise ValueError(
                    "Translation search is not implemented for subtomogram averaging!"
                )

        if self.average_over_tilts and self.n_tilts_pose_search % 2 == 0:
            raise ValueError(
                "`n_tilts_pose_search` must be odd to use `average_over_tilts`!"
            )

        if self.n_tilts_pose_search > self.n_tilts:
            raise ValueError("`n_tilts_pose_search` must be smaller than `n_tilts`!")

        if self.pose_estimation == "fixed":
            # "poses" include translations
            self.use_gt_trans = True
        if self.no_trans:
            self.t_extent = 0.0
        if self.t_extent == 0.0:
            self.t_ngrid = 1

    @property
    def file_dict(self) -> dict[str, Any]:
        """Organizing the parameter values for use in human-readable formats."""
        configs = super().file_dict

        # cnn
        cnn_params = {
            "use_conf_encoder": self.use_conf_encoder,
            "depth_cnn": self.depth_cnn,
            "channels_cnn": self.channels_cnn,
            "kernel_size_cnn": self.kernel_size_cnn,
        }
        # conformational encoder
        conf_regressor_params = {
            "std_z_init": self.std_z_init,
            "variational_het": self.variational_het,
        }
        # hypervolume
        hypervolume_params = {
            "explicit_volume": self.explicit_volume,
            "pe_type": self.pe_type,
            "pe_type_conf": self.pe_type_conf,
            "use_conf_encoder": self.use_conf_encoder,
            "lr_conf_table": self.lr_conf_table,
            "lr_conf_encoder": self.lr_conf_encoder,
            "beta_conf": self.beta_conf,
            "trans_l1_regularizer": self.trans_l1_regularizer,
            "l2_smoothness_regularizer": self.l2_smoothness_regularizer,
        }
        # pose search
        if self.pose_estimation != "fixed":
            ps_params = {
                "n_iter": self.n_iter,
                "n_kept_poses": self.n_kept_poses,
                "no_trans_search_at_pose_search": self.no_trans_search_at_pose_search,
                "n_tilts_pose_search": self.n_tilts_pose_search,
                "average_over_tilts": self.average_over_tilts,
            }
        else:
            ps_params = dict()

        configs["model_args"] = dict(
            hidden_layers=self.hidden_layers,
            hidden_dim=self.hidden_dim,
            **configs["model_args"],
            cnn_params=cnn_params,
            hypervolume_params=hypervolume_params,
            conf_regressor_params=conf_regressor_params,
        )
        if ps_params:
            configs["model_args"]["ps_params"] = ps_params
        configs["train_args"] = dict(
            **configs["train_args"],
            n_imgs_pose_search=self.n_imgs_pose_search,
            epochs_sgd=self.epochs_sgd,
        )

        for param_k in (
            cnn_params.keys()
            | conf_regressor_params.keys()
            | hypervolume_params.keys()
            | ps_params.keys()
        ):
            if param_k in configs["model_args"]:
                del configs["model_args"][param_k]

        return configs


class SGDPoseSearchTrainer(ReconstructionModelTrainer):
    """An engine for training the reconstruction model on particle data.

    Attributes
    ----------
    configs (TrainingConfigurations):   Values of all parameters that can be
                                        set by the user.

    particle_count (int):  The number of picked particles in the data.
    pretraining (bool):     Whether we are in the pretraining stage.
    epoch (int):    Which training epoch the model is in.

    logger (logging.Logger):    Utility for printing and writing information
                                about the model as it is running.
    """

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

    configs: SGDPoseSearchConfigurations
    config_cls = SGDPoseSearchConfigurations
    label = "cDRGN v4 training"

    def make_output_mask(self) -> CircularMask:
        if self.configs.output_mask == "circ":
            radius = self.configs.max_freq or self.lattice.D // 2
            output_mask = CircularMask(self.lattice, radius)

        elif self.configs.output_mask == "frequency_marching":
            output_mask = FrequencyMarchingMask(
                self.lattice,
                radius=self.configs.l_start_fm,
                radius_max=self.lattice.D // 2,
                add_one_every=self.configs.add_one_frequency_every,
            )

        else:
            raise NotImplementedError

        return output_mask

    def make_reconstruction_model(self, weights=None) -> nn.Module:
        output_mask = self.make_output_mask()

        if self.configs.zdim > 0:
            self.logger.info(
                "Heterogeneous reconstruction with " f"zdim = {self.configs.zdim}"
            )
        else:
            self.logger.info("Homogeneous reconstruction")

        cfgs = self.get_configs()
        model_args = cfgs["model_args"] if "model_args" in cfgs else cfgs

        if model_args["pe_dim"] is None:
            model_args["hypervolume_params"]["pe_dim"] = self.lattice.D // 2
        else:
            model_args["hypervolume_params"]["pe_dim"] = model_args["pe_dim"]

        model_args["hypervolume_params"].update(
            dict(
                hidden_dim=self.configs.hidden_dim,
                hidden_layers=self.configs.hidden_layers,
                pe_type=self.configs.pe_type,
                feat_sigma=self.configs.feat_sigma,
                volume_domain=self.configs.volume_domain,
            )
        )
        model_args["conf_regressor_params"].update(
            dict(
                zdim=self.configs.zdim,
            )
        )
        if "ps_params" in model_args:
            model_args["ps_params"].update(
                dict(
                    base_healpy=model_args["base_healpy"],
                    t_extent=model_args["t_extent"],
                    t_ngrid=model_args["t_ngrid"],
                    t_xshift=model_args["t_xshift"],
                    t_yshift=model_args["t_yshift"],
                )
            )

        model = DRGNai(
            self.lattice,
            output_mask,
            self.particle_count,
            self.image_count,
            model_args["cnn_params"],
            model_args["conf_regressor_params"],
            model_args["hypervolume_params"],
            resolution_encoder=self.configs.resolution_encoder,
            no_trans=self.configs.no_trans,
            use_gt_poses=self.configs.pose_estimation == "fixed",
            use_gt_trans=self.configs.use_gt_trans,
            will_use_point_estimates=self.epochs_sgd >= 1,
            ps_params=model_args["ps_params"] if "ps_params" in model_args else None,
            verbose_time=self.configs.verbose_time,
            pretrain_with_gt_poses=self.configs.pretrain_with_gt_poses,
            n_tilts_pose_search=self.configs.n_tilts_pose_search,
        )

        all_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
        self.logger.info(f"{all_params} parameters in model")

        return model

    @property
    def epochs_pose_search(self) -> int:
        if self.configs.n_imgs_pose_search > 0:
            epochs_pose_search = max(
                2, self.configs.n_imgs_pose_search // (self.particle_count + 1)
            )
        else:
            epochs_pose_search = 0

        return epochs_pose_search

    @property
    def epochs_sgd(self) -> int:
        if self.configs.epochs_sgd is None:
            epochs_sgd = self.configs.num_epochs - self.epochs_pose_search
        else:
            epochs_sgd = self.configs.epochs_sgd

        return epochs_sgd

    def __init__(self, configs: dict[str, Any], outdir: str) -> None:
        super().__init__(configs, outdir)
        self.configs: SGDPoseSearchConfigurations

        if self.configs.num_epochs is None:
            self.configs.num_epochs = self.epochs_sgd + self.epochs_pose_search

        # tensorboard writer
        self.summaries_dir = os.path.join(self.outdir, "summaries")
        os.makedirs(self.summaries_dir, exist_ok=True)
        self.writer = SummaryWriter(self.summaries_dir)
        self.logger.info("Will write tensorboard summaries " f"in {self.summaries_dir}")
        self.reconstruction_model.output_mask.binary_mask = (
            self.reconstruction_model.output_mask.binary_mask.cpu()
        )

        # TODO: Replace with DistributedDataParallel
        if self.configs.multigpu and torch.cuda.device_count() > 1:
            self.reconstruction_model = MyDataParallel(self.reconstruction_model)

        # pose table
        if self.configs.pose_estimation != "fixed":
            if self.epochs_sgd > 0:
                pose_table_params = [
                    {"params": list(self.reconstruction_model.pose_table.parameters())}
                ]
                self.optimizers["pose_table"] = self.optim_types[
                    self.configs.pose_table_optim_type
                ](pose_table_params, lr=self.configs.pose_learning_rate)
                self.optimizer_types["pose_table"] = self.configs.pose_table_optim_type

        # conformations
        if self.configs.zdim > 0:
            if self.configs.use_conf_encoder:
                conf_encoder_params = [
                    {
                        "params": (
                            list(self.reconstruction_model.conf_cnn.parameters())
                            + list(
                                self.reconstruction_model.conf_regressor.parameters()
                            )
                        )
                    }
                ]

                self.optimizers["conf_encoder"] = self.optim_types[
                    self.configs.conf_encoder_optim_type
                ](
                    conf_encoder_params,
                    lr=self.configs.lr_conf_encoder,
                    weight_decay=self.configs.weight_decay,
                )
                self.optimizer_types[
                    "conf_encoder"
                ] = self.configs.conf_encoder_optim_type

            else:
                conf_table_params = [
                    {"params": list(self.reconstruction_model.conf_table.parameters())}
                ]

                self.optimizers["conf_table"] = self.optim_types[
                    self.configs.conf_table_optim_type
                ](conf_table_params, lr=self.configs.lr_conf_table)
                self.optimizer_types["conf_table"] = self.configs.conf_table_optim_type

        self.optimized_modules = []
        # initialization from a previous checkpoint
        if self.configs.load:
            for key in self.optimizers:
                self.optimizers[key].load_state_dict(
                    self.checkpoint["optimizers_state_dict"][key]
                )

        # dataloaders
        self.data_generators = {"hps": None, "known": None, "sgd": None}

        self.data_generators["hps"] = make_dataloader(
            self.data,
            batch_size=self.configs.batch_size_hps,
            num_workers=self.configs.num_workers,
            shuffler_size=self.configs.shuffler_size,
        )

        self.data_generators["known"] = make_dataloader(
            self.data,
            batch_size=self.configs.batch_size_known_poses,
            num_workers=self.configs.num_workers,
            shuffler_size=self.configs.shuffler_size,
        )
        self.data_generators["sgd"] = make_dataloader(
            self.data,
            batch_size=self.configs.batch_size_sgd,
            num_workers=self.configs.num_workers,
            shuffler_size=self.configs.shuffler_size,
        )

        epsilon = 1e-8
        # booleans
        self.log_latents = False
        self.pose_only = True
        self.use_point_estimates = False
        self.first_switch_to_point_estimates = True
        self.first_switch_to_point_estimates_conf = True

        if self.configs.load is not None:
            if self.start_epoch >= self.epochs_pose_search:
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

        self.in_dict_last = None
        self.y_pred_last = None
        self.mask_particles_seen_at_last_epoch = np.zeros(self.particle_count)
        self.mask_tilts_seen_at_last_epoch = np.zeros(self.image_count)

        # counters
        self.run_times = {phase: [] for phase in self.run_phases}
        self.batch_idx = None
        self.cur_loss = None
        self.end_time = None

        self.predicted_logvar = (
            np.empty((self.particle_count, self.configs.zdim))
            if self.configs.zdim > 0 and self.configs.variational_het
            else None
        )

    def begin_epoch(self):
        self.configs: SGDPoseSearchConfigurations
        self.mask_particles_seen_at_last_epoch = np.zeros(self.particle_count)
        self.mask_tilts_seen_at_last_epoch = np.zeros(self.image_count)
        self.optimized_modules = ["hypervolume"]

        self.pose_only = self.configs.zdim == 0
        self.pose_only |= self.current_epoch == 0
        self.pose_only |= (
            self.epoch_images_seen is None
            or self.epoch_images_seen < self.configs.pose_only_phase
        )

        if self.configs.pose_estimation != "fixed":
            self.use_point_estimates = self.current_epoch >= max(
                0, self.epochs_pose_search
            )

        if self.current_epoch == 0:
            self.logger.info("Pretraining")
            self.data_iterator = self.data_generators["known"] or self.data_iterator

        # HPS
        elif self.in_pose_search_step:
            n_max_particles = self.particle_count
            self.logger.info(f"Will use pose search on {n_max_particles} particles")
            self.data_iterator = self.data_generators["hps"] or self.data_iterator

        # SGD
        elif self.use_point_estimates:
            if self.first_switch_to_point_estimates:
                self.first_switch_to_point_estimates = False
                self.logger.info("Switched to autodecoding poses")

                if self.configs.pose_estimation == "refine":
                    self.logger.info("Initializing pose table from ground truth")

                    poses_gt = cryodrgn.utils.load_pkl(self.configs.poses)
                    if poses_gt[0].ndim == 3:
                        # contains translations
                        rotmat_gt = torch.tensor(poses_gt[0]).float()
                        trans_gt = torch.tensor(poses_gt[1]).float()
                        trans_gt *= self.resolution

                        if self.ind is not None:
                            rotmat_gt = rotmat_gt[self.ind]
                            trans_gt = trans_gt[self.ind]

                    else:
                        rotmat_gt = torch.tensor(poses_gt).float()
                        trans_gt = None

                        if self.ind is not None:
                            rotmat_gt = rotmat_gt[self.ind]

                    self.reconstruction_model.pose_table.initialize(rotmat_gt, trans_gt)

                else:
                    self.logger.info(
                        "Initializing pose table from hierarchical pose search"
                    )
                    self.reconstruction_model.pose_table.initialize(
                        self.predicted_rots, self.predicted_trans
                    )

                self.reconstruction_model.to(self.device)

            self.logger.info(
                "Will use latent optimization on " f"{self.particle_count} particles"
            )

            self.data_iterator = self.data_generators["sgd"] or self.data_iterator
            self.optimized_modules.append("pose_table")

        # GT poses
        elif self.configs.pose_estimation == "fixed":
            self.data_iterator = self.data_generators["known"] or self.data_iterator

        # conformations
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
                        self.reconstruction_model.conf_table.initialize(
                            cryodrgn.utils.load_pkl(self.configs.initial_conf)
                        )

                    self.reconstruction_model.to(self.device)

                self.optimized_modules.append("conf_table")

        for key in self.run_times.keys():
            self.run_times[key] = []

        self.end_time = time.time()

    def end_epoch(self) -> None:
        # update output mask -- epoch-based scaling
        if hasattr(self.reconstruction_model.output_mask, "update_epoch"):
            if self.use_point_estimates:
                self.reconstruction_model.output_mask.update_epoch(
                    self.configs.n_frequencies_per_epoch
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

            ctf_local = cryodrgn.ctf.compute_ctf(
                freqs, *torch.split(ctf_params_local[:, 1:], 1, 1)
            ).view(batch_size, self.resolution, self.resolution)

        else:
            ctf_local = None

        return ctf_local

    def pretrain(self) -> None:
        self.begin_epoch()
        super().pretrain()

    def train_batch(self, batch: dict[str, torch.Tensor]) -> tuple:
        if self.configs.verbose_time:
            torch.cuda.synchronize()
            self.run_times["dataloading"].append(time.time() - self.end_time)

        # update output mask -- image-based scaling
        if hasattr(self.reconstruction_model.output_mask, "update"):
            if self.in_pose_search_step:
                self.reconstruction_model.output_mask.update(self.total_images_seen)

        if self.in_pose_search_step:
            self.reconstruction_model.ps_params["l_min"] = self.configs.l_start

            if self.configs.output_mask == "circ":
                self.reconstruction_model.ps_params["l_max"] = self.configs.l_end
            else:
                self.reconstruction_model.ps_params["l_max"] = min(
                    self.reconstruction_model.output_mask.current_radius,
                    self.configs.l_end,
                )

        if "tilt_indices" not in batch:
            batch["tilt_indices"] = batch["indices"]
        else:
            batch["tilt_indices"] = batch["tilt_indices"].reshape(-1)

        if self.configs.verbose_time:
            torch.cuda.synchronize()
        start_time_gpu = time.time()

        for key in batch.keys():
            batch[key] = batch[key].to(self.device)

        if self.configs.verbose_time:
            torch.cuda.synchronize()
            self.run_times["to_gpu"].append(time.time() - start_time_gpu)

        # zero grad
        for key in self.optimized_modules:
            self.optimizers[key].zero_grad()

        latent_variables_dict, y_pred, y_gt_processed = self.forward_pass(**batch)

        if isinstance(self.reconstruction_model, MyDataParallel):
            self.reconstruction_model.module.is_in_pose_search_step = False
        else:
            self.reconstruction_model.is_in_pose_search_step = False

        # loss
        if self.configs.verbose_time:
            torch.cuda.synchronize()

        start_time_loss = time.time()
        total_loss, all_losses = self.loss(
            y_pred, y_gt_processed, latent_variables_dict
        )

        if self.configs.verbose_time:
            torch.cuda.synchronize()
            self.run_times["loss"].append(time.time() - start_time_loss)

        # backward pass
        if self.configs.verbose_time:
            torch.cuda.synchronize()

        start_time_backward = time.time()
        total_loss.backward()

        for key in self.optimized_modules:
            if self.optimizer_types[key] == "adam":
                self.optimizers[key].step()

            elif self.optimizer_types[key] == "lbfgs":

                def closure():
                    self.optimizers[key].zero_grad()
                    (
                        _latent_variables_dict,
                        _y_pred,
                        _y_gt_processed,
                    ) = self.forward_pass(**batch)
                    _loss, _ = self.loss(
                        _y_pred, _y_gt_processed, _latent_variables_dict
                    )
                    _loss.backward()

                    return _loss

                self.optimizers[key].step(closure)

            else:
                raise NotImplementedError

        if self.configs.verbose_time:
            torch.cuda.synchronize()
            self.run_times["backward"].append(time.time() - start_time_backward)

        # detach
        rot_pred, trans_pred = None, None
        if self.will_make_checkpoint:
            self.in_dict_last = batch
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

            if self.use_cuda:
                batch["indices"] = batch["indices"].cpu()
                if batch["tilt_indices"] is not None:
                    batch["tilt_indices"] = batch["tilt_indices"].cpu()

            # keep track of predicted variables
            self.mask_particles_seen_at_last_epoch[batch["indices"]] = 1
            self.mask_tilts_seen_at_last_epoch[batch["tilt_indices"]] = 1
            if self.configs.pose_estimation != "fixed":
                self.predicted_rots[batch["tilt_indices"]] = rot_pred.reshape(-1, 3, 3)
            if not self.configs.no_trans and self.configs.pose_estimation != "fixed":
                self.predicted_trans[batch["tilt_indices"]] = trans_pred.reshape(-1, 2)

            if self.configs.zdim > 0 and self.predicted_conf is not None:
                self.predicted_conf[batch["indices"]] = conf_pred
                if self.configs.variational_het and self.predicted_logvar is not None:
                    self.predicted_logvar[batch["indices"]] = logvar_pred

        else:
            self.run_times["to_cpu"].append(0.0)

        batch_size = len(batch["indices"])
        if self.in_pretraining:
            self.accum_losses["total"] += total_loss.item() * batch_size

        all_losses["total"] = total_loss
        return (
            all_losses,
            batch["tilt_indices"],
            batch["indices"],
            rot_pred,
            trans_pred,
            latent_variables_dict["z"] if "z" in latent_variables_dict else None,
            None,
        )

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

    def forward_pass(self, y, y_real, tilt_indices, indices, rots=None, trans=None):
        if self.configs.verbose_time:
            torch.cuda.synchronize()

        start_time_ctf = time.time()
        if tilt_indices is not None:
            ctf_local = self.get_ctfs_at(tilt_indices)
        else:
            ctf_local = self.get_ctfs_at(indices)

        if self.configs.subtomo_averaging:
            ctf_local = ctf_local.reshape(-1, self.image_count, *ctf_local.shape[1:])

        if self.configs.verbose_time:
            torch.cuda.synchronize()
            self.run_times["ctf"].append(time.time() - start_time_ctf)

        # forward pass
        if "hypervolume" in self.optimized_modules:
            self.reconstruction_model.hypervolume.train()
        else:
            self.reconstruction_model.hypervolume.eval()

        if hasattr(self.reconstruction_model, "conf_cnn"):
            if hasattr(self.reconstruction_model, "conf_regressor"):
                if "conf_encoder" in self.optimized_modules:
                    self.reconstruction_model.conf_cnn.train()
                    self.reconstruction_model.conf_regressor.train()
                else:
                    self.reconstruction_model.conf_cnn.eval()
                    self.reconstruction_model.conf_regressor.eval()

        if hasattr(self.reconstruction_model, "pose_table"):
            if "pose_table" in self.optimized_modules:
                self.reconstruction_model.pose_table.train()
            else:
                self.reconstruction_model.pose_table.eval()

        if hasattr(self.reconstruction_model, "conf_table"):
            if "conf_table" in self.optimized_modules:
                self.reconstruction_model.conf_table.train()
            else:
                self.reconstruction_model.conf_table.eval()

        if isinstance(self.reconstruction_model, MyDataParallel):
            self.reconstruction_model.module.pose_only = self.pose_only
            self.reconstruction_model.module.use_point_estimates = (
                self.use_point_estimates
            )
            self.reconstruction_model.module.pretrain = self.in_pretraining
            self.reconstruction_model.module.is_in_pose_search_step = (
                self.in_pose_search_step
            )
            self.reconstruction_model.module.use_point_estimates_conf = (
                not self.configs.use_conf_encoder
            )

        else:
            self.reconstruction_model.pose_only = self.pose_only
            self.reconstruction_model.use_point_estimates = self.use_point_estimates
            self.reconstruction_model.pretrain = self.in_pretraining
            self.reconstruction_model.is_in_pose_search_step = self.in_pose_search_step
            self.reconstruction_model.use_point_estimates_conf = (
                not self.configs.use_conf_encoder
            )

        if self.configs.subtomo_averaging:
            tilt_indices = tilt_indices.reshape(y.shape[0:2])

        in_dict = {
            "y": y,
            "y_real": y_real,
            "indices": indices,
            "tilt_indices": tilt_indices,
            "ctf": ctf_local,
        }
        if rots is not None:
            in_dict["R"] = rots
        if trans is not None:
            in_dict["t"] = trans

        if in_dict["tilt_indices"] is None:
            in_dict["tilt_indices"] = in_dict["indices"]
        else:
            in_dict["tilt_indices"] = in_dict["tilt_indices"].reshape(-1)

        out_dict = self.reconstruction_model(in_dict)
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

        if self.configs.subtomo_averaging and self.configs.dose_exposure_correction:
            mask = self.reconstruction_model.output_mask.binary_mask
            a_pix = self.ctf_params[0, 0]

            dose_filters = self.data.get_dose_filters(
                in_dict["tilt_indices"].reshape(-1), self.lattice, a_pix
            ).reshape(*y_pred.shape[:2], -1)

            y_pred *= dose_filters[..., mask]

        return latent_variables_dict, y_pred, y_gt_processed

    def loss(self, y_pred, y_gt, latent_variables_dict):
        """
        y_pred: [batch_size(, n_tilts), n_pts]
        y: [batch_size(, n_tilts), n_pts]
        """
        all_losses = {}

        # data loss
        data_loss = F.mse_loss(y_pred, y_gt)
        all_losses["gen"] = data_loss.item()
        total_loss = data_loss

        # KL divergence
        if self.use_kl_divergence:
            kld_conf = kl_divergence_conf(latent_variables_dict)
            total_loss += self.configs.beta_conf * kld_conf / self.resolution**2
            all_losses["kld"] = kld_conf.item()

        # L1 regularization for translations
        if self.use_trans_l1_regularizer and self.use_point_estimates:
            trans_l1_loss = l1_regularizer(latent_variables_dict["t"])
            total_loss += self.configs.trans_l1_regularizer * trans_l1_loss
            all_losses["L1-Reg"] = trans_l1_loss.item()

        # L2 smoothness prior
        if self.use_l2_smoothness_regularizer:
            smoothness_loss = l2_frequency_bias(
                y_pred,
                self.lattice.freqs2d,
                self.reconstruction_model.output_mask.binary_mask,
                self.resolution,
            )
            total_loss += self.configs.l2_smoothness_regularizer * smoothness_loss
            all_losses["L2-Smooth"] = smoothness_loss.item()

        return total_loss, all_losses

    def get_configs(self) -> dict[str, Any]:
        configs = super().get_configs()

        configs["model_args"]["hypervolume_params"]["l_extent"] = self.lattice.extent
        if self.configs.pose_estimation != "fixed":
            tilt_fx = (
                self.data.get_tilting_func() if self.configs.subtomo_averaging else None
            )
            configs["model_args"]["ps_params"]["tilting_func"] = tilt_fx

        return configs

    def save_epoch_data(self):
        summary.make_img_summary(
            self.writer,
            self.in_dict_last,
            self.y_pred_last,
            self.reconstruction_model.output_mask,
            self.current_epoch,
        )

        # conformation
        if self.configs.zdim > 0:
            labels = None

            if self.configs.labels is not None:
                labels = cryodrgn.utils.load_pkl(self.configs.labels)

                if self.ind is not None:
                    labels = labels[self.ind]

            if self.mask_particles_seen_at_last_epoch is not None:
                mask_idx = self.mask_particles_seen_at_last_epoch > 0.5
            else:
                mask_idx = np.ones((self.particle_count,), dtype=bool)

            predicted_conf = self.predicted_conf[mask_idx]
            labels = labels[mask_idx] if labels is not None else None
            logvar = (
                self.predicted_logvar[mask_idx]
                if self.predicted_logvar is not None
                else None
            )
            summary.make_conf_summary(
                self.writer,
                predicted_conf,
                self.current_epoch,
                labels,
                pca=None,
                logvar=logvar,
                palette_type=self.configs.palette_type,
            )

        # pose
        rotmat_gt = None
        trans_gt = None
        shift = not self.configs.no_trans

        if self.mask_particles_seen_at_last_epoch is not None:
            mask_tilt_idx = self.mask_tilts_seen_at_last_epoch > 0.5
        else:
            mask_tilt_idx = np.ones((self.image_count,), dtype=bool)

        if self.configs.poses is not None:
            poses_gt = cryodrgn.utils.load_pkl(self.configs.poses)

            if poses_gt[0].ndim == 3:
                # contains translations
                rotmat_gt = torch.tensor(poses_gt[0]).float()
                trans_gt = torch.tensor(poses_gt[1]).float() * self.resolution

                if self.ind is not None:
                    rotmat_gt = rotmat_gt[self.ind]
                    trans_gt = trans_gt[self.ind]

            else:
                rotmat_gt = torch.tensor(poses_gt).float()
                trans_gt = None
                assert not shift, "Shift activated but trans not given in gt"

                if self.ind is not None:
                    rotmat_gt = rotmat_gt[self.ind]

            rotmat_gt = rotmat_gt[mask_tilt_idx]
            trans_gt = trans_gt[mask_tilt_idx] if trans_gt is not None else None

        if self.configs.pose_estimation != "fixed":
            predicted_rots = self.predicted_rots[mask_tilt_idx]

            if not self.configs.no_trans:
                predicted_trans = (
                    self.predicted_trans[mask_tilt_idx]
                    if self.predicted_trans is not None
                    else None
                )
            else:
                predicted_trans = None

            summary.make_pose_summary(
                self.writer,
                predicted_rots,
                predicted_trans,
                rotmat_gt,
                trans_gt,
                self.current_epoch,
                shift=shift,
            )

        self.save_latents()
        self.save_volume()
        self.save_model()

    def print_batch_summary(self, losses: dict[str, float]) -> None:
        kld_str = f"kld={losses['kld']:.4f}, " if "kld" in losses else ""
        self.logger.info(
            f"### Epoch [{self.current_epoch}/{self.configs.num_epochs}], "
            f"Batch [{self.epoch_batch_count}] "
            f"({self.epoch_images_seen}/{self.image_count} images); "
            f"gen loss={losses['gen']:.4g}, {kld_str}"
        )
        if self.reconstruction_model.trans_search_factor is not None:
            losses["Trans. Search Factor"] = getattr(
                self.reconstruction_model, "trans_search_factor"
            )

        summary.make_scalar_summary(self.writer, losses, self.total_images_seen)
        if self.configs.verbose_time:
            for key in self.run_times.keys():
                self.logger.info(
                    f"{key} time: {np.mean(np.array(self.run_times[key]))}"
                )

    def save_latents(self):
        """Write model's latent variables to file."""
        out_pose = os.path.join(self.outdir, f"pose.{self.current_epoch}.pkl")

        if self.configs.pose_estimation != "fixed":
            if self.configs.no_trans:
                with open(out_pose, "wb") as f:
                    pickle.dump(self.predicted_rots, f)
            else:
                if self.predicted_trans is not None:
                    out_trans = self.predicted_trans / self.model_resolution
                else:
                    out_trans = None
                with open(out_pose, "wb") as f:
                    pickle.dump((self.predicted_rots, out_trans), f)

        if self.configs.zdim > 0:
            out_conf = os.path.join(self.outdir, f"z.{self.current_epoch}.pkl")
            with open(out_conf, "wb") as f:
                pickle.dump(self.predicted_conf, f)

    def save_volume(self):
        """Write reconstructed volume to file."""
        out_mrc = os.path.join(self.outdir, f"reconstruct.{self.epoch_lbl}.mrc")
        self.reconstruction_model.hypervolume.eval()
        if hasattr(self.reconstruction_model, "conf_cnn"):
            if hasattr(self.reconstruction_model, "conf_regressor"):
                self.reconstruction_model.conf_cnn.eval()
                self.reconstruction_model.conf_regressor.eval()

        if hasattr(self.reconstruction_model, "pose_table"):
            self.reconstruction_model.pose_table.eval()
        if hasattr(self.reconstruction_model, "conf_table"):
            self.reconstruction_model.conf_table.eval()

        if self.configs.zdim > 0:
            zval = self.predicted_conf[0].reshape(-1)
        else:
            zval = None

        vol = -1.0 * self.reconstruction_model.eval_volume(
            norm=self.data.norm, zval=zval
        )
        write_mrc(out_mrc, np.array(vol, dtype=np.float32))

    def save_model(self):
        """Write model state to file."""
        out_weights = os.path.join(self.outdir, f"weights.{self.current_epoch}.pkl")
        optimizers_state_dict = {}
        for key in self.optimizers.keys():
            optimizers_state_dict[key] = self.optimizers[key].state_dict()

        saved_objects = {
            "epoch": self.current_epoch,
            "model_state_dict": (
                self.reconstruction_model.module.state_dict()
                if isinstance(self.reconstruction_model, MyDataParallel)
                else self.reconstruction_model.state_dict()
            ),
            "hypervolume_state_dict": (
                self.reconstruction_model.module.hypervolume.state_dict()
                if isinstance(self.reconstruction_model, MyDataParallel)
                else self.reconstruction_model.hypervolume.state_dict()
            ),
            "hypervolume_params": (
                self.reconstruction_model.hypervolume.get_building_params()
            ),
            "optimizers_state_dict": optimizers_state_dict,
        }

        if hasattr(self.reconstruction_model.output_mask, "current_radius"):
            saved_objects["output_mask_radius"] = getattr(
                self.reconstruction_model.output_mask, "current_radius"
            )

        torch.save(saved_objects, out_weights)

    @property
    def in_pose_search_step(self) -> bool:
        in_pose_search = False

        if self.configs.pose_estimation != "fixed":
            in_pose_search = 1 <= self.current_epoch <= self.epochs_pose_search

        return in_pose_search
