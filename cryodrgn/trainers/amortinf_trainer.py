import os
import pickle
from collections import OrderedDict
import numpy as np
from typing import Any
import time

import torch
from torch import nn
import torch.nn.functional as F
from torch.utils.tensorboard import SummaryWriter

import cryodrgn.utils
from cryodrgn import ctf, mrc
from cryodrgn.dataset import make_dataloader
from cryodrgn.trainers import summary
from cryodrgn.models.losses import kl_divergence_conf, l1_regularizer, l2_frequency_bias
from cryodrgn.models.amortized_inference import DRGNai, MyDataParallel
from cryodrgn.masking import CircularMask, FrequencyMarchingMask
from cryodrgn.trainers._base import ModelTrainer, ModelConfigurations


class AmortizedInferenceConfigurations(ModelConfigurations):

    __slots__ = (
        "batch_size_known_poses",
        "batch_size_hps",
        "batch_size_sgd",
        "pose_table_optimizer_type",
        "conf_table_optimizer_type",
        "conf_encoder_optimizer_type",
        "lr_pose_table",
        "lr_conf_table",
        "lr_conf_encoder",
        "n_imgs_pose_search",
        "epochs_sgd",
        "pose_only_phase",
        "output_mask",
        "add_one_frequency_every",
        "n_frequencies_per_epoch",
        "max_freq",
        "l_start_fm",
        "beta_conf",
        "trans_l1_regularizer",
        "l2_smoothness_regularizer",
        # conformations
        "variational_het",
        "std_z_init",
        "use_conf_encoder",
        "depth_cnn",
        "channels_cnn",
        "kernel_size_cnn",
        "resolution_encoder",
        "initial_conf",
        # hypervolume
        "explicit_volume",
        "pe_type_conf",
        # pre-training
        "pretrain_with_gt_poses",
        # pose search
        "n_iter",
        "n_tilts_pose_search",
        "average_over_tilts",
        "no_trans_search_at_pose_search",
        "n_kept_poses",
        # subtomogram averaging
        "palette_type",
    )
    default_values = OrderedDict(
        {
            # data loading
            "batch_size_known_poses": 32,
            "batch_size_hps": 8,
            "batch_size_sgd": 256,
            # optimizers
            "pose_table_optimizer_type": "adam",
            "conf_table_optimizer_type": "adam",
            "conf_encoder_optimizer_type": "adam",
            "lr_pose_table": 1.0e-3,
            "lr_conf_table": 1.0e-2,
            "lr_conf_encoder": 1.0e-4,
            # scheduling
            "n_imgs_pose_search": 500000,
            "epochs_sgd": 100,
            "pose_only_phase": 0,
            # masking
            "output_mask": "circ",
            "add_one_frequency_every": 100000,
            "n_frequencies_per_epoch": 10,
            "max_freq": None,
            "l_start_fm": 12,
            # loss
            "beta_conf": 0.0,
            "trans_l1_regularizer": 0.0,
            "l2_smoothness_regularizer": 0.0,
            # conformations
            "variational_het": False,
            "std_z_init": 0.1,
            "use_conf_encoder": False,
            "depth_cnn": 5,
            "channels_cnn": 32,
            "kernel_size_cnn": 3,
            "resolution_encoder": None,
            "initial_conf": None,
            # hypervolume
            "explicit_volume": False,
            "pe_type_conf": None,
            # pre-training
            "pretrain_with_gt_poses": False,
            # pose search
            "n_iter": 4,
            "n_tilts_pose_search": 11,
            "average_over_tilts": False,
            "no_trans_search_at_pose_search": False,
            "n_kept_poses": 8,
            # others
            "palette_type": None,
        }
    )

    quick_config = OrderedDict(
        {
            "capture_setup": {
                "spa": dict(),
                "et": {
                    "subtomo_averaging": True,
                    "shuffler_size": 0,
                    "num_workers": 0,
                    "t_extent": 0.0,
                    "batch_size_known_poses": 8,
                    "batch_size_sgd": 32,
                    "n_imgs_pose_search": 150000,
                    "pose_only_phase": 50000,
                    "lr_pose_table": 1.0e-5,
                },
            },
            "reconstruction_type": {"homo": {"z_dim": 0}, "het": dict()},
            "pose_estimation": {
                "abinit": dict(),
                "refine": {"refine_gt_poses": True, "lr_pose_table": 1.0e-4},
                "fixed": {"use_gt_poses": True},
            },
            "conf_estimation": {
                None: dict(),
                "autodecoder": dict(),
                "refine": dict(),
                "encoder": {"use_conf_encoder": True},
            },
        }
    )

    def __init__(self, config_vals: dict[str, Any]) -> None:
        super().__init__(config_vals)

        if "model" in config_vals and config_vals["model"] != "amort":
            raise ValueError(
                f"Mismatched model {config_vals['model']} "
                "for AmortizedInferenceTrainer!"
            )

        if "explicit_volume" in config_vals:
            if config_vals["explicit_volume"] and config_vals["z_dim"] >= 1:
                raise ValueError(
                    "Explicit volumes do not support heterogeneous reconstruction."
                )

        if "dataset" in config_vals:
            if config_vals["dataset"] is None:
                if config_vals["particles"] is None:
                    raise ValueError(
                        "As dataset was not specified, please " "specify particles!"
                    )
                if config_vals["ctf"] is None:
                    raise ValueError("As dataset wasn't specified, please specify ctf!")

        if "hypervolume_optimizer_type" in config_vals:
            if config_vals["hypervolume_optimizer_type"] not in {"adam"}:
                raise ValueError(
                    "Invalid value "
                    f"`{config_vals['hypervolume_optimizer_type']}` "
                    "for hypervolume_optimizer_type!"
                )

        if "pose_table_optimizer_type" in config_vals:
            if config_vals["pose_table_optimizer_type"] not in {"adam", "lbfgs"}:
                raise ValueError(
                    "Invalid value "
                    f"`{config_vals['pose_table_optimizer_type']}` "
                    "for pose_table_optimizer_type!"
                )

        if "conf_table_optimizer_type" in config_vals:
            if config_vals["conf_table_optimizer_type"] not in {"adam", "lbfgs"}:
                raise ValueError(
                    "Invalid value "
                    f"`{config_vals['conf_table_optimizer_type']}` "
                    "for conf_table_optimizer_type!"
                )

        if "conf_encoder_optimizer_type" in config_vals:
            if config_vals["conf_encoder_optimizer_type"] not in {"adam"}:
                raise ValueError(
                    "Invalid value "
                    f"`{config_vals['conf_encoder_optimizer_type']}` "
                    "for conf_encoder_optimizer_type!"
                )

        if "output_mask" in config_vals:
            if config_vals["output_mask"] not in {"circ", "frequency_marching"}:
                raise ValueError(
                    f"Invalid value {config_vals['output_mask']} for output_mask!"
                )

        if "pe_type" in config_vals and config_vals["pe_type"] not in {"gaussian"}:
            raise ValueError(f"Invalid value {config_vals['pe_type']} for pe_type!")

        if "pe_type_conf" in config_vals:
            if config_vals["pe_type_conf"] not in {None, "geom"}:
                raise ValueError(
                    f"Invalid value {config_vals['pe_type_conf']} for pe_type_conf!"
                )

        if "hypervolume_domain" in config_vals:
            if config_vals["hypervolume_domain"] not in {"hartley"}:
                raise ValueError(
                    f"Invalid value {config_vals['hypervolume_domain']} "
                    "for hypervolume_domain."
                )

        if self.volume_domain is None:
            self.volume_domain = "hartley"

        if "n_imgs_pose_search" in config_vals:
            if config_vals["n_imgs_pose_search"] < 0:
                raise ValueError("n_imgs_pose_search must be greater than 0!")

        if "use_conf_encoder" in config_vals and "initial_conf" in config_vals:
            if config_vals["use_conf_encoder"] and config_vals["initial_conf"]:
                raise ValueError(
                    "Conformations cannot be initialized when also using an encoder!"
                )

        if "use_gt_trans" in config_vals and "pose" in config_vals:
            if config_vals["use_gt_trans"] and config_vals["pose"] is None:
                raise ValueError(
                    "Poses must be specified to use ground-truth translations!"
                )

        if "refine_gt_poses" in config_vals and config_vals["refine_gt_poses"]:
            config_vals["n_imgs_pose_search"] = 0
            if "pose" not in config_vals or config_vals["pose"] is None:
                raise ValueError("Initial poses must be specified to be refined!")

        if "subtomo_averaging" in config_vals:
            if config_vals["subtomo_averaging"]:
                # TODO: Implement conformation encoder for subtomogram averaging.
                if "use_conf_encoder" in config_vals:
                    if config_vals["use_conf_encoder"]:
                        raise ValueError(
                            "Conformation encoder is not implemented "
                            "for subtomogram averaging!"
                        )

                # TODO: Implement translation search for subtomogram averaging.
                if not (
                    "use_gt_poses" in config_vals
                    and config_vals["use_gt_poses"]
                    or "use_gt_trans" in config_vals
                    and config_vals["use_gt_trans"]
                    or "t_extent" in config_vals
                    and config_vals["t_extent"] == 0.0
                ):
                    raise ValueError(
                        "Translation search is not implemented "
                        "for subtomogram averaging!"
                    )

                if (
                    "average_over_tilts" in config_vals
                    and config_vals["average_over_tilts"]
                    and "n_tilts_pose_search" in config_vals
                    and config_vals["n_tilts_pose_search"] % 2 == 0
                ):
                    raise ValueError(
                        "`n_tilts_pose_search` must be odd to use `average_over_tilts`!"
                    )

                if "n_tilts_pose_search" in config_vals and "n_tilts" in config_vals:
                    if config_vals["n_tilts_pose_search"] > config_vals["n_tilts"]:
                        raise ValueError(
                            "`n_tilts_pose_search` must be smaller than `n_tilts`!"
                        )

        if "use_gt_poses" in config_vals and config_vals["use_gt_poses"]:
            # "poses" include translations
            config_vals["use_gt_trans"] = True
            if "pose" in config_vals and config_vals["pose"] is None:
                raise ValueError("Ground truth poses must be specified!")

        if "no_trans" in config_vals and config_vals["no_trans"]:
            config_vals["t_extent"] = 0.0
        if "t_extent" in config_vals and config_vals["t_extent"] == 0.0:
            config_vals["t_n_grid"] = 1


class AmortizedInferenceTrainer(ModelTrainer):
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

    config_cls = AmortizedInferenceConfigurations

    def make_volume_model(self) -> nn.Module:
        self.configs: AmortizedInferenceConfigurations

        # output mask
        if self.configs.output_mask == "circ":
            radius = self.configs.max_freq or self.lattice.D // 2
            output_mask = CircularMask(self.lattice, radius)

        elif self.configs.output_mask == "frequency_marching":
            output_mask = FrequencyMarchingMask(
                self.lattice,
                self.lattice.D // 2,
                radius=self.configs.l_start_fm,
                add_one_every=self.configs.add_one_frequency_every,
            )

        else:
            raise NotImplementedError

        # cnn
        cnn_params = {
            "conf": self.configs.use_conf_encoder,
            "depth_cnn": self.configs.depth_cnn,
            "channels_cnn": self.configs.channels_cnn,
            "kernel_size_cnn": self.configs.kernel_size_cnn,
        }

        # conformational encoder
        if self.configs.z_dim > 0:
            self.logger.info(
                "Heterogeneous reconstruction with " f"z_dim = {self.configs.z_dim}"
            )
        else:
            self.logger.info("Homogeneous reconstruction")

        conf_regressor_params = {
            "z_dim": self.configs.z_dim,
            "std_z_init": self.configs.std_z_init,
            "variational": self.configs.variational_het,
        }

        # hypervolume
        hyper_volume_params = {
            "explicit_volume": self.configs.explicit_volume,
            "n_layers": self.configs.hidden_layers,
            "hidden_dim": self.configs.hidden_dim,
            "pe_type": self.configs.pe_type,
            "pe_dim": self.configs.pe_dim,
            "feat_sigma": self.configs.feat_sigma,
            "domain": self.configs.volume_domain,
            "extent": self.lattice.extent,
            "pe_type_conf": self.configs.pe_type_conf,
        }

        # pose search
        if self.epochs_pose_search > 0:
            ps_params = {
                "l_min": self.configs.l_start,
                "l_max": self.configs.l_end,
                "t_extent": self.configs.t_extent,
                "t_n_grid": self.configs.t_ngrid,
                "niter": self.configs.n_iter,
                "nkeptposes": self.configs.n_kept_poses,
                "base_healpy": self.configs.base_healpy,
                "t_xshift": self.configs.t_xshift,
                "t_yshift": self.configs.t_yshift,
                "no_trans_search_at_pose_search": self.configs.no_trans_search_at_pose_search,
                "n_tilts_pose_search": self.configs.n_tilts_pose_search,
                "tilting_func": (
                    self.data.get_tilting_func()
                    if self.configs.subtomo_averaging
                    else None
                ),
                "average_over_tilts": self.configs.average_over_tilts,
            }
        else:
            ps_params = None

        return DRGNai(
            self.lattice,
            output_mask,
            self.particle_count,
            self.image_count,
            cnn_params,
            conf_regressor_params,
            hyper_volume_params,
            resolution_encoder=self.configs.resolution_encoder,
            no_trans=self.configs.no_trans,
            use_gt_poses=self.configs.use_gt_poses,
            use_gt_trans=self.configs.use_gt_trans,
            will_use_point_estimates=self.configs.epochs_sgd >= 1,
            ps_params=ps_params,
            verbose_time=self.configs.verbose_time,
            pretrain_with_gt_poses=self.configs.pretrain_with_gt_poses,
            n_tilts_pose_search=self.configs.n_tilts_pose_search,
        )

    @property
    def epochs_pose_search(self) -> int:
        if self.configs.n_imgs_pose_search > 0:
            epochs_pose_search = max(
                2, self.configs.n_imgs_pose_search // self.particle_count + 1
            )
        else:
            epochs_pose_search = 0

        return epochs_pose_search

    def __init__(self, configs: dict[str, Any]) -> None:
        super().__init__(configs)
        self.configs: AmortizedInferenceConfigurations
        self.model = self.volume_model

        self.batch_size_known_poses = self.configs.batch_size_known_poses * self.n_prcs
        self.batch_size_hps = self.configs.batch_size_hps * self.n_prcs
        self.batch_size_sgd = self.configs.batch_size_sgd * self.n_prcs

        # tensorboard writer
        self.summaries_dir = os.path.join(self.configs.outdir, "summaries")
        os.makedirs(self.summaries_dir, exist_ok=True)
        self.writer = SummaryWriter(self.summaries_dir)
        self.logger.info("Will write tensorboard summaries " f"in {self.summaries_dir}")

        # TODO: Replace with DistributedDataParallel
        if self.n_prcs > 1:
            self.model = MyDataParallel(self.volume_model)

        self.model.output_mask.binary_mask = self.model.output_mask.binary_mask.cpu()
        self.optimizers = {"hypervolume": self.volume_optimizer}
        self.optimizer_types = {"hypervolume": self.configs.volume_optim_type}

        # pose table
        if not self.configs.use_gt_poses:
            if self.configs.epochs_sgd > 0:
                pose_table_params = [
                    {"params": list(self.model.pose_table.parameters())}
                ]

                self.optimizers["pose_table"] = self.optim_types[
                    self.configs.pose_table_optimizer_type
                ](pose_table_params, lr=self.configs.lr_pose_table)
                self.optimizer_types[
                    "pose_table"
                ] = self.configs.pose_table_optimizer_type

        # conformations
        if self.configs.z_dim > 0:
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
                    weight_decay=self.configs.weight_decay,
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
        self.data_generators = {"hps": None, "known": None, "sgd": None}

        # dataloaders
        if self.batch_size_hps != self.configs.batch_size:
            self.data_generators["hps"] = make_dataloader(
                self.data,
                batch_size=self.batch_size_hps,
                num_workers=self.configs.num_workers,
                shuffler_size=self.configs.shuffler_size,
            )

        if self.batch_size_known_poses != self.configs.batch_size:
            self.data_generators["known"] = make_dataloader(
                self.data,
                batch_size=self.batch_size_known_poses,
                num_workers=self.configs.num_workers,
                shuffler_size=self.configs.shuffler_size,
            )
        if self.batch_size_sgd != self.configs.batch_size:
            self.data_generators["sgd"] = make_dataloader(
                self.data,
                batch_size=self.batch_size_sgd,
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
            not self.configs.z_dim == 0
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

        self.num_epochs = self.epochs_pose_search + self.configs.epochs_sgd
        if self.configs.load:
            self.num_epochs += self.start_epoch

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
            np.empty((self.particle_count, self.configs.z_dim))
            if self.configs.z_dim > 0 and self.configs.variational_het
            else None
        )

    def begin_epoch(self):
        self.configs: AmortizedInferenceConfigurations
        self.mask_particles_seen_at_last_epoch = np.zeros(self.particle_count)
        self.mask_tilts_seen_at_last_epoch = np.zeros(self.image_count)
        self.optimized_modules = ["hypervolume"]

        self.pose_only = (
                self.total_images_seen < self.configs.pose_only_phase
                or self.configs.z_dim == 0
        )

        if not self.configs.use_gt_poses:
            self.use_point_estimates = self.current_epoch >= max(
                0, self.epochs_pose_search
            )

        # HPS
        if self.is_in_pose_search_step:
            n_max_particles = self.particle_count
            self.logger.info(f"Will use pose search on {n_max_particles} particles")
            self.data_iterator = self.data_generators["hps"] or self.data_iterator

        # SGD
        elif self.use_point_estimates:
            if self.first_switch_to_point_estimates:
                self.first_switch_to_point_estimates = False
                self.logger.info("Switched to autodecoding poses")

                if self.configs.refine_gt_poses:
                    self.logger.info("Initializing pose table from ground truth")

                    poses_gt = cryodrgn.utils.load_pkl(self.configs.pose)
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

                    self.model.pose_table.initialize(rotmat_gt, trans_gt)

                else:
                    self.logger.info(
                        "Initializing pose table from hierarchical pose search"
                    )
                    self.model.pose_table.initialize(
                        self.predicted_rots, self.predicted_trans
                    )

                self.model.to(self.device)

            self.logger.info(
                "Will use latent optimization on " f"{self.particle_count} particles"
            )

            self.data_iterator = self.data_generators["sgd"] or self.data_iterator
            self.optimized_modules.append("pose_table")

        # GT poses
        else:
            assert self.configs.use_gt_poses
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
                        self.model.conf_table.initialize(
                            cryodrgn.utils.load_pkl(self.configs.initial_conf)
                        )

                    self.model.to(self.device)

                self.optimized_modules.append("conf_table")

        for key in self.run_times.keys():
            self.run_times[key] = []

        self.end_time = time.time()

    def end_epoch(self) -> None:
        # update output mask -- epoch-based scaling
        if hasattr(self.model.output_mask, "update_epoch") and self.use_point_estimates:
            self.model.output_mask.update_epoch(self.configs.n_frequencies_per_epoch)

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

    def train_batch(self, batch: tuple) -> None:
        y_gt, tilt_ind, ind = batch

        if self.configs.verbose_time:
            torch.cuda.synchronize()
            self.run_times["dataloading"].append(time.time() - self.end_time)

        # update output mask -- image-based scaling
        if hasattr(self.model.output_mask, "update") and self.is_in_pose_search_step:
            self.model.output_mask.update(self.total_images_seen)

        if self.is_in_pose_search_step:
            self.model.ps_params["l_min"] = self.configs.l_start

            if self.configs.output_mask == "circ":
                self.model.ps_params["l_max"] = self.configs.l_end
            else:
                self.model.ps_params["l_max"] = min(
                    self.model.output_mask.current_radius, self.configs.l_end
                )

        if tilt_ind is not None:
            ind_tilt = tilt_ind.reshape(-1)
            ind_tilt = ind_tilt.to(self.device)
        else:
            ind_tilt = None

        # move to gpu
        if self.configs.verbose_time:
            torch.cuda.synchronize()
        start_time_gpu = time.time()

        y_gt = y_gt.to(self.device)
        ind = ind.to(self.device)
        if self.configs.verbose_time:
            torch.cuda.synchronize()
            self.run_times["to_gpu"].append(time.time() - start_time_gpu)

        # zero grad
        for key in self.optimized_modules:
            self.optimizers[key].zero_grad()

        in_dict = {
            "y": y_gt,
            "index": ind,
            "tilt_index": tilt_ind,
        }

        # forward pass
        latent_variables_dict, y_pred, y_gt_processed = self.forward_pass(
            y_gt, tilt_ind, ind
        )

        if self.n_prcs > 1:
            self.model.module.is_in_pose_search_step = False
        else:
            self.model.is_in_pose_search_step = False

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
        # self.cur_loss += total_loss.item() * len(ind)

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
                    ) = self.forward_pass(in_dict)
                    _loss, _ = self.loss(
                        _y_pred, _y_gt_processed, _latent_variables_dict
                    )
                    _loss.backward()
                    return _loss.item()

                self.optimizers[key].step(closure)

            else:
                raise NotImplementedError

        if self.configs.verbose_time:
            torch.cuda.synchronize()
            self.run_times["backward"].append(time.time() - start_time_backward)

        # detach
        if self.will_make_checkpoint:
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

            # log
            if self.use_cuda:
                ind = ind.cpu()
                if ind_tilt is not None:
                    ind_tilt = ind_tilt.cpu()

            self.mask_particles_seen_at_last_epoch[ind] = 1
            self.mask_tilts_seen_at_last_epoch[ind_tilt] = 1
            self.predicted_rots[ind_tilt] = rot_pred.reshape(-1, 3, 3)

            if not self.configs.no_trans:
                self.predicted_trans[ind_tilt] = trans_pred.reshape(-1, 2)

            if self.configs.z_dim > 0:
                self.predicted_conf[ind] = conf_pred

                if self.configs.variational_het:
                    self.predicted_logvar[ind] = logvar_pred

        else:
            self.run_times["to_cpu"].append(0.0)

        self.end_time = time.time()

    def detach_latent_variables(self, latent_variables_dict):
        rot_pred = latent_variables_dict["R"].detach().cpu().numpy()
        trans_pred = (
            latent_variables_dict["t"].detach().cpu().numpy()
            if not self.configs.no_trans
            else None
        )

        conf_pred = (
            latent_variables_dict["z"].detach().cpu().numpy()
            if self.configs.z_dim > 0 and "z" in latent_variables_dict
            else None
        )

        logvar_pred = (
            latent_variables_dict["z_logvar"].detach().cpu().numpy()
            if self.configs.z_dim > 0 and "z_logvar" in latent_variables_dict
            else None
        )

        return rot_pred, trans_pred, conf_pred, logvar_pred

    def forward_pass(self, y_gt, tilt_ind, ind):
        if self.configs.verbose_time:
            torch.cuda.synchronize()

        start_time_ctf = time.time()
        if tilt_ind is not None:
            ctf_local = self.get_ctfs_at(tilt_ind)
        else:
            ctf_local = self.get_ctfs_at(ind)

        if self.configs.subtomo_averaging:
            ctf_local = ctf_local.reshape(-1, self.image_count, *ctf_local.shape[1:])

        if self.configs.verbose_time:
            torch.cuda.synchronize()
            self.run_times["ctf"].append(time.time() - start_time_ctf)

        # forward pass
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

        if self.configs.subtomo_averaging:
            tilt_ind = tilt_ind.reshape(y_gt.shape[0:2])

        in_dict = {
            "y": y_gt,
            "index": ind,
            "tilt_index": tilt_ind,
            "ctf": ctf_local,
        }

        if in_dict["tilt_index"] is None:
            in_dict["tilt_index"] = in_dict["index"]
        else:
            in_dict["tilt_index"] = in_dict["tilt_index"].reshape(-1)

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

        if self.configs.subtomo_averaging and self.configs.dose_exposure_correction:
            mask = self.model.output_mask.binary_mask
            a_pix = self.ctf_params[0, 0]

            dose_filters = self.data.get_dose_filters(
                in_dict["tilt_index"].reshape(-1), self.lattice, a_pix
            ).reshape(*y_pred.shape[:2], -1)

            y_pred *= dose_filters[..., mask]

        return latent_variables_dict, y_pred, y_gt_processed

    def loss(self, y_pred, y_gt, latent_variables_dict):
        """
        y_pred: [batch_size(, n_tilts), n_pts]
        y_gt: [batch_size(, n_tilts), n_pts]
        """
        all_losses = {}

        # data loss
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
                self.model.output_mask.binary_mask,
                self.resolution,
            )
            total_loss += self.configs.l2_smoothness_regularizer * smoothness_loss
            all_losses["L2 Smoothness Loss"] = smoothness_loss.item()

        return total_loss, all_losses

    def make_epoch_summary(self):
        summary.make_img_summary(
            self.writer,
            self.in_dict_last,
            self.y_pred_last,
            self.model.output_mask,
            self.current_epoch,
        )

        # conformation
        pca = None
        if self.configs.z_dim > 0:
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

            pca = summary.make_conf_summary(
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

        if self.configs.pose is not None:
            poses_gt = cryodrgn.utils.load_pkl(self.configs.pose)

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

        predicted_rots = self.predicted_rots[mask_tilt_idx]
        predicted_trans = (
            self.predicted_trans[mask_tilt_idx]
            if self.predicted_trans is not None
            else None
        )

        summary.make_pose_summary(
            self.writer,
            predicted_rots,
            predicted_trans,
            rotmat_gt,
            trans_gt,
            self.current_epoch,
            shift=shift,
        )

        return pca

    def make_batch_summary(self) -> None:
        self.logger.info(
            f"# [Train Epoch: {self.current_epoch}/{self.num_epochs - 1}] "
            f"[{self.current_epoch_particles_count}"
            f"/{self.particle_count} particles]"
        )

        if hasattr(self.model.output_mask, "current_radius"):
            self.current_losses["Mask Radius"] = self.model.output_mask.current_radius
        if self.model.trans_search_factor is not None:
            self.current_losses["Trans. Search Factor"] = self.model.trans_search_factor

        summary.make_scalar_summary(
            self.writer,
            self.current_losses,
            self.total_images_seen
        )

        if self.configs.verbose_time:
            for key in self.run_times.keys():
                self.logger.info(
                    f"{key} time: {np.mean(np.array(self.run_times[key]))}"
                )

    def save_latents(self):
        """Write model's latent variables to file."""
        out_pose = os.path.join(self.configs.outdir, f"pose.{self.current_epoch}.pkl")

        if self.configs.no_trans:
            with open(out_pose, "wb") as f:
                pickle.dump(self.predicted_rots, f)
        else:
            with open(out_pose, "wb") as f:
                pickle.dump((self.predicted_rots, self.predicted_trans), f)

        if self.configs.z_dim > 0:
            out_conf = os.path.join(
                self.configs.outdir, f"conf.{self.current_epoch}.pkl"
            )
            with open(out_conf, "wb") as f:
                pickle.dump(self.predicted_conf, f)

    def save_volume(self):
        """Write reconstructed volume to file."""
        out_mrc = os.path.join(
            self.configs.outdir, f"reconstruct.{self.current_epoch}.mrc"
        )

        self.model.hypervolume.eval()
        if hasattr(self.model, "conf_cnn"):
            if hasattr(self.model, "conf_regressor"):
                self.model.conf_cnn.eval()
                self.model.conf_regressor.eval()

        if hasattr(self.model, "pose_table"):
            self.model.pose_table.eval()
        if hasattr(self.model, "conf_table"):
            self.model.conf_table.eval()

        if self.configs.z_dim > 0:
            zval = self.predicted_conf[0].reshape(-1)
        else:
            zval = None

        vol = -1.0 * self.model.eval_volume(self.data.norm, zval=zval)
        mrc.write(out_mrc, vol.astype(np.float32))

    # TODO: weights -> model and reconstruct -> volume for output labels?
    def save_model(self):
        """Write model state to file."""
        out_weights = os.path.join(
            self.configs.outdir, f"weights.{self.current_epoch}.pkl"
        )

        optimizers_state_dict = {}
        for key in self.optimizers.keys():
            optimizers_state_dict[key] = self.optimizers[key].state_dict()

        saved_objects = {
            "epoch": self.current_epoch,
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

        if hasattr(self.model.output_mask, "current_radius"):
            saved_objects["output_mask_radius"] = self.model.output_mask.current_radius

        torch.save(saved_objects, out_weights)

    @property
    def is_in_pose_search_step(self) -> bool:
        in_pose_search = False

        if not self.configs.use_gt_poses:
            in_pose_search = 0 <= self.current_epoch < self.epochs_pose_search

        return in_pose_search
