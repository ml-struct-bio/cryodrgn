import os
import shutil
import pickle
import yaml
import logging
from datetime import datetime as dt
from collections import OrderedDict
import numpy as np
import time

import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader
from torch.utils.tensorboard import SummaryWriter

import cryodrgn
from cryodrgn import mrc
from cryodrgn import utils
from cryodrgn import dataset_fast
from cryodrgn import dataset
from cryodrgn import ctf
from cryodrgn.lattice import Lattice
from cryodrgn.losses import kl_divergence_conf, l1_regularizer, l2_frequency_bias
from cryodrgn.models.amortized_inference import DRGNai, MyDataParallel
from cryodrgn.mask import CircularMask, FrequencyMarchingMask
from cryodrgn.models.config import ModelConfigurations


class AmortizedInferenceTrainer:
    """An engine for training the reconstruction model on particle data.

    Attributes
    ----------
    configs (TrainingConfigurations):   Values of all parameters that can be
                                        set by the user.

    n_particles_dataset (int):  The number of picked particles in the data.
    pretraining (bool):     Whether we are in the pretraining stage.
    epoch (int):    Which training epoch the model is in.

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

    def __init__(self, config_vals: dict) -> None:
        self.logger = logging.getLogger(__name__)
        self.configs = AmortizedInferenceConfigurations(config_vals)

        # output directory
        if os.path.exists(self.configs.outdir):
            self.logger.warning(
                "Output directory already exists." "Renaming the old one."
            )

            newdir = self.configs.outdir + "_old"
            if os.path.exists(newdir):
                self.logger.warning("Must delete the previously " "saved old output.")
                shutil.rmtree(newdir)

            os.rename(self.configs.outdir, newdir)

        os.makedirs(self.configs.outdir)
        self.logger.addHandler(
            logging.FileHandler(os.path.join(self.configs.outdir, "training.log"))
        )

        n_gpus = torch.cuda.device_count()
        self.logger.info(f"Number of available gpus: {n_gpus}")
        self.n_prcs = max(n_gpus, 1)

        self.batch_size_known_poses = self.configs.batch_size_known_poses * self.n_prcs
        self.batch_size_hps = self.configs.batch_size_hps * self.n_prcs
        self.batch_size_sgd = self.configs.batch_size_sgd * self.n_prcs

        np.random.seed(self.configs.seed)
        torch.manual_seed(self.configs.seed)

        # set the device
        self.use_cuda = torch.cuda.is_available()
        self.device = torch.device("cuda:0" if self.use_cuda else "cpu")
        self.logger.info(f"Use cuda {self.use_cuda}")

        # tensorboard writer
        self.summaries_dir = os.path.join(self.configs.outdir, "summaries")
        os.makedirs(self.summaries_dir, exist_ok=True)
        self.writer = SummaryWriter(self.summaries_dir)
        self.logger.info("Will write tensorboard summaries " f"in {self.summaries_dir}")

        # load the optional index used to filter particles
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

        # load the particles
        self.logger.info("Creating dataset")
        if self.configs.fast_dataloading:
            if not self.configs.subtomogram_averaging:
                self.data = dataset_fast.ImageDataset(
                    self.configs.particles,
                    max_threads=self.configs.max_threads,
                    lazy=True,
                    poses_gt_pkl=self.configs.pose,
                    resolution_input=self.configs.resolution_encoder,
                    window_r=self.configs.window_radius_gt_real,
                    datadir=self.configs.datadir,
                    no_trans=self.configs.no_trans,
                )

            else:
                self.data = dataset_fast.TiltSeriesData(
                    self.configs.particles,
                    self.configs.n_tilts,
                    self.configs.angle_per_tilt,
                    window_r=self.configs.window_radius_gt_real,
                    datadir=self.configs.datadir,
                    max_threads=self.configs.max_threads,
                    dose_per_tilt=self.configs.dose_per_tilt,
                    device=self.device,
                    poses_gt_pkl=self.configs.pose,
                    tilt_axis_angle=self.configs.tilt_axis_angle,
                    ind=self.index,
                    no_trans=self.configs.no_trans,
                )

        else:
            self.data = dataset.MRCData(
                self.configs.particles,
                max_threads=self.configs.max_threads,
                ind=self.index,
                lazy=self.configs.lazy,
                relion31=self.configs.relion31,
                poses_gt_pkl=self.configs.pose,
                resolution_input=self.configs.resolution_encoder,
                window_r=self.configs.window_radius_gt_real,
                datadir=self.configs.datadir,
                no_trans=self.configs.no_trans,
            )

        self.n_particles_dataset = self.data.N
        self.n_tilts_dataset = (
            self.data.Nt if self.configs.subtomogram_averaging else self.data.N
        )
        self.resolution = self.data.D

        # load ctf
        if self.configs.ctf is not None:
            self.logger.info(f"Loading ctf params from {self.configs.ctf}")

            ctf_params = ctf.load_ctf_for_training(
                self.resolution - 1, self.configs.ctf
            )

            if self.index is not None:
                self.logger.info("Filtering dataset")
                ctf_params = ctf_params[self.index]

            assert ctf_params.shape == (self.n_tilts_dataset, 8)
            if self.configs.subtomogram_averaging:
                ctf_params = np.concatenate(
                    (ctf_params, self.data.ctfscalefactor.reshape(-1, 1)),
                    axis=1,  # type: ignore
                )
                self.data.voltage = float(ctf_params[0, 4])

            self.ctf_params = torch.tensor(ctf_params)
            self.ctf_params = self.ctf_params.to(self.device)

            if self.configs.subtomogram_averaging:
                self.data.voltage = float(self.ctf_params[0, 4])

        else:
            self.ctf_params = None

        # lattice
        self.logger.info("Building lattice")
        self.lattice = Lattice(self.resolution, extent=0.5, device=self.device)

        # output mask
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

        # pose search
        ps_params = None
        self.epochs_pose_search = 0

        if self.configs.n_imgs_pose_search > 0:
            self.epochs_pose_search = max(
                2, self.configs.n_imgs_pose_search // self.n_particles_dataset + 1
            )

        if self.epochs_pose_search > 0:
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
                "n_tilts_pose_search": self.configs.n_tilts_pose_search,
                "tilting_func": (
                    self.configs.data.get_tilting_func()
                    if self.configs.subtomogram_averaging
                    else None
                ),
                "average_over_tilts": self.configs.average_over_tilts,
            }

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

        self.model = DRGNai(
            self.lattice,
            self.output_mask,
            self.n_particles_dataset,
            self.n_tilts_dataset,
            cnn_params,
            conf_regressor_params,
            hyper_volume_params,
            resolution_encoder=self.configs.resolution_encoder,
            no_trans=self.configs.no_trans,
            use_gt_poses=self.configs.use_gt_poses,
            use_gt_trans=self.configs.use_gt_trans,
            will_use_point_estimates=will_use_point_estimates,
            ps_params=ps_params,
            verbose_time=self.configs.verbose_time,
            pretrain_with_gt_poses=self.configs.pretrain_with_gt_poses,
            n_tilts_pose_search=self.configs.n_tilts_pose_search,
        )

        # TODO: auto-loading from last weights file if load=True?
        # initialization from a previous checkpoint
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
            self.start_epoch = -1

        # move to gpu and parallelize
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

        # hypervolume
        hyper_volume_params = [{"params": list(self.model.hypervolume.parameters())}]

        self.optimizers["hypervolume"] = self.optim_types[
            self.configs.hypervolume_optimizer_type
        ](hyper_volume_params, lr=self.configs.lr)
        self.optimizer_types["hypervolume"] = self.configs.hypervolume_optimizer_type

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

        # initialization from a previous checkpoint
        if self.configs.load:
            checkpoint = torch.load(self.configs.load)

            for key in self.optimizers:
                self.optimizers[key].load_state_dict(
                    checkpoint["optimizers_state_dict"][key]
                )

        # dataloaders
        if self.configs.fast_dataloading:
            self.data_generator_pose_search = dataset_fast.make_dataloader(
                self.data,
                batch_size=self.batch_size_hps,
                num_workers=self.configs.num_workers,
                shuffler_size=self.configs.shuffler_size,
            )
            self.data_generator = dataset_fast.make_dataloader(
                self.data,
                batch_size=self.batch_size_known_poses,
                num_workers=self.configs.num_workers,
                shuffler_size=self.configs.shuffler_size,
            )
            self.data_generator_latent_optimization = dataset_fast.make_dataloader(
                self.data,
                batch_size=self.batch_size_sgd,
                num_workers=self.configs.num_workers,
                shuffler_size=self.configs.shuffler_size,
            )

        else:
            self.data_generator_pose_search = DataLoader(
                self.data,
                batch_size=self.batch_size_hps,
                shuffle=self.configs.shuffle,
                num_workers=self.configs.num_workers,
                drop_last=False,
            )
            self.data_generator = DataLoader(
                self.data,
                batch_size=self.batch_size_known_poses,
                shuffle=self.configs.shuffle,
                num_workers=self.configs.num_workers,
                drop_last=False,
            )
            self.data_generator_latent_optimization = DataLoader(
                self.data,
                batch_size=self.batch_size_sgd,
                shuffle=self.configs.shuffle,
                num_workers=self.configs.num_workers,
                drop_last=False,
            )

        # save configurations
        self.configs.write(os.path.join(self.configs.outdir, "train-configs.yaml"))

        epsilon = 1e-8
        # booleans
        self.log_latents = False
        self.pose_only = True
        self.pretraining = False
        self.is_in_pose_search_step = False
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

        self.n_particles_pretrain = (
            self.configs.n_imgs_pretrain
            if self.configs.n_imgs_pretrain >= 0
            else self.n_particles_dataset
        )

        # placeholders for predicted latent variables,
        # last input/output batch, losses
        self.in_dict_last = None
        self.y_pred_last = None

        self.predicted_rots = np.empty((self.n_tilts_dataset, 3, 3))
        self.predicted_trans = (
            np.empty((self.n_tilts_dataset, 2)) if not self.configs.no_trans else None
        )
        self.predicted_conf = (
            np.empty((self.n_particles_dataset, self.configs.z_dim))
            if self.configs.z_dim > 0
            else None
        )

        self.predicted_logvar = (
            np.empty((self.n_particles_dataset, self.configs.z_dim))
            if self.configs.z_dim > 0 and self.configs.variational_het
            else None
        )

        self.mask_particles_seen_at_last_epoch = np.zeros(self.n_particles_dataset)
        self.mask_tilts_seen_at_last_epoch = np.zeros(self.n_tilts_dataset)

        # counters
        self.epoch = 0
        self.run_times = {phase: [] for phase in self.run_phases}
        self.current_epoch_particles_count = 0
        self.total_batch_count = 0
        self.total_particles_count = 0
        self.batch_idx = 0
        self.cur_loss = None

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
            np.zeros((self.n_particles_dataset, self.configs.z_dim))
            if self.configs.z_dim > 0
            else None
        )

        self.total_batch_count = 0
        self.total_particles_count = 0

        self.epoch = self.start_epoch - 1
        for epoch in range(self.start_epoch, self.num_epochs):
            te = dt.now()

            self.mask_particles_seen_at_last_epoch = np.zeros(self.n_particles_dataset)
            self.mask_tilts_seen_at_last_epoch = np.zeros(self.n_tilts_dataset)

            self.epoch += 1
            self.current_epoch_particles_count = 0
            self.optimized_modules = ["hypervolume"]

            self.pose_only = (
                self.total_particles_count < self.configs.pose_only_phase
                or self.configs.z_dim == 0
                or epoch < 0
            )
            self.pretraining = self.epoch < 0

            if not self.configs.use_gt_poses:
                self.is_in_pose_search_step = 0 <= epoch < self.epochs_pose_search
                self.use_point_estimates = epoch >= max(0, self.epochs_pose_search)

            n_max_particles = self.n_particles_dataset
            data_generator = self.data_generator

            # pre-training
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

                    if self.configs.refine_gt_poses:
                        self.logger.info("Initializing pose table from ground truth")

                        poses_gt = utils.load_pkl(self.configs.pose)
                        if poses_gt[0].ndim == 3:
                            # contains translations
                            rotmat_gt = torch.tensor(poses_gt[0]).float()
                            trans_gt = torch.tensor(poses_gt[1]).float()
                            trans_gt *= self.resolution

                            if self.index is not None:
                                rotmat_gt = rotmat_gt[self.index]
                                trans_gt = trans_gt[self.index]

                        else:
                            rotmat_gt = torch.tensor(poses_gt).float()
                            trans_gt = None

                            if self.index is not None:
                                rotmat_gt = rotmat_gt[self.index]

                        self.model.pose_table.initialize(rotmat_gt, trans_gt)

                    else:
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
                assert self.configs.use_gt_poses

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
                                utils.load_pkl(self.configs.initial_conf)
                            )

                        self.model.to(self.device)

                    self.optimized_modules.append("conf_table")

            will_make_summary = (
                (
                    self.configs.log_heavy_interval
                    and epoch % self.configs.log_heavy_interval == 0
                )
                or self.is_in_pose_search_step
                or self.pretraining
            )
            self.log_latents = will_make_summary

            if will_make_summary:
                self.logger.info("Will make a full summary at the end of this epoch")

            for key in self.run_times.keys():
                self.run_times[key] = []

            end_time = time.time()
            self.cur_loss = 0

            # inner loop
            for batch_idx, in_dict in enumerate(data_generator):
                self.batch_idx = batch_idx

                # with torch.autograd.detect_anomaly():
                self.train_step(in_dict, end_time=end_time)
                if self.configs.verbose_time:
                    torch.cuda.synchronize()

                end_time = time.time()

                if self.current_epoch_particles_count > n_max_particles:
                    break

            total_loss = self.cur_loss / n_max_particles
            self.logger.info(
                f"# =====> SGD Epoch: {self.epoch} "
                f"finished in {dt.now() - te}; "
                f"total loss = {format(total_loss, '.6f')}"
            )

            # image and pose summary
            if will_make_summary:
                self.make_heavy_summary()
                self.save_latents()
                self.save_volume()
                self.save_model()

            # update output mask -- epoch-based scaling
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

        # update output mask -- image-based scaling
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

        if "tilt_index" not in in_dict:
            in_dict["tilt_index"] = in_dict["index"]
        else:
            in_dict["tilt_index"] = in_dict["tilt_index"].reshape(-1)

        ind_tilt = in_dict["tilt_index"]
        self.total_batch_count += 1
        batch_size = len(y_gt)
        self.total_particles_count += batch_size
        self.current_epoch_particles_count += batch_size

        # move to gpu
        if self.configs.verbose_time:
            torch.cuda.synchronize()
        start_time_gpu = time.time()

        for key in in_dict.keys():
            in_dict[key] = in_dict[key].to(self.device)
        if self.configs.verbose_time:
            torch.cuda.synchronize()
            self.run_times["to_gpu"].append(time.time() - start_time_gpu)

        # zero grad
        for key in self.optimized_modules:
            self.optimizers[key].zero_grad()

        # forward pass
        latent_variables_dict, y_pred, y_gt_processed = self.forward_pass(in_dict)

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
        self.cur_loss += total_loss.item() * len(ind)

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

            # log
            if self.use_cuda:
                ind = ind.cpu()
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

        # scalar summary
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
            if self.configs.z_dim > 0 and "z" in latent_variables_dict
            else None
        )

        logvar_pred = (
            latent_variables_dict["z_logvar"].detach().cpu().numpy()
            if self.configs.z_dim > 0 and "z_logvar" in latent_variables_dict
            else None
        )

        return rot_pred, trans_pred, conf_pred, logvar_pred

    def forward_pass(self, in_dict):
        if self.configs.verbose_time:
            torch.cuda.synchronize()

        start_time_ctf = time.time()
        ctf_local = self.get_ctfs_at(in_dict["tilt_index"])

        if self.configs.subtomogram_averaging:
            ctf_local = ctf_local.reshape(-1, self.n_tilts, *ctf_local.shape[1:])

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

        if self.configs.subtomogram_averaging:
            in_dict["tilt_index"] = in_dict["tilt_index"].reshape(
                *in_dict["y"].shape[0:2]
            )

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

        if self.configs.subtomogram_averaging and self.configs.dose_exposure_correction:
            mask = self.output_mask.binary_mask
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
                self.output_mask.binary_mask,
                self.resolution,
            )
            total_loss += self.configs.l2_smoothness_regularizer * smoothness_loss
            all_losses["L2 Smoothness Loss"] = smoothness_loss.item()

        return total_loss, all_losses

    def make_heavy_summary(self):
        cryodrgn.summary.make_img_summary(
            self.writer,
            self.in_dict_last,
            self.y_pred_last,
            self.output_mask,
            self.epoch,
        )

        # conformation
        pca = None
        if self.configs.z_dim > 0:
            labels = None

            if self.configs.labels is not None:
                labels = utils.load_pkl(self.configs.labels)

                if self.index is not None:
                    labels = labels[self.index]

            if self.mask_particles_seen_at_last_epoch is not None:
                mask_idx = self.mask_particles_seen_at_last_epoch > 0.5
            else:
                mask_idx = np.ones((self.n_particles_dataset,), dtype=bool)

            predicted_conf = self.predicted_conf[mask_idx]
            labels = labels[mask_idx] if labels is not None else None
            logvar = (
                self.predicted_logvar[mask_idx]
                if self.predicted_logvar is not None
                else None
            )

            pca = cryodrgn.summary.make_conf_summary(
                self.writer,
                predicted_conf,
                self.epoch,
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
            mask_tilt_idx = np.ones((self.n_tilts_dataset,), dtype=bool)

        if self.configs.pose is not None:
            poses_gt = utils.load_pkl(self.configs.pose)

            if poses_gt[0].ndim == 3:
                # contains translations
                rotmat_gt = torch.tensor(poses_gt[0]).float()
                trans_gt = torch.tensor(poses_gt[1]).float() * self.resolution

                if self.index is not None:
                    rotmat_gt = rotmat_gt[self.index]
                    trans_gt = trans_gt[self.index]

            else:
                rotmat_gt = torch.tensor(poses_gt).float()
                trans_gt = None
                assert not shift, "Shift activated but trans not given in gt"

                if self.index is not None:
                    rotmat_gt = rotmat_gt[self.index]

            rotmat_gt = rotmat_gt[mask_tilt_idx]
            trans_gt = trans_gt[mask_tilt_idx] if trans_gt is not None else None

        predicted_rots = self.predicted_rots[mask_tilt_idx]
        predicted_trans = (
            self.predicted_trans[mask_tilt_idx]
            if self.predicted_trans is not None
            else None
        )

        cryodrgn.summary.make_pose_summary(
            self.writer,
            predicted_rots,
            predicted_trans,
            rotmat_gt,
            trans_gt,
            self.epoch,
            shift=shift,
        )

        return pca

    def make_light_summary(self, all_losses):
        self.logger.info(
            f"# [Train Epoch: {self.epoch}/{self.num_epochs - 1}] "
            f"[{self.current_epoch_particles_count}"
            f"/{self.n_particles_dataset} particles]"
        )

        if hasattr(self.output_mask, "current_radius"):
            all_losses["Mask Radius"] = self.output_mask.current_radius

        if self.model.trans_search_factor is not None:
            all_losses["Trans. Search Factor"] = self.model.trans_search_factor

        cryodrgn.summary.make_scalar_summary(
            self.writer, all_losses, self.total_particles_count
        )

        if self.configs.verbose_time:
            for key in self.run_times.keys():
                self.logger.info(
                    f"{key} time: {np.mean(np.array(self.run_times[key]))}"
                )

    def save_latents(self):
        """Write model's latent variables to file."""
        out_pose = os.path.join(self.configs.outdir, f"pose.{self.epoch}.pkl")

        if self.configs.no_trans:
            with open(out_pose, "wb") as f:
                pickle.dump(self.predicted_rots, f)
        else:
            with open(out_pose, "wb") as f:
                pickle.dump((self.predicted_rots, self.predicted_trans), f)

        if self.configs.z_dim > 0:
            out_conf = os.path.join(self.configs.outdir, f"conf.{self.epoch}.pkl")
            with open(out_conf, "wb") as f:
                pickle.dump(self.predicted_conf, f)

    def save_volume(self):
        """Write reconstructed volume to file."""
        out_mrc = os.path.join(self.configs.outdir, f"reconstruct.{self.epoch}.mrc")

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
        out_weights = os.path.join(self.configs.outdir, f"weights.{self.epoch}.pkl")

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


class AmortizedInferenceConfigurations(ModelConfigurations):

    __slots__ = (
        "particles",
        "ctf",
        "pose",
        "dataset",
        "server",
        "datadir",
        "ind",
        "labels",
        "relion31",
        "no_trans",
        "use_gt_poses",
        "refine_gt_poses",
        "use_gt_trans",
        "load",
        "initial_conf",
        "log_interval",
        "log_heavy_interval",
        "verbose_time",
        "shuffle",
        "lazy",
        "num_workers",
        "max_threads",
        "fast_dataloading",
        "shuffler_size",
        "batch_size_known_poses",
        "batch_size_hps",
        "batch_size_sgd",
        "hypervolume_optimizer_type",
        "pose_table_optimizer_type",
        "conf_table_optimizer_type",
        "conf_encoder_optimizer_type",
        "lr",
        "lr_pose_table",
        "lr_conf_table",
        "lr_conf_encoder",
        "wd",
        "n_imgs_pose_search",
        "epochs_sgd",
        "pose_only_phase",
        "output_mask",
        "add_one_frequency_every",
        "n_frequencies_per_epoch",
        "max_freq",
        "window_radius_gt_real",
        "l_start_fm",
        "beta_conf",
        "trans_l1_regularizer",
        "l2_smoothness_regularizer",
        "variational_het",
        "z_dim",
        "std_z_init",
        "use_conf_encoder",
        "depth_cnn",
        "channels_cnn",
        "kernel_size_cnn",
        "resolution_encoder",
        "explicit_volume",
        "hypervolume_layers",
        "hypervolume_dim",
        "pe_type",
        "pe_dim",
        "feat_sigma",
        "hypervolume_domain",
        "pe_type_conf",
        "n_imgs_pretrain",
        "pretrain_with_gt_poses",
        "l_start",
        "l_end",
        "n_iter",
        "t_extent",
        "t_n_grid",
        "t_x_shift",
        "t_y_shift",
        "no_trans_search_at_pose_search",
        "n_kept_poses",
        "base_healpy",
        "subtomogram_averaging",
        "n_tilts",
        "dose_per_tilt",
        "angle_per_tilt",
        "n_tilts_pose_search",
        "average_over_tilts",
        "tilt_axis_angle",
        "dose_exposure_correction",
        "seed",
    )
    defaults = OrderedDict(
        {
            # dataset
            "particles": None,
            "ctf": None,
            "pose": None,
            "dataset": None,
            "server": None,
            "datadir": None,
            "ind": None,
            "labels": None,
            "relion31": False,
            "no_trans": False,
            # initialization
            "use_gt_poses": False,
            "refine_gt_poses": False,
            "use_gt_trans": False,
            "load": None,
            "initial_conf": None,
            # logging
            "log_interval": 10000,
            "log_heavy_interval": 5,
            "verbose_time": False,
            # data loading
            "shuffle": True,
            "lazy": False,
            "num_workers": 2,
            "max_threads": 16,
            "fast_dataloading": False,
            "shuffler_size": 32768,
            "batch_size_known_poses": 32,
            "batch_size_hps": 8,
            "batch_size_sgd": 256,
            # optimizers
            "hypervolume_optimizer_type": "adam",
            "pose_table_optimizer_type": "adam",
            "conf_table_optimizer_type": "adam",
            "conf_encoder_optimizer_type": "adam",
            "lr": 1.0e-4,
            "lr_pose_table": 1.0e-3,
            "lr_conf_table": 1.0e-2,
            "lr_conf_encoder": 1.0e-4,
            "wd": 0.0,
            # scheduling
            "n_imgs_pose_search": 500000,
            "epochs_sgd": 100,
            "pose_only_phase": 0,
            # masking
            "output_mask": "circ",
            "add_one_frequency_every": 100000,
            "n_frequencies_per_epoch": 10,
            "max_freq": None,
            "window_radius_gt_real": 0.85,
            "l_start_fm": 12,
            # loss
            "beta_conf": 0.0,
            "trans_l1_regularizer": 0.0,
            "l2_smoothness_regularizer": 0.0,
            # conformations
            "variational_het": False,
            "z_dim": 4,
            "std_z_init": 0.1,
            "use_conf_encoder": False,
            "depth_cnn": 5,
            "channels_cnn": 32,
            "kernel_size_cnn": 3,
            "resolution_encoder": None,
            # hypervolume
            "explicit_volume": False,
            "hypervolume_layers": 3,
            "hypervolume_dim": 256,
            "pe_type": "gaussian",
            "pe_dim": 64,
            "feat_sigma": 0.5,
            "hypervolume_domain": "hartley",
            "pe_type_conf": None,
            # pre-training
            "n_imgs_pretrain": 10000,
            "pretrain_with_gt_poses": False,
            # pose search
            "l_start": 12,
            "l_end": 32,
            "n_iter": 4,
            "t_extent": 20.0,
            "t_n_grid": 7,
            "t_x_shift": 0.0,
            "t_y_shift": 0.0,
            "no_trans_search_at_pose_search": False,
            "n_kept_poses": 8,
            "base_healpy": 2,
            # subtomogram averaging
            "subtomogram_averaging": False,
            "n_tilts": 11,
            "dose_per_tilt": 2.93,
            "angle_per_tilt": 3.0,
            "n_tilts_pose_search": 11,
            "average_over_tilts": False,
            "tilt_axis_angle": 0.0,
            "dose_exposure_correction": True,
            # others
            "seed": -1,
            "palette_type": None,
        }
    )

    quick_defaults = {
        "capture_setup": {
            "spa": dict(),
            "et": {
                "subtomogram_averaging": True,
                "fast_dataloading": True,
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

    def __init__(self, config_vals: dict):
        super().__init__(config_vals)

        if self.seed < 0:
            self.seed = np.random.randint(0, 10000)

        # process the quick_config parameter
        if self.quick_config is not None:
            for key, value in self.quick_config.items():
                if key not in self.quick_defaults:
                    raise ValueError(
                        "Unrecognized parameter " f"shortcut field `{key}`!"
                    )

                if value not in self.quick_defaults[key]:
                    raise ValueError(
                        "Unrecognized parameter shortcut label "
                        f"`{value}` for field `{key}`!"
                    )

                for _key, _value in self.quick_defaults[key][value].items():
                    if _key not in self.defaults:
                        raise ValueError(
                            "Unrecognized configuration " f"parameter `{key}`!"
                        )

                    # given parameters have priority
                    if _key not in config_vals:
                        setattr(self, _key, _value)

        if not self.outdir:
            raise ValueError("Must specify an outdir!")

        if self.explicit_volume and self.z_dim >= 1:
            raise ValueError(
                "Explicit volumes do not support " "heterogeneous reconstruction."
            )

        if self.dataset is None:
            if self.particles is None:
                raise ValueError(
                    "As dataset was not specified, please " "specify particles!"
                )

            if self.ctf is None:
                raise ValueError("As dataset was not specified, please " "specify ctf!")

        if self.hypervolume_optimizer_type not in {"adam"}:
            raise ValueError(
                "Invalid value "
                f"`{self.hypervolume_optimizer_type}` "
                "for hypervolume_optimizer_type!"
            )

        if self.pose_table_optimizer_type not in {"adam", "lbfgs"}:
            raise ValueError(
                "Invalid value "
                f"`{self.pose_table_optimizer_type}` "
                "for pose_table_optimizer_type!"
            )

        if self.conf_table_optimizer_type not in {"adam", "lbfgs"}:
            raise ValueError(
                "Invalid value "
                f"`{self.conf_table_optimizer_type}` "
                "for conf_table_optimizer_type!"
            )

        if self.conf_encoder_optimizer_type not in {"adam"}:
            raise ValueError(
                "Invalid value "
                f"`{self.conf_encoder_optimizer_type}` "
                "for conf_encoder_optimizer_type!"
            )

        if self.output_mask not in {"circ", "frequency_marching"}:
            raise ValueError("Invalid value " f"{self.output_mask} for output_mask!")

        if self.pe_type not in {"gaussian"}:
            raise ValueError(f"Invalid value {self.pe_type} for pe_type!")

        if self.pe_type_conf not in {None, "geom"}:
            raise ValueError(f"Invalid value {self.pe_type_conf} " "for pe_type_conf!")

        if self.hypervolume_domain not in {"hartley"}:
            raise ValueError(
                f"Invalid value {self.hypervolume_domain} " "for hypervolume_domain."
            )

        if self.n_imgs_pose_search < 0:
            raise ValueError("n_imgs_pose_search must be greater than 0!")

        if self.use_conf_encoder and self.initial_conf:
            raise ValueError(
                "Conformations cannot be initialized " "when using an encoder!"
            )

        if self.use_gt_trans and self.pose is None:
            raise ValueError("Poses must be specified to use GT translations!")

        if self.refine_gt_poses:
            self.n_imgs_pose_search = 0

            if self.pose is None:
                raise ValueError("Initial poses must be specified " "to be refined!")

        if self.subtomogram_averaging:
            self.fast_dataloading = True

            # TODO: Implement conformation encoder for subtomogram averaging.
            if self.use_conf_encoder:
                raise ValueError(
                    "Conformation encoder is not implemented "
                    "for subtomogram averaging!"
                )

            # TODO: Implement translation search for subtomogram averaging.
            if not (self.use_gt_poses or self.use_gt_trans or self.t_extent == 0.0):
                raise ValueError(
                    "Translation search is not implemented "
                    "for subtomogram averaging!"
                )

            if self.average_over_tilts and self.n_tilts_pose_search % 2 == 0:
                raise ValueError(
                    "n_tilts_pose_search must be odd " "to use average_over_tilts!"
                )

            if self.n_tilts_pose_search > self.n_tilts:
                raise ValueError("n_tilts_pose_search must be " "smaller than n_tilts!")
        if self.use_gt_poses:
            # "poses" include translations
            self.use_gt_trans = True

            if self.pose is None:
                raise ValueError("Ground truth poses must be specified!")

        if self.no_trans:
            self.t_extent = 0.0
        if self.t_extent == 0.0:
            self.t_n_grid = 1

        if self.dataset:
            with open(os.environ.get("DRGNAI_DATASETS"), "r") as f:
                paths = yaml.safe_load(f)

            self.particles = paths[self.dataset]["particles"]
            self.ctf = paths[self.dataset]["ctf"]

            if self.pose is None and "pose" in paths[self.dataset]:
                self.pose = paths[self.dataset]["pose"]
            if "datadir" in paths[self.dataset]:
                self.datadir = paths[self.dataset]["datadir"]
            if "labels" in paths[self.dataset]:
                self.labels = paths[self.dataset]["labels"]
            if self.ind is None and "ind" in paths[self.dataset]:
                self.ind = paths[self.dataset]["ind"]
            if "dose_per_tilt" in paths[self.dataset]:
                self.dose_per_tilt = paths[self.dataset]["dose_per_tilt"]
