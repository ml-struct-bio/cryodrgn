"""Base classes for model training engines."""
import argparse
import os
import shutil
import sys
from collections import OrderedDict
from abc import ABC, abstractmethod
from typing import Any
from typing_extensions import Self
import yaml
from datetime import datetime as dt
import logging

import numpy as np
import torch
from torch import nn
import cryodrgn.utils
from cryodrgn import ctf
from cryodrgn.dataset import ImageDataset, TiltSeriesData, make_dataloader
from cryodrgn.lattice import Lattice


class BaseConfigurations(ABC):
    """Base class for sets of model configuration parameters."""

    # a parameter belongs to this set if and only if it has a default value
    # defined in these variables, ordering makes e.g. printing easier for user
    __slots__ = ("outdir", "verbose", "seed")
    default_values = OrderedDict(
        {"outdir": os.getcwd(), "verbose": 0, "seed": np.random.randint(0, 10000)}
    )
    quick_config = OrderedDict()

    def __init__(self, config_vals: dict[str, Any]) -> None:
        if "test_installation" in config_vals and config_vals["test_installation"]:
            print("Installation was successful!")
            sys.exit()

        if "verbose" in config_vals:
            if not isinstance(config_vals["verbose"], int):
                raise ValueError(
                    f"Given verbosity `{config_vals['verbose']}` is not an integer!"
                )
            if config_vals["verbose"] < 0:
                raise ValueError(
                    f"Given verbosity `{config_vals['verbose']}` is not positive!"
                )

        if "seed" in config_vals and not isinstance(config_vals["seed"], int):
            raise ValueError(
                "Configuration `seed` must be given as an integer, "
                f"given `{config_vals['seed']}` instead!"
            )

        if "outdir" in config_vals and config_vals["outdir"] is not None:
            config_vals["outdir"] = os.path.abspath(config_vals["outdir"])

        # process the quick_config parameter
        if "quick_config" in config_vals and config_vals["quick_config"] is not None:
            for key, value in config_vals["quick_config"].items():
                if key not in self.quick_config:
                    raise ValueError(
                        "Unrecognized parameter " f"shortcut field `{key}`!"
                    )

                if value not in self.quick_config[key]:
                    raise ValueError(
                        "Unrecognized parameter shortcut label "
                        f"`{value}` for field `{key}`!"
                    )

                for _key, _value in self.quick_config[key][value].items():
                    if _key not in self.defaults:
                        raise ValueError(
                            "Unrecognized configuration " f"parameter `{key}`!"
                        )

                    # given parameters have priority
                    if _key not in config_vals:
                        config_vals[_key] = _value

        for key in config_vals:
            if key not in self.defaults:
                raise ValueError(f"Unrecognized configuration parameter `{key}`!")

        # an attribute is created for every entry in the defaults dictionary
        for key, value in self.defaults.items():
            if key in config_vals:
                setattr(self, key, config_vals[key])

            # if not in given parameters, use defaults
            else:
                setattr(self, key, value)

    @classmethod
    @property
    def defaults(cls) -> dict[str, Any]:
        def_values = cls.default_values

        for c in cls.__bases__:
            if hasattr(c, "defaults"):
                def_values.update(c.defaults)

        return def_values

    def __iter__(self):
        return iter((par, getattr(self, par)) for par in self.defaults)

    def __str__(self):
        return "\n".join([f"{par}{str(val):>20}" for par, val in self])

    def write(self, fl: str) -> None:
        """Saving configurations to file using the original order."""

        with open(fl, "w") as f:
            yaml.dump(dict(self), f, default_flow_style=False, sort_keys=False)


class BaseTrainer(ABC):
    """Abstract base class for reconstruction model training engines.

    Attributes
    ----------
    outdir(str):    Path to where experiment output is stored.
    quick_config(dict):    Configuration shortcuts.
    """

    config_cls = BaseConfigurations

    @classmethod
    @property
    def parameters(cls) -> list:
        """The user-set parameters governing the behaviour of this model."""
        return list(cls.config_cls.defaults.keys())

    @classmethod
    def parse_args(cls, args: argparse.Namespace) -> Self:
        """Utility for initializing using a namespace as opposed to a dictionary."""
        return cls(
            {
                par: (
                    getattr(args, par)
                    if hasattr(args, par)
                    else cls.config_cls.defaults[par]
                )
                for par in tuple(cls.config_cls.defaults)
            }
        )

    def __init__(self, configs: dict[str, Any]) -> None:
        if "load" in configs and configs["load"] == "latest":
            configs = self.get_latest_configs()

        self.configs = self.config_cls(configs)
        self.outdir = self.configs.outdir
        self.verbose = self.configs.verbose
        np.random.seed(self.configs.seed)
        self.logger = logging.getLogger(__name__)

        # TODO: more sophisticated management of existing output folders
        if os.path.exists(self.outdir):
            self.logger.warning("Output folder already exists, renaming the old one.")
            newdir = self.outdir + "_old"

            if os.path.exists(newdir):
                self.logger.warning("Must delete the previously saved old output.")
                shutil.rmtree(newdir)

            os.rename(self.outdir, newdir)
        os.makedirs(self.outdir)

        if self.configs.verbose:
            self.logger.setLevel(logging.DEBUG)
        self.logger.addHandler(
            logging.FileHandler(os.path.join(self.outdir, "training.log"))
        )
        self.logger.info(str(configs))


class ModelConfigurations(BaseConfigurations):

    __slots__ = (
        "model",
        "particles",
        "ctf",
        "dataset",
        "datadir",
        "ind",
        "log_interval",
        "verbose",
        "load",
        "load_poses",
        "initial_conf",
        "checkpoint",
        "z_dim",
        "invert_data",
        "lazy",
        "window",
        "window_r",
        "shuffler_size",
        "max_threads",
        "tilt",
        "tilt_deg",
        "num_epochs",
        "batch_size",
        "weight_decay",
        "learning_rate",
        "data_norm",
        "multigpu",
        "pretrain",
        "enc_layers",
        "enc_dim",
        "encode_mode",
        "enc_mask",
        "use_real",
        "dec_layers",
        "dec_dim",
        "pe_type",
        "feat_sigma",
        "pe_dim",
        "hypervolume_domain",
        "activation",
        "subtomo_averaging",
        "no_trans",
    )
    default_values = OrderedDict(
        {
            "model": "amort",
            "particles": None,
            "ctf": None,
            "dataset": None,
            "datadir": None,
            "ind": None,
            "log_interval": 1000,
            "verbose": False,
            "load": None,
            "load_poses": None,
            "checkpoint": 1,
            "z_dim": None,
            "invert_data": True,
            "lazy": False,
            "window": True,
            "window_r": 0.85,
            "shuffler_size": 0,
            "max_threads": 16,
            "tilt": None,
            "tilt_deg": 45,
            "num_epochs": 30,
            "batch_size": 8,
            "weight_decay": 0,
            "learning_rate": 1e-4,
            "data_norm": None,
            "multigpu": False,
            "pretrain": 10000,
            "enc_layers": 3,
            "enc_dim": 256,
            "encode_mode": "resid",
            "enc_mask": None,
            "use_real": False,
            "dec_layers": 3,
            "dec_dim": 256,
            "pe_type": "gaussian",
            "feat_sigma": 0.5,
            "pe_dim": 64,
            "hypervolume_domain": "hartley",
            "activation": "relu",
            "subtomo_averaging": False,
            "no_trans": False,
        }
    )

    def __init__(self, config_vals: dict[str, Any]) -> None:
        if "model" in config_vals:
            if config_vals["model"] not in {"hps", "amort"}:
                raise ValueError(
                    f"Given model `{config_vals['model']}` not in currently supported "
                    f"model types:\n"
                    f"`amort` (cryoDRGN v4 amortized inference pose estimation)\n"
                    f"`hps` (cryoDRGN v3 hierarchical pose estimation)\n"
                )

        if "ind" in config_vals:
            if isinstance(config_vals["ind"], str):
                if not os.path.exists(config_vals["ind"]):
                    raise ValueError(
                        f"Subset indices file {config_vals['ind']} does not exist!"
                    )

        if "dataset" in config_vals and config_vals["dataset"]:
            paths_file = os.environ.get("DRGNAI_DATASETS")
            paths = cryodrgn.utils.load_yaml(paths_file)

            if config_vals["dataset"] not in paths:
                raise ValueError(
                    f"Given dataset label `{config_vals['dataset']}` not in list of "
                    f"datasets specified at `{paths_file}`!"
                )
            use_paths = paths[config_vals["dataset"]]

            config_vals["particles"] = use_paths["particles"]
            for k in ["ctf", "pose", "datadir", "labels", "ind", "dose_per_tilt"]:
                if k not in config_vals and k in use_paths:
                    if not os.path.exists(use_paths[k]):
                        raise ValueError(
                            f"Given {k} file `{use_paths[k]}` does not exist!"
                        )
                    config_vals[k] = use_paths[k]

        super().__init__(config_vals)


class ModelTrainer(BaseTrainer, ABC):

    config_cls = ModelConfigurations

    # options for optimizers to use
    optim_types = {"adam": torch.optim.Adam, "lbfgs": torch.optim.LBFGS}

    @abstractmethod
    def make_model(self, configs: dict[str, Any]) -> nn.Module:
        pass

    def __init__(self, configs: dict[str, Any]) -> None:
        super().__init__(configs)
        self.configs: ModelConfigurations

        # set the device
        torch.manual_seed(self.configs.seed)
        self.use_cuda = torch.cuda.is_available()
        self.device = torch.device("cuda" if self.use_cuda else "cpu")

        if self.use_cuda:
            self.logger.info(f"Using GPU {self.device}")
            self.n_prcs = max(torch.cuda.device_count(), 1)
            self.logger.info(f"Number of available gpus: {self.n_prcs}")
        else:
            self.logger.warning(f"No GPUs detected, using {self.device} instead!")
            self.n_prcs = 1

        # load index filter
        if self.configs.ind is not None:
            if isinstance(self.configs.ind, int):
                self.logger.info(f"Keeping the first {self.configs.ind} particles")
                self.ind = np.arange(self.configs.ind)
            else:
                self.logger.info(f"Filtering image dataset with {configs['ind']}")
                self.ind = cryodrgn.utils.load_yaml(self.configs.ind)
        else:
            self.ind = None

        # load dataset
        self.logger.info(f"Loading dataset from {self.configs.particles}")

        if self.configs.tilt is not None:
            self.tilt = torch.tensor(
                cryodrgn.utils.xrot(configs["tilt_deg"]).astype(np.float32),
                device=self.device,
            )
        else:
            self.tilt = None

        if not self.configs.subtomo_averaging:
            self.data = ImageDataset(
                mrcfile=self.configs.particles,
                norm=self.configs.data_norm,
                invert_data=self.configs.invert_data,
                ind=self.configs.ind,
                window=self.configs.window,
                keepreal=self.configs.use_real,
                datadir=self.configs.datadir,
                window_r=self.configs.window_r,
            )
        else:
            self.data = TiltSeriesData(
                tiltstar=self.configs.particles,
                ind=self.configs.ind,
                ntilts=self.configs.n_tilts,
                angle_per_tilt=self.configs.angle_per_tilt,
                window_r=self.configs.window_radius_gt_real,
                datadir=self.configs.datadir,
                max_threads=self.configs.max_threads,
                dose_per_tilt=self.configs.dose_per_tilt,
                device=self.device,
                poses_gt_pkl=self.configs.pose,
                tilt_axis_angle=self.configs.tilt_axis_angle,
                no_trans=self.configs.no_trans,
            )

        self.image_count = self.data.N
        self.particle_count = self.data.N
        self.resolution = self.data.D

        if self.configs.ctf:
            self.logger.info(f"Loading ctf params from {self.configs.ctf}")
            ctf_params = ctf.load_ctf_for_training(
                self.resolution - 1, self.configs.ctf
            )

            if self.ind is not None:
                ctf_params = ctf_params[self.ind]
            assert ctf_params.shape == (self.image_count, 8), ctf_params.shape

            if self.configs.subtomo_averaging:
                ctf_params = np.concatenate(
                    (ctf_params, self.data.ctfscalefactor.reshape(-1, 1)),
                    axis=1,  # type: ignore
                )
                self.data.voltage = float(ctf_params[0, 4])

            self.ctf_params = torch.tensor(ctf_params, device=self.device)
            if self.configs.subtomo_averaging:
                self.data.voltage = float(self.ctf_params[0, 4])

        else:
            self.ctf_params = None

        # lattice
        self.logger.info("Building lattice...")
        self.lattice = Lattice(self.resolution, extent=0.5, device=self.device)

        self.logger.info("Initializing model...")
        self.model = self.make_model()
        self.model.to(self.device)
        self.logger.info(self.model)
        param_count = sum(p.numel() for p in self.model.parameters() if p.requires_grad)
        self.logger.info(f"{param_count} parameters in model")
        self.current_epoch = None

        # TODO: auto-loading from last weights file if load=True?
        # initialization from a previous checkpoint
        if self.configs.load:
            self.logger.info(f"Loading checkpoint from {self.configs.load}")
            checkpoint = torch.load(self.configs.load)
            state_dict = checkpoint["model_state_dict"]

            if "base_shifts" in state_dict:
                state_dict.pop("base_shifts")

            self.logger.info(self.model.load_state_dict(state_dict, strict=False))

            if "output_mask_radius" in checkpoint:
                self.output_mask.update_radius(checkpoint["output_mask_radius"])

            for key in self.optimizers:
                self.optimizers[key].load_state_dict(
                    checkpoint["optimizers_state_dict"][key]
                )

            self.start_epoch = checkpoint["epoch"] + 1
            self.pretrain_iter = 0
        else:
            self.start_epoch = 1

        # save configuration
        self.configs.write(os.path.join(self.outdir, "train-configs.yaml"))

        self.predicted_rots = (
            np.eye(3).reshape(1, 3, 3).repeat(self.image_count, axis=0)
        )
        self.predicted_trans = (
            np.zeros((self.image_count, 2)) if not self.configs.no_trans else None
        )
        self.predicted_conf = (
            np.zeros((self.particle_count, self.configs.z_dim))
            if self.configs.z_dim > 0
            else None
        )

        self.data_iterator = make_dataloader(
            self.data,
            batch_size=self.configs.batch_size,
            shuffler_size=self.configs.shuffler_size,
        )

    def train(self) -> None:
        t0 = dt.now()

        self.pretrain()
        self.current_epoch = self.start_epoch
        self.logger.info("--- Training Starts Now ---")

        self.total_batch_count = 0
        self.total_particles_count = 0
        self.conf_search_particles = 0

        while self.current_epoch <= self.final_epoch:
            will_make_summary = (
                (
                    self.current_epoch == self.final_epoch
                    or self.configs.log_heavy_interval
                    and self.current_epoch % self.configs.log_heavy_interval == 0
                )
                or self.is_in_pose_search_step
                or self.pretraining
            )
            self.log_latents = will_make_summary

            if will_make_summary:
                self.logger.info("Will make a full summary at the end of this epoch")

            te = dt.now()
            self.cur_loss = 0
            self.train_epoch()

            total_loss = self.cur_loss / self.train_particles
            self.logger.info(
                f"# =====> SGD Epoch: {self.epoch} "
                f"finished in {dt.now() - te}; "
                f"total loss = {format(total_loss, '.6f')}"
            )

            # image and pose summary
            if will_make_summary:
                self.make_summary()

            self.current_epoch += 1

        t_total = dt.now() - t0
        self.logger.info(
            f"Finished in {t_total} ({t_total / self.num_epochs} per epoch)"
        )

    def pretrain(self):
        """Pretrain the decoder using random initial poses."""
        particles_seen = 0
        loss = None
        self.logger.info(f"Using random poses for {self.configs.pretrain} iterations")

        for batch in self.data_iterator:
            particles_seen += len(batch[0])

            batch = (
                (batch[0].to(self.device), None)
                if self.configs.tilt is None
                else (batch[0].to(self.device), batch[1].to(self.device))
            )
            loss = self.pretrain_step(batch)

            if particles_seen % self.configs.log_interval == 0:
                self.logger.info(
                    f"[Pretrain Iteration {particles_seen}] loss={loss:4f}"
                )

            if particles_seen > self.configs.pretrain_iter:
                break

        # reset model after pretraining
        if self.configs.reset_optim_after_pretrain:
            self.logger.info(">> Resetting optim after pretrain")
            self.optim = torch.optim.Adam(
                self.model.parameters(),
                lr=self.configs.learning_rate,
                weight_decay=self.configs.weight_decay,
            )

        return loss

    @abstractmethod
    def train_epoch(self):
        pass

    @abstractmethod
    def pretrain_step(self, batch, **pretrain_kwargs):
        pass

    @abstractmethod
    def make_summary(self):
        pass
