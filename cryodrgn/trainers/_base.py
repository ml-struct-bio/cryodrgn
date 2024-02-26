"""Base classes for model training engines."""
import os
import shutil
import sys
from collections import OrderedDict
import numpy as np
from abc import ABC, abstractmethod
from typing import Any
from typing_extensions import Self
import yaml
from datetime import datetime as dt
import logging
from cryodrgn import ctf
from cryodrgn import dataset
from cryodrgn import utils
from cryodrgn.lattice import Lattice
import torch
from torch import nn


class BaseTrainer(ABC):
    """Abstract base class for reconstruction model training engines.

    Attributes
    ----------
    outdir(str):    Path to where experiment output is stored.
    quick_config(dict):    Configuration shortcuts.
    """

    # a configuration parameter is defined for this engine if and only if it has a
    # default value defined here; ordering makes e.g. printing easier for user
    config_defaults: OrderedDict[str, Any] = OrderedDict()

    @classmethod
    @property
    def parameters(cls) -> list[str]:
        return list(cls.config_defaults.keys())

    @classmethod
    def parse_args(cls, args) -> Self:
        return cls(
            {
                par: (
                    getattr(args, par)
                    if hasattr(args, par)
                    else cls.config_defaults[par]
                )
                for par in tuple(cls.config_defaults)
            }
        )

    @classmethod
    def parse_configs(cls, configs: dict[str, Any]) -> dict[str, Any]:
        if "test_installation" in configs and configs["test_installation"]:
            print("Installation was successful!")
            sys.exit()

        if "outdir" in configs and configs["outdir"] is not None:
            configs["outdir"] = os.path.abspath(configs["outdir"])
        else:
            configs["outdir"] = os.getcwd()

        # process the quick_config parameter
        if "quick_config" in configs and configs["quick_config"] is not None:
            for key, value in configs["quick_config"].items():
                if key not in cls.quick_config:
                    raise ValueError(
                        "Unrecognized parameter " f"shortcut field `{key}`!"
                    )

                if value not in cls.quick_config[key]:
                    raise ValueError(
                        "Unrecognized parameter shortcut label "
                        f"`{value}` for field `{key}`!"
                    )

                for _key, _value in cls.quick_config[key][value].items():
                    if _key not in cls.config_defaults:
                        raise ValueError(
                            "Unrecognized configuration " f"parameter `{key}`!"
                        )

                    # given parameters have priority
                    if _key not in configs:
                        configs[_key] = _value

        return configs

    def __init__(self, configs: dict[str, Any]) -> None:
        configs = self.parse_configs(configs)
        self.logger = logging.getLogger(__name__)
        self.outdir = configs["outdir"]

        # output directory
        if os.path.exists(self.outdir):
            self.logger.warning("Output folder already exists, renaming the old one.")

            newdir = self.outdir + "_old"
            if os.path.exists(newdir):
                self.logger.warning("Must delete the previously saved old output.")
                shutil.rmtree(newdir)

            os.rename(self.outdir, newdir)

        os.makedirs(self.outdir, exist_ok=True)

    def __iter__(self):
        return iter((par, getattr(self, par)) for par in self.parameters)

    def __str__(self):
        return "\n".join([f"{par}{str(val):>20}" for par, val in self])

    def write(self, fl: str) -> None:
        """Saving configurations to file using the original order."""

        with open(fl, "w") as f:
            yaml.dump(dict(self), f, default_flow_style=False, sort_keys=False)


class ModelTrainer(BaseTrainer, ABC):

    config_defaults = OrderedDict({"model": "amort", "verbose": 0, "seed": -1})

    quick_config = {
        "capture_setup": {
            "spa": dict(),
            "et": {
                "subtomogram_averaging": True,
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

    # options for optimizers to use
    optim_types = {"adam": torch.optim.Adam, "lbfgs": torch.optim.LBFGS}

    @classmethod
    def parse_configs(cls, configs: dict[str, Any]) -> dict[str, Any]:
        if "model" in configs:
            if configs["model"] not in {"hps", "amort"}:
                raise ValueError(
                    f"Given model `{configs['model']}` not in currently supported "
                    f"model types:\n"
                    f"`amort` (cryoDRGN v4 amortized inference pose estimation)\n"
                    f"`hps` (cryoDRGN v3 hierarchical pose estimation)\n"
                )
        else:
            configs["model"] = cls.config_defaults["model"]

        if "verbose" in configs:
            if not isinstance(configs["verbose"], int) or configs["verbose"] < 0:
                raise ValueError(
                    f"Given verbosity `{configs['verbose']}` not a positive integer!"
                )
        else:
            configs["verbose"] = cls.config_defaults["verbose"]

        if "seed" in configs and not isinstance(configs["seed"], int):
            raise ValueError(
                "Configuration `seed` must be given as an integer, "
                f"given `{configs['seed']}` instead!"
            )
        else:
            configs["seed"] = cls.config_defaults["seed"]

        if "ind" in configs:
            if isinstance(configs["ind"], str) and not os.path.exists(configs["ind"]):
                raise ValueError(
                    f"Given configuration subset file {configs['ind']} does not exist!"
                )

        if configs["dataset"]:
            paths_file = os.environ.get("DRGNAI_DATASETS")
            paths = utils.load_yaml(paths_file)

            if configs["dataset"] not in paths:
                raise ValueError(
                    f"Given dataset label `{configs['dataset']}` not in list of "
                    f"datasets specified at `{paths_file}`!"
                )
            use_paths = paths[configs["dataset"]]

            configs["particles"] = use_paths["particles"]
            for k in ["ctf", "pose", "datadir", "labels", "ind", "dose_per_tilt"]:
                if k not in configs and k in use_paths:
                    if not os.path.exists(use_paths[k]):
                        raise ValueError(
                            f"Given {k} file `{use_paths[k]}` does not exist!"
                        )
                    configs[k] = use_paths[k]

        return super().parse_configs(configs)

    @abstractmethod
    def make_model(self, configs: dict[str, Any]) -> nn.Module:
        pass

    def __init__(self, configs: dict[str, Any]) -> None:
        if "load" in configs and configs["load"] == "latest":
            configs = self.get_latest_configs()
            self.logger.info(str(configs))

        super().__init__(configs)

        self.verbose = configs["verbose"]
        if self.verbose:
            self.logger.setLevel(logging.DEBUG)

        self.logger.addHandler(
            logging.FileHandler(os.path.join(self.outdir, "training.log"))
        )

        # set the random seed
        if configs["seed"] >= 0:
            self.seed = configs["seed"]
        else:
            self.seed = np.random.randint(0, 10000)

        np.random.seed(self.seed)
        torch.manual_seed(self.seed)

        # set the device
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
        if configs["ind"] is not None:
            if isinstance(configs["ind"], int):
                self.logger.info(f"Keeping the first {configs['ind']} particles")
                self.ind = np.arange(configs["ind"])
            else:
                self.logger.info(f"Filtering image dataset with {configs['ind']}")
                self.ind = utils.load_yaml(configs["ind"])
        else:
            self.ind = None

        # load dataset
        self.logger.info(f"Loading dataset from {configs['particles']}")

        if "tilt" in configs and configs["tilt"] is not None:
            self.tilt = torch.tensor(
                utils.xrot(configs["tilt_deg"]).astype(np.float32), device=self.device
            )
        else:
            self.tilt = None

        if "subtomogram_averaging" in configs and configs["subtomogram_averaging"]:
            self.subtomo_averaging = True
        else:
            self.subtomo_averaging = False

        if self.subtomo_averaging:
            self.data = dataset.ImageDataset(
                mrcfile=configs["particles"],
                norm=configs["data_norm"],
                invert_data=configs["invert_data"],
                ind=configs["ind"],
                window=configs["window"],
                keepreal=configs["use_real"],
                datadir=configs["datadir"],
                window_r=configs["window_r"],
            )
        else:
            self.data = dataset.TiltSeriesData(
                tiltstar=configs["particles"],
                ind=configs["ind"],
                ntilts=configs["n_tilts"],
                angle_per_tilt=configs["angle_per_tilt"],
                window_r=configs["window_radius_gt_real"],
                datadir=configs["datadir"],
                max_threads=configs["max_threads"],
                dose_per_tilt=configs["dose_per_tilt"],
                device=self.device,
                poses_gt_pkl=configs["pose"],
                tilt_axis_angle=configs["tilt_axis_angle"],
                no_trans=configs["no_trans"],
            )

        self.Nimg = self.data.N
        self.resolution = self.data.D

        if configs["ctf"]:
            self.logger.info(f"Loading ctf params from {configs['ctf']}")
            ctf_params = ctf.load_ctf_for_training(self.resolution - 1, configs["ctf"])

            if self.ind is not None:
                ctf_params = ctf_params[self.ind]
            assert ctf_params.shape == (self.Nimg, 8), ctf_params.shape

            if self.subtomo_averaging:
                ctf_params = np.concatenate(
                    (ctf_params, self.data.ctfscalefactor.reshape(-1, 1)),
                    axis=1,  # type: ignore
                )
                self.data.voltage = float(ctf_params[0, 4])

            self.ctf_params = torch.tensor(ctf_params, device=self.device)
            if self.subtomo_averaging:
                self.data.voltage = float(self.ctf_params[0, 4])

        else:
            self.ctf_params = None

        # lattice
        self.logger.info("Building lattice...")
        self.lattice = Lattice(self.resolution, extent=0.5, device=self.device)

        self.logger.info("Initializing model...")
        self.model = self.make_model(configs)
        self.model.to(self.device)
        self.logger.info(self.model)
        param_count = sum(p.numel() for p in self.model.parameters() if p.requires_grad)
        self.logger.info(f"{param_count} parameters in model")
        self.current_epoch = None

        # TODO: auto-loading from last weights file if load=True?
        # initialization from a previous checkpoint
        if configs["load"]:
            self.logger.info(f"Loading checkpoint from {configs['load']}")
            checkpoint = torch.load(configs["load"])
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
        configs.write(os.path.join(self.outdir, "train-configs.yaml"))

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
    def pretrain_step(self):
        pass

    @abstractmethod
    def make_summary(self):
        pass
