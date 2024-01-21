"""Storing and retrieving model parameter settings set by the user and from defaults."""

from datetime import datetime
import argparse
from typing import Optional, Any
import cryodrgn
import os
import sys
import yaml
import numpy as np
from collections import OrderedDict
from abc import ABC
from cryodrgn.utils import load_yaml


class _BaseConfigurations(ABC):
    """Base class for sets of model configuration parameters."""

    # the base parameters for all sets: whether we are testing the installation,
    # where the output is located and which configuration shortcuts were used
    __slots__: tuple[str] = ("outdir", "quick_config")

    # any other parameters belong to this set if and only if they have a default value
    # defined in this dictionary, ordering makes e.g. printing easier for user
    defaults: OrderedDict[str, Any] = OrderedDict()

    @classmethod
    @property
    def parameters(cls) -> list[str]:
        return list(cls.defaults.keys())

    @classmethod
    def parse_args(cls, args):
        return cls(
            {
                par: (getattr(args, par) if hasattr(args, par) else cls.defaults[par])
                for par in tuple(cls.defaults) + cls.__slots__
            }
        )

    def __init__(self, config_vals: dict[str, Any]) -> None:
        if "test_installation" in config_vals and config_vals["test_installation"]:
            print("Installation was successful!")
            sys.exit()

        if "outdir" not in config_vals:
            raise ValueError("`config_vals` must have a `outdir` entry!")
        if "quick_config" not in config_vals:
            raise ValueError("`config_vals` must have a `quick_config` entry!")

        for key in set(config_vals) - {"outdir", "quick_config"}:
            if key not in self.defaults:
                raise ValueError("Unrecognized configuration " f"parameter `{key}`!")

        self.outdir = config_vals["outdir"]
        self.quick_config = config_vals["quick_config"]

        # an attribute is created for every entry in the defaults dictionary
        for key, value in self.defaults.items():
            if key in config_vals:
                setattr(self, key, config_vals[key])

            # if not in given parameters, use defaults
            else:
                setattr(self, key, value)

    def __iter__(self):
        return iter((par, getattr(self, par)) for par in self.parameters)

    def __str__(self):
        return "\n".join([f"{par}{str(val):>20}" for par, val in self])

    def write(self, fl: str) -> None:
        """Saving configurations to file using the original order."""

        with open(fl, "w") as f:
            yaml.dump(dict(self), f, default_flow_style=False, sort_keys=False)


class AbInitioConfigurations(_BaseConfigurations):

    __slots__ = ("particles", "ctf", "pose", "dataset", "log_interval")
    defaults = OrderedDict(
        {
            "particles": None,
            "ctf": None,
            "pose": None,
            "dataset": None,
            "log_interval": 1000,
        }
    )

    def __init__(self, config_vals: dict[str, Any]) -> None:
        super().__init__(config_vals)

        if self.dataset is None:
            if self.particles is None:
                raise ValueError(
                    "As dataset was not specified, please " "specify particles!"
                )

            if self.ctf is None:
                raise ValueError("As dataset was not specified, please " "specify ctf!")


class AmortizedInferenceConfigurations(_BaseConfigurations):

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


def load_configs(outdir: str) -> dict:
    cfg_fl = os.path.join(outdir, "configs.yaml")

    if not os.path.isfile(cfg_fl):
        raise ValueError(
            f"Folder `{outdir}` does not contain a cryoDRGN configuration file!"
        )

    return load_yaml(cfg_fl)


def save(config: dict, filename: Optional[str] = None, folder: Optional[str] = None):
    filename = filename or "config.yaml"
    if folder is not None:
        filename = os.path.join(folder, filename)

    # Add extra useful information to incoming config dict
    if "version" not in config:
        config["version"] = cryodrgn.__version__
    if "time" not in config:
        config["time"] = datetime.now()
    if "cmd" not in config:
        config["cmd"] = sys.argv

    cryodrgn.utils.save_yaml(config, filename)
    return filename


def update_config_v1(config: dict) -> dict:
    arg = "feat_sigma"
    if arg not in config["model_args"]:
        assert config["model_args"]["pe_type"] != "gaussian"
        config["model_args"][arg] = None

    # older version used relu
    config["model_args"].setdefault("activation", "relu")
    config["model_args"].setdefault("tilt_params", {})

    return config


def overwrite_config(config: dict, args: argparse.Namespace) -> dict:
    args_dict = vars(args)

    if hasattr(args, "norm") and args.norm is not None:
        config["dataset_args"]["norm"] = args.norm
    if hasattr(args, "D") and args.D is not None:
        config["lattice_args"]["D"] = args.D + 1
    if hasattr(args, "l_extent") and args.l_extent is not None:
        config["lattice_args"]["extent"] = args.l_extent

    # Overwrite any arguments that are not None
    for arg in (
        "qlayers",
        "qdim",
        "zdim",
        "encode_mode",
        "players",
        "pdim",
        "enc_mask",
        "pe_type",
        "feat_sigma",
        "pe_dim",
        "domain",
        "activation",
    ):
        # Set default to None to maintain backwards compatibility
        if arg in ("pe_dim", "feat_sigma") and arg not in config["model_args"]:
            assert (
                args_dict[arg] is None
            ), f"Should not reach here. Something is wrong: {arg}"
            config["model_args"][arg] = None
            continue

        # Set default activation to ReLU to maintain backwards compatibility with v0.3.1 and earlier
        if arg == "activation" and arg not in config["model_args"]:
            assert (
                args_dict[arg] == "relu"
            ), f"Should not reach here. Something is wrong: {arg}"
            config["model_args"]["activation"] = "relu"
            continue

        if arg in args_dict and args_dict[arg] is not None:
            config["model_args"][arg] = args_dict[arg]

    return config
