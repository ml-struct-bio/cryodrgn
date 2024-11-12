"""Storing and retrieving model parameter settings set by the user and from defaults."""

from datetime import datetime
import argparse
from typing import Optional
import cryodrgn.utils
import os
import sys


def load_configs(outdir: str) -> dict:
    cfg_fl = os.path.join(outdir, "configs.yaml")

    if not os.path.isfile(cfg_fl):
        raise ValueError(
            f"Folder `{outdir}` does not contain a cryoDRGN configuration file!"
        )

    return cryodrgn.utils.load_yaml(cfg_fl)


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
