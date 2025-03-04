"""Storing and retrieving model parameter settings set by the user and from defaults."""

import os
import sys
from datetime import datetime
import argparse
from typing import Optional
import warnings
import cryodrgn.utils


def load(config):
    if isinstance(config, str):
        ext = os.path.splitext(config)[-1]
        if ext == ".pkl":
            warnings.warn(
                "Loading configuration from a .pkl file is deprecated. Please "
                "save/load configuration as a .yaml file instead."
            )
            return cryodrgn.utils.load_pkl(config)
        elif ext in (".yml", ".yaml"):
            return cryodrgn.utils.load_yaml(config)
        else:
            raise RuntimeError(f"Unrecognized config extension {ext}")
    else:
        return config


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
    config = load(config)
    model_args = config["model_args"] if "model_args" in config else config
    dataset_args = config["dataset_args"] if "dataset_args" in config else config
    lattice_args = config["lattice_args"] if "lattice_args" in config else config

    if hasattr(args, "norm") and args.norm is not None:
        dataset_args["norm"] = args.norm
    if hasattr(args, "D") and args.D is not None:
        lattice_args["D"] = args.D + 1
    if hasattr(args, "l_extent") and args.l_extent is not None:
        lattice_args["extent"] = args.l_extent

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
        if arg in ("pe_dim", "feat_sigma") and arg not in model_args:
            assert (
                not hasattr(args, arg) or getattr(args, arg) is None
            ), f"Should not reach here. Something is wrong: {arg}"
            model_args[arg] = None
            continue

        # Set default activation to ReLU to maintain backwards compatibility
        # with v0.3.1 and earlier
        if arg == "activation" and arg not in model_args:
            assert (
                not hasattr(args, arg) or getattr(args, arg) == "relu"
            ), f"Should not reach here. Something is wrong: {arg}"
            model_args["activation"] = "relu"
            continue

        if hasattr(args, arg) and getattr(args, arg) is not None:
            model_args[arg] = getattr(args, arg)

    if "zdim" in model_args:
        model_args["z_dim"] = model_args["zdim"]

    return config
