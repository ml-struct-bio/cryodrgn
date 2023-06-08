from datetime import datetime
import os.path
import sys
from typing import Optional
import warnings
import cryodrgn
from cryodrgn import utils


def load(config):
    if isinstance(config, str):
        ext = os.path.splitext(config)[-1]
        if ext == ".pkl":
            warnings.warn(
                "Loading configuration from a .pkl file is deprecated. Please save/load configuration"
                "as a .yaml file instead."
            )
            return utils.load_pkl(config)
        elif ext in (".yml", ".yaml"):
            return utils.load_yaml(config)
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

    utils.save_yaml(config, filename)
    return filename


def update_config_v1(config):
    config = load(config)
    arg = "feat_sigma"
    if arg not in config["model_args"]:
        assert config["model_args"]["pe_type"] != "gaussian"
        config["model_args"][arg] = None

    # older version used relu
    config["model_args"].setdefault("activation", "relu")
    config["model_args"].setdefault("tilt_params", {})
    return config


def overwrite_config(config, args):
    config = load(config)
    if args.norm is not None:
        config["dataset_args"]["norm"] = args.norm
    v = vars(args)
    if "D" in v and args.D is not None:
        config["lattice_args"]["D"] = args.D + 1
    if "l_extent" in v and args.l_extent is not None:
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
            assert v[arg] is None, f"Should not reach here. Something is wrong: {arg}"
            config["model_args"][arg] = None
            continue
        # Set default activation to ReLU to maintain backwards compatibility with v0.3.1 and earlier
        if arg == "activation" and arg not in config["model_args"]:
            assert v[arg] == "relu"
            config["model_args"]["activation"] = "relu"
            continue
        if v[arg] is not None:
            config["model_args"][arg] = v[arg]
    return config
