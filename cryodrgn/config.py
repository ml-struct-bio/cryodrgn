"""Tools for working with cryoDRGN configuration parameters saved to .yaml files."""

from datetime import datetime
import os.path
import sys
from typing import Optional, Union
import cryodrgn
from cryodrgn import utils


def load(config: Union[str, dict]) -> dict:
    if isinstance(config, str):
        ext = os.path.splitext(config)[-1]
        if ext == ".pkl":
            raise RuntimeError(
                "Loading configuration from a .pkl file is deprecated. Please "
                "save/load configuration as a .yaml file instead."
            )
        elif ext in (".yml", ".yaml"):
            return utils.load_yaml(config)
        else:
            raise RuntimeError(f"Unrecognized config extension {ext}")
    else:
        return config


def save(
    config: dict, filename: Optional[str] = None, folder: Optional[str] = None
) -> str:
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


def update_config_v1(config: Union[str, dict]) -> dict:
    config = load(config)
    arg = "feat_sigma"
    if arg not in config["model_args"]:
        assert config["model_args"]["pe_type"] != "gaussian"
        config["model_args"][arg] = None

    # older version used relu
    config["model_args"].setdefault("activation", "relu")
    config["model_args"].setdefault("tilt_params", {})

    return config
