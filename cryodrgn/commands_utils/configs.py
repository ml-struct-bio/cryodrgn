"""Helper functions for loading and manipulating saved model parameter settings."""

from pathlib import Path
import yaml
import pickle

NORMAL_EXCEPTIONS = (
    EOFError,
    BufferError,
    pickle.UnpicklingError,
    ImportError,
    IndexError,
    AttributeError,
)


def check_open_config(d: Path, config: str, version: str) -> tuple[dict, str]:
    version_code = str(version)

    if version == "2":
        try:
            with open(Path(d, config), "rb") as f:
                cfg = pickle.load(f)
        except NORMAL_EXCEPTIONS:
            cfg = dict()
            version_code = "2C"

    else:
        try:
            with open(Path(d, config), "r") as f:
                cfg = yaml.safe_load(f)

        except NORMAL_EXCEPTIONS:
            cfg = dict()
            version_code = f"{version}C"

    if version_code[-1] != "C":
        if version == "4":
            if "particles" not in cfg and "dataset" not in cfg:
                version_code = "N"

        else:
            if (
                "cmd" not in cfg
                or "cryodrgn" not in cfg["cmd"][0]
                or cfg["cmd"][1] not in {"train_nn", "train_vae"}
                or "dataset_args" not in cfg
            ):
                version_code = "N"

    return cfg, version_code
