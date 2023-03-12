"""
Display config information of a cryoDRGN job
"""

import argparse
import os
import os.path
from pprint import pprint
import logging
import warnings
from cryodrgn import utils

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument(
        "workdir", type=os.path.abspath, help="Directory with cryoDRGN results"
    )
    return parser


def main(args):
    warnings.warn(
        "The view_config command is deprecated."
        "Please save configuration in yaml format and view the config.yaml file directly.",
        DeprecationWarning,
    )
    config_yaml = f"{args.workdir}/config.yaml"
    config_pkl = f"{args.workdir}/config.pkl"
    if os.path.exists(config_yaml):
        cfg = utils.load_yaml(config_yaml)
    elif os.path.exists(config_pkl):
        cfg = utils.load_pkl(config_pkl)
    else:
        raise RuntimeError(f"A configuration file was not found at {args.workdir}")

    if "version" in cfg:
        logger.info(f'Version: {cfg["version"]}')
    if "time" in cfg:
        logger.info(f'Creation time: {cfg["time"]}')
    if "cmd" in cfg:
        logger.info("Command:")
        print(" ".join(cfg["cmd"]))
    logger.info("Config:")
    pprint(cfg)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
