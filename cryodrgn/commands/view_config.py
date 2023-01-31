"""
Display config information of a cryoDRGN job
"""

import argparse
import os
import pickle
from pprint import pprint
import logging

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument(
        "workdir", type=os.path.abspath, help="Directory with cryoDRGN results"
    )
    return parser


def main(args):
    f = open(f"{args.workdir}/config.pkl", "rb")
    cfg = pickle.load(f)
    try:
        meta = pickle.load(f)
        logger.info(f'Version: {meta["version"]}')
        logger.info(f'Creation time: {meta["time"]}')
        logger.info("Command:")
        print(" ".join(meta["cmd"]))
    except:  # noqa: E722
        pass
    logger.info("Config:")
    pprint(cfg)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
