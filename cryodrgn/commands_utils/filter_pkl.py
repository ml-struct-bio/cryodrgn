"""Filter cryodrgn data stored in a .pkl file"""

import argparse
import os
import pickle
import logging
import numpy as np

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument("input", help="Input data (.pkl)")
    parser.add_argument("--ind", help="Array of selected indices (.pkl)")
    parser.add_argument(
        "--first", type=int, help="Alternatively, save the first N datapoints"
    )
    parser.add_argument("-o", type=os.path.abspath, help="Output data (.pkl)")
    return parser


def load_pkl(x):
    return pickle.load(open(x, "rb"))


def main(args):
    x = load_pkl(args.input)
    if args.first:
        assert args.ind is None
        ind = np.arange(args.first)
    else:
        assert args.first is None
        ind = load_pkl(args.ind)

    # pose.pkl contains a tuple of rotations and translations
    if type(x) == tuple:
        logger.info("Detected pose.pkl")
        logger.info(f"Old shape: {[xx.shape for xx in x]}")
        x = (xx[ind] for xx in x)
        x = tuple(x)
        logger.info(f"New shape: {[xx.shape for xx in x]}")

    # all other cryodrgn pkls
    else:
        logger.info(f"Old shape: {x.shape}")
        x = x[ind]
        logger.info(f"New shape: {x.shape}")
    logger.info(f"Saving {args.o}")
    pickle.dump(x, open(args.o, "wb"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
