"""Concatenate arrays from multiple .pkl files"""

import argparse
import logging
import numpy as np
from cryodrgn.utils import load_pkl, save_pkl

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input", nargs="+", help="Input .pkl files")
    parser.add_argument("-o", required=True, help="Output .pkl file")


def main(args: argparse.Namespace) -> None:
    x = [load_pkl(f) for f in args.input]

    if any(isinstance(i, np.ndarray) for i in x) and not all(
        isinstance(i, np.ndarray) for i in x
    ):
        raise ValueError(
            "All input files must contain numpy arrays of the same shape or pose tuples"
        )

    if isinstance(x[0], tuple):  # pose tuples
        r = [xx[0] for xx in x]
        t = [xx[1] for xx in x]
        r2 = np.concatenate(r)
        t2 = np.concatenate(t)
        logger.info(r2.shape)
        logger.info(t2.shape)
        x2 = (r2, t2)
    else:
        for i in x:
            logger.info(i.shape)
        x2 = np.concatenate(x)
        logger.info(x2.shape)

    save_pkl(x2, args.o)
