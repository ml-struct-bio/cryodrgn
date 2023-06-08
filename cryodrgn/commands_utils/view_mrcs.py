"""View images in a particle stack"""

import argparse
import os
import logging
import matplotlib.pyplot as plt
import os.path
from cryodrgn import analysis, utils
from cryodrgn.source import ImageSource

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument(
        "mrcs", help="Input particles or volume (.mrc, .mrcs, .star, .cs, or .txt)"
    )
    parser.add_argument(
        "--datadir",
        help="Optionally provide path to input .mrcs if loading from a .star or .cs file",
    )
    parser.add_argument("--invert", action="store_true")
    parser.add_argument(
        "--ind",
        type=os.path.abspath,
        metavar="PKL",
        help="Filter particle stack by these indices",
    )
    parser.add_argument("-o", help="Optionally, specify image to save")
    return parser


def main(args):
    ind = None
    if args.ind is not None:
        logger.info(f"Filtering image dataset with {args.ind}")
        ind = utils.load_pkl(args.ind).astype(int)

    src = ImageSource.from_file(args.mrcs, datadir=args.datadir, indices=ind)
    logger.info("{n} {L}x{L} images".format(n=len(src), L=src.images(0).shape[-1]))
    stack = [src.images(x).squeeze(dim=0) for x in range(25)]
    if args.invert:
        stack = [-1 * x for x in stack]
    analysis.plot_projections(stack)
    if args.o:
        plt.savefig(args.o)
        logger.info(f"Wrote {args.o}")
    else:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
