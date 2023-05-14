"""View images in a particle stack"""

import argparse
import os
import logging
import matplotlib.pyplot as plt
from cryodrgn import analysis, dataset, mrc, utils

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
    stack = dataset.load_particles(args.mrcs, lazy=True, datadir=args.datadir)

    if args.ind is not None:
        logger.info(f"Filtering image dataset with {args.ind}")
        ind = utils.load_pkl(args.ind).astype(int)
        stack = [stack[i] for i in ind]

    logger.info("{} {}x{} images".format(len(stack), *stack[0].get().shape))
    stack = [stack[x].get() for x in range(25)]
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
