"""Select particle or volume data based on (kmeans) cluster labels"""

import argparse
import os
import logging
from cryodrgn import analysis, utils

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument("labels", help="Input labels.pkl")
    parser.add_argument("--sel", nargs="+", type=int, help="Ids of clusters to select")
    parser.add_argument(
        "-o", type=os.path.abspath, help="Output particle index selection (.pkl)"
    )

    group = parser.add_argument_group(
        "Get original particle selection (if trained on a subset of the dataset with --ind)"
    )
    group.add_argument("--parent-ind", type=os.path.abspath, help="Parent index .pkl")
    group.add_argument(
        "--N-orig", type=int, help="Number of particles in original dataset"
    )
    return parser


def main(args):
    labels = utils.load_pkl(args.labels)
    logger.info(f"{len(labels)} particles")
    logger.info(f"Selecting clusters {args.sel}")
    ind = analysis.get_ind_for_cluster(labels, args.sel)
    logger.info(f"Selected {len(ind)} particles")
    logger.info(ind)
    if args.parent_ind is not None:
        logger.info("Converting to original indices")
        parent_ind = utils.load_pkl(args.parent_ind)
        assert args.N_orig
        ind = analysis.convert_original_indices(ind, args.N_orig, parent_ind)
        logger.info(ind)
    utils.save_pkl(ind, args.o)
    logger.info(f"Saved {args.o}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
