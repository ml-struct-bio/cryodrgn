"""Select particle or volume data based on (kmeans) cluster labels.

Example usage
-------------
$ cryodrgn_utils select_clusters labels.pkl --sel 1 2 3 -o selected.pkl
$ cryodrgn_utils select_clusters labels.pkl --sel 1 2 3 -o selected.pkl \
    --parent-ind parent.pkl --N-orig 103547

See also
--------
`cryodrgn_utils select_random` for selecting a random subset of particles

"""
import argparse
import os
import logging
from cryodrgn import analysis, utils

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("labels", help="Input labels.pkl")
    parser.add_argument(
        "--sel", nargs="+", type=int, help="1-indexed cluster IDs to select"
    )
    parser.add_argument(
        "-o", type=os.path.abspath, help="Output particle index selection (.pkl)"
    )

    group = parser.add_argument_group(
        "Get original particle selection (if "
        "trained on a subset of the dataset with --ind)"
    )
    group.add_argument("--parent-ind", type=os.path.abspath, help="Parent index .pkl")
    group.add_argument(
        "--N-orig", type=int, help="Number of particles in original dataset"
    )


def main(args: argparse.Namespace) -> None:
    if any(i < 1 for i in args.sel):
        raise ValueError("Cluster indices must be 1-indexed!")

    labels = utils.load_pkl(args.labels)
    logger.info(f"{len(labels)} particles")
    logger.info(f"Selecting clusters {args.sel}")

    ind = analysis.get_ind_for_cluster(labels, [i - 1 for i in args.sel])
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
