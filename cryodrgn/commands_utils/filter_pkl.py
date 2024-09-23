"""Filter cryoDRGN data stored in a .pkl file, writing to a new .pkl file.

Example usage
-------------
$ cryodrgn_utils filter_pkl pose.pkl --ind my_indices.pkl -o filtered_pose.pkl

# Use first <n> data points instead of a set of filtering indices
$ cryodrgn_utils filter_pkl ctf.pkl --first 10000 -o ctf_first-10k.pkl

"""
import argparse
import os
import logging
import numpy as np
from cryodrgn.utils import load_pkl, save_pkl

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    """Specifies the command-line interface used by `cryodrgn_utils filter_pkl`."""

    parser.add_argument("input", help="Input data file (.pkl)")
    parser.add_argument("--ind", help="Array of selected indices (.pkl)")
    parser.add_argument(
        "--first", type=int, help="Alternatively, save the first N datapoints"
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=os.path.abspath,
        help="Output data file (.pkl)",
    )


def main(args: argparse.Namespace) -> None:
    """Running the command `cryodrgn_utils filter_pkl` (see `add_args` above)."""

    if not (args.ind is None) ^ (args.first is None):
        raise ValueError(
            "Must provide exactly one of `--ind` (for filtering using a set of indices)"
            " or `--first` (for filtering out everything but the first <x> datapoints)"
        )
    indices = load_pkl(args.ind) if args.ind is not None else np.arange(args.first)
    input_data = load_pkl(args.input)

    # Pickled pose files (e.g. pose.pkl) can contain a tuple of rotations, translations
    if isinstance(input_data, tuple):
        logger.info("Detected pose.pkl")
        logger.info(f"Old shape: {[xx.shape for xx in input_data]}")
        input_data = tuple(x[indices] for x in input_data)
        logger.info(f"New shape: {[x.shape for x in input_data]}")

    # Handles all other pickled files cryoDRGN can encounter
    else:
        logger.info(f"Old shape: {input_data.shape}")
        input_data = input_data[indices]
        logger.info(f"New shape: {input_data.shape}")

    logger.info(f"Saving {args.output}")
    save_pkl(input_data, args.output)
