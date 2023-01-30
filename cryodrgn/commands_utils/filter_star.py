"""
Filter a .star file
"""

import argparse
import os
import warnings
import logging
from cryodrgn.commands_utils import write_star

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument("input", help="Input .star file")
    parser.add_argument("--ind", required=True, help="Array of selected indices (.pkl)")
    parser.add_argument("-o", type=os.path.abspath, help="Output .star file")
    return parser


def main(args):
    warning_msg = "cryodrgn_utils filter_star is deprecated. Please use cryodrgn_utils write_star instead."
    warnings.warn(warning_msg, DeprecationWarning)
    logger.warning(f"WARNING: {warning_msg}")

    args = write_star.add_args(argparse.ArgumentParser()).parse_args(
        [
            args.input,
            "-o",
            args.o,
            "--ind",
            args.ind,
        ]
    )
    write_star.main(args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
