"""Filter a particle stack using given indices to produce a new stack file.

Example usage
-------------
cryodrgn_utils filter_mrcs particles.mrcs --ind good-particles.pkl -o new-particles.mrcs

# Will save output to particles_indices004.mrcs
cryodrgn_utils filter_mrcs particles.mrcs --ind indices004.pkl

"""
import os
import argparse
import logging
from cryodrgn import utils
from cryodrgn.source import ImageSource

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input", help="Input particles (.mrcs, .txt, .star, .cs)")
    parser.add_argument("--ind", required=True, help="Selected indices array (.pkl)")
    parser.add_argument("--outfile", "-o", help="Output .mrcs file")


def main(args: argparse.Namespace) -> None:
    ind = utils.load_pkl(args.ind)
    src = ImageSource.from_file(args.input, lazy=True, indices=ind)
    logger.info(
        f"Loaded {src.orig_n} particles which were "
        f"filtered down to {len(src)} particles."
    )

    if args.outfile is None:
        outfile = (
            f"{os.path.splitext(args.input)[0]}"
            f"_{os.path.basename(os.path.splitext(args.ind)[0])}.mrcs"
        )
    else:
        outfile = str(args.outfile)

    src.write(outfile)
