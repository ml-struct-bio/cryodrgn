"""Filter a particle stack"""

import argparse
import logging
from cryodrgn import utils
from cryodrgn.source import ImageSource, write_mrc

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser):
    parser.add_argument("input", help="Input particles (.mrcs, .txt, .star, .cs)")
    parser.add_argument("--ind", required=True, help="Selected indices array (.pkl)")
    parser.add_argument("-o", help="Output .mrcs file")
    return parser


def main(args: argparse.Namespace):
    ind = utils.load_pkl(args.ind)
    src = ImageSource.from_file(args.input, lazy=True, indices=ind)
    logger.info(f"Loaded {len(src)} particles")
    write_mrc(args.o, src)
