"""Filter a particle stack"""

import argparse
import logging
from cryodrgn import utils
from cryodrgn.mrc import MRCFile
from cryodrgn.source import ImageSource

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument("input", help="Input particles (.mrcs, .txt, .star, .cs)")
    parser.add_argument("--ind", required=True, help="Selected indices array (.pkl)")
    parser.add_argument("-o", help="Output .mrcs file")
    return parser


def main(args):
    ind = utils.load_pkl(args.ind)
    src = ImageSource.from_file(args.input, lazy=True, indices=ind)
    logger.info(f"Loaded {len(src)} particles")
    MRCFile.write(args.o, src)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
