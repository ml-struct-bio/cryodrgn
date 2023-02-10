"""Filter a particle stack"""

import argparse
import logging
import numpy as np
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
    x = ImageSource.from_file(args.input, lazy=True)
    logger.info(f"Loaded {len(x)} particles")
    ind = utils.load_pkl(args.ind)
    x = np.array([x.images(i) for i in ind])
    logger.info(f"New dimensions: {x.shape}")
    MRCFile.write(args.o, x)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
