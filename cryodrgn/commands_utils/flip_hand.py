"""Flip handedness of an .mrc file"""

import argparse
import logging
import numpy as np
from cryodrgn.mrc import MRCFile, MRCHeader
from cryodrgn.source import ImageSource

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument("input", help="Input volume (.mrc)")
    parser.add_argument("-o", help="Output volume (.mrc)")
    return parser


def main(args):
    assert args.input.endswith(".mrc"), "Input volume must be .mrc file"
    assert args.o.endswith(".mrc"), "Output volume must be .mrc file"
    header = MRCHeader.parse(args.input)
    data = np.array(ImageSource.from_mrcs(args.input).images())
    assert isinstance(data, np.ndarray)
    data = data[::-1]
    MRCFile.write(args.o, data, header=header)
    logger.info(f"Wrote {args.o}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
