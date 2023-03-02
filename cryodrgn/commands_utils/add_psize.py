"""Add pixel size to header of .mrc file"""

import argparse
import logging
from cryodrgn.mrc import MRCFile, MRCHeader
from cryodrgn.source import ImageSource

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument("input", help="Input volume (.mrc)")
    parser.add_argument(
        "--Apix", type=float, default=1, help="Angstrom/pixel (default: %(default)s)"
    )
    parser.add_argument("-o", help="Output volume (.mrc)")
    return parser


def main(args):
    assert args.input.endswith(".mrc"), "Input volume must be .mrc file"
    assert args.o.endswith(".mrc"), "Output volume must be .mrc file"
    header = MRCHeader.parse(args.input)
    header.update_apix(args.Apix)

    src = ImageSource.from_file(args.input)
    MRCFile.write(args.o, src, header=header)

    logger.info(f"Wrote {args.o}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
