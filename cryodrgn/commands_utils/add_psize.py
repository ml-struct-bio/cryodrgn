"""Add pixel size to the header of .mrc file containing a volume.

Example usage
-------------
# Overwrite given file with new Apix in header
$ cryodrgn_utils add_psize my_volume.mrc --Apix 1.73 -o my_volume.mrc

"""
import argparse
import logging
from cryodrgn.mrcfile import parse_mrc, write_mrc

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input", help="Input volume (.mrc)")
    parser.add_argument(
        "--Apix", type=float, default=1, help="Angstrom/pixel (default: %(default)s)"
    )
    parser.add_argument("-o", help="Output volume (.mrc)")


def main(args: argparse.Namespace) -> None:
    assert args.input.endswith(".mrc"), "Input volume must be .mrc file"
    assert args.o.endswith(".mrc"), "Output volume must be .mrc file"

    vol, header = parse_mrc(args.input)
    header.apix = args.Apix
    write_mrc(args.o, vol, header=header)
    logger.info(f"Wrote {args.o}")
