"""Invert the contrast of an .mrc file

Example usage
-------------
# Writes to vol_000_inverted.mrc
$ cryodrgn_utils invert_contrast vol_000.mrc

# Manually specify an output file
$ cryodrgn_utils invert_contrast vol_000.mrc -o outputs/vol-inv.mrc

"""
import os
import argparse
import logging
from cryodrgn.mrcfile import parse_mrc, write_mrc

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input", help="Input volume (.mrc)")
    parser.add_argument(
        "--outmrc", "-o", type=os.path.abspath, help="Output volume (.mrc)"
    )


def main(args: argparse.Namespace) -> None:
    if not args.input.endswith(".mrc"):
        raise ValueError(f"Input volume {args.input} is not a .mrc file!")
    outmrc = args.outmrc or args.input.replace(".mrc", "_inverted.mrc")
    if os.path.dirname(outmrc):
        os.makedirs(os.path.dirname(outmrc), exist_ok=True)
    if not outmrc.endswith(".mrc"):
        raise ValueError(f"Output volume {outmrc} is not a .mrc file!")

    vol, header = parse_mrc(args.input)
    write_mrc(outmrc, -vol, header=header)
    logger.info(f"Wrote {outmrc}")
