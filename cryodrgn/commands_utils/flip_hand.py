"""Flip handedness of an .mrc file

Example usage
-------------
# Writes to vol_000_flipped.mrc
$ cryodrgn_utils flip_hand vol_000.mrc

# You can also specify an output file manually
$ cryodrgn_utils flip_hand vol_000.mrc -o vol-flipped.mrc

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
    outmrc = args.outmrc or args.input.replace(".mrc", "_flipped.mrc")
    if os.path.dirname(outmrc):
        os.makedirs(os.path.dirname(outmrc), exist_ok=True)
    if not outmrc.endswith(".mrc"):
        raise ValueError(f"Output volume {outmrc} is not a .mrc file!")

    vol, header = parse_mrc(args.input)
    vol = vol[::-1]
    write_mrc(outmrc, vol, header=header)
    logger.info(f"Wrote {outmrc}")
