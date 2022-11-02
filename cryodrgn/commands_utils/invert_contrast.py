"""Invert the contrast of an .mrc file"""

import argparse

from cryodrgn import mrc, utils

log = utils.log


def add_args(parser):
    parser.add_argument("input", help="Input volume (.mrc)")
    parser.add_argument("-o", help="Output volume (.mrc)")
    return parser


def main(args):
    assert args.input.endswith(".mrc"), "Input volume must be .mrc file"
    assert args.o.endswith(".mrc"), "Output volume must be .mrc file"
    x, h = mrc.parse_mrc(args.input)
    x *= -1
    mrc.write(args.o, x, header=h)
    log(f"Wrote {args.o}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
