"""Flip handedness of an .mrc file

Example usages
--------------
# writes to vol_000_flipped.mrc
$ cryodrgn_utils flip_hand vol_000.mrc

$ cryodrgn_utils flip_hand vol_000.mrc -o vol-flipped.mrc

"""
import os
import argparse
import logging
import numpy as np
from cryodrgn.mrc import MRCFile
from cryodrgn.source import ImageSource

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument("input", help="Input volume (.mrc)")
    parser.add_argument(
        "--outmrc", "-o", type=os.path.abspath, help="Output volume (.mrc)"
    )
    return parser


def main(args):
    if not args.input.endswith(".mrc"):
        raise ValueError(f"Input volume {args.input} is not a .mrc file!")
    outmrc = args.outmrc or args.input.replace(".mrc", "_flipped.mrc")
    if os.path.dirname(outmrc):
        os.makedirs(os.path.dirname(outmrc), exist_ok=True)
    if not outmrc.endswith(".mrc"):
        raise ValueError(f"Output volume {outmrc} is not a .mrc file!")

    src = ImageSource.from_file(args.input)
    # Note: Proper flipping (compatible with legacy implementation) only happens when chunksize is equal to src.n
    MRCFile.write(
        outmrc,
        src,
        transform_fn=lambda data, indices: np.array(data.cpu())[::-1],
        chunksize=src.n,
    )
    logger.info(f"Wrote {outmrc}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
