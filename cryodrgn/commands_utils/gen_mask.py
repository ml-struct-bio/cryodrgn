"""Automated generation of masking filters for reconstructed volumes.

Example usages
--------------


"""
import os
import argparse
import logging
import numpy as np
from scipy.ndimage.morphology import distance_transform_edt
from scipy.ndimage.morphology import binary_dilation
from cryodrgn.mrc import MRCFile
from cryodrgn.commands.analyze_landscape import view_slices

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument(
        "input", type=os.path.abspath, help="Input .mrc file for volume to be masked"
    )
    parser.add_argument(
        "output", type=os.path.abspath, help="Output .mrc file for generated mask"
    )

    parser.add_argument(
        "--thresh",
        type=float,
        help="Density value to threshold for masking (default: half of max density value)",
    )
    parser.add_argument(
        "--dilate",
        type=int,
        default=3,
        help="Dilate initial mask by this amount (default: %(default)s pixels)",
    )
    parser.add_argument(
        "--dist",
        type=int,
        default=10,
        help="Width of cosine edge (default: %(default)s pix)",
    )
    parser.add_argument(
        "--png-output",
        "-p",
        type=os.path.abspath,
        help="Output .png file for slice plots",
    )

    return parser


def main(args: argparse.Namespace) -> None:
    vol, header = MRCFile.parse(args.input)
    thresh = np.percentile(vol, 99.99) / 2 if args.thresh is None else args.thresh
    logger.info(f"Threshold: {thresh}")
    x = (vol >= thresh).astype(bool)

    if args.dilate:
        logger.info(f"Dilate initial mask by: {args.dilate}")
        x = binary_dilation(x, iterations=args.dilate)

    logger.info(f"Width of cosine edge: {args.dist}")
    if args.dist:
        y = distance_transform_edt(~x.astype(bool))
        y[y > args.dist] = args.dist
        z = np.cos(np.pi * y / args.dist / 2)
    else:
        z = x

    MRCFile.write(args.output, z.astype(np.float32), header=header)
    if args.png_output:
        view_slices(z, out_png=args.png_output, D=vol.shape[0])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    main(add_args(parser).parse_args())
