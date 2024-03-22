"""Automated generation of masking filters for 3D volumes.

Example usages
--------------
$ cryodrgn_utils gen_mask volume16.mrc mask16.mrc
$ cryodrgn_utils gen_mask volume16.mrc mask16.mrc -p slices.png
$ cryodrgn_utils gen_mask volume16.mrc mask16.mrc --dist 15

"""
import os
import argparse
import logging
import numpy as np
from scipy.ndimage import distance_transform_edt, binary_dilation
from cryodrgn.mrc import MRCFile
from cryodrgn.commands.analyze_landscape import view_slices

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
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
        default=25,
        help="Dilate initial mask by this amount (default: %(default)s angstroms)",
    )
    parser.add_argument(
        "--dist",
        type=int,
        default=15,
        help="Width of cosine edge (default: %(default)s angstroms)",
    )
    parser.add_argument(
        "--Apix",
        type=float,
        help="use this A/px value instead of inferring from the input header",
    )
    parser.add_argument(
        "--png-output",
        "-p",
        type=os.path.abspath,
        help="Output .png file for slice plots",
    )


def main(args: argparse.Namespace) -> None:
    vol, header = MRCFile.parse(args.input)
    apix = args.Apix or header.get_apix()
    thresh = np.percentile(vol, 99.99) / 2 if args.thresh is None else args.thresh
    logger.info(f"A/px={apix:.5g}; Threshold={thresh:.5g}")

    x = (vol >= thresh).astype(bool)

    if args.dilate:
        dilate_val = int(args.dilate // apix)
        logger.info(f"Dilate initial mask by: {dilate_val} px")
        x = binary_dilation(x, iterations=dilate_val)
    else:
        logger.info("no mask dilation applied")

    dist_val = args.dist / apix
    logger.info(f"Width of cosine edge: {dist_val} px")
    if dist_val:
        y = distance_transform_edt(~x.astype(bool))
        y[y > dist_val] = dist_val
        z = np.cos(np.pi * y / dist_val / 2)
    else:
        z = x

    MRCFile.write(args.output, z.astype(np.float32), header=header)
    if args.png_output:
        view_slices(z, out_png=args.png_output)
