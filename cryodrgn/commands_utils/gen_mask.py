"""Creating masking filters for 3D volumes using threshold dilation with a cosine edge.

Example usage
-------------
$ cryodrgn_utils gen_mask volume16.mrc mask16.mrc
$ cryodrgn_utils gen_mask volume16.mrc mask16.mrc -p slices.png
$ cryodrgn_utils gen_mask volume16.mrc mask16.mrc --dist 15

"""
import os
import argparse
import logging
import numpy as np
import torch
from cryodrgn.masking import cosine_dilation_mask
from cryodrgn.mrcfile import write_mrc, parse_mrc
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
        "--threshold",
        type=float,
        help="Density value to use as the initial threshold for masking "
        "(default: half of max density value)",
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
        help="Use this A/px value instead of inferring from the input header",
    )
    parser.add_argument(
        "--png-output",
        "-p",
        type=os.path.abspath,
        help="Output .png file for slice plots",
    )


def main(args: argparse.Namespace) -> None:
    vol, header = parse_mrc(args.input)
    vol = torch.Tensor(vol)
    apix = args.Apix or header.apix

    mask = cosine_dilation_mask(
        vol,
        threshold=args.threshold,
        dilation=args.dilate,
        edge_dist=args.dist,
        apix=apix,
    )
    header.apix = apix
    write_mrc(args.output, mask.astype(np.float32), header=header)

    if args.png_output:
        view_slices(mask, out_png=args.png_output)
