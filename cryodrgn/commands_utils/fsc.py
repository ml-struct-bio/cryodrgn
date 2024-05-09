"""Compute Fourier shell correlation between two volumes.

Example usages
--------------
$ cryodrgn_utils fsc volume1.mrc volume2.mrc
$ cryodrgn_utils fsc vol1.mrc vol2.mrc -o fsc.txt -p
$ cryodrgn_utils fsc vol1.mrc vol2.mrc --mask test-mask.mrc -o fsc.txt -p fsc-plot.png

"""
import os
import argparse
import logging
import numpy as np
import pandas as pd
import torch
from typing import Optional
from cryodrgn import fft
from cryodrgn.source import ImageSource
from cryodrgn.commands_utils.plot_fsc import create_fsc_plot

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("volumes", nargs=2, help="volumes to compare")

    parser.add_argument(
        "--mask",
        type=os.path.abspath,
        help="if given, apply the mask in this file before calculating FSCs",
    )
    parser.add_argument(
        "--plot",
        "-p",
        type=os.path.abspath,
        nargs="?",
        const=True,
        default=None,
        help="also plot the FSC curve: optionally supply a .png file name instead of "
        "generating one automatically",
    )
    parser.add_argument(
        "--Apix",
        type=float,
        default=1.0,
        help="Ang/pixels to use when printing the resolutions at thresholds",
    )
    parser.add_argument(
        "--outtxt",
        "-o",
        type=os.path.abspath,
        help=(
            "if given, a file to save the FSC values, "
            "with each space-delimited row as <resolution> <fsc_val>; "
            "otherwise print these to screen"
        ),
    )


def calculate_fsc(
    vol1: torch.Tensor, vol2: torch.Tensor, mask_file: Optional[str] = None
) -> pd.DataFrame:

    if mask_file:
        mask = ImageSource.from_file(mask_file)
        mask = mask.images()
        vol1 *= mask
        vol2 *= mask

    D = vol1.shape[0]
    x = np.arange(-D // 2, D // 2)
    x2, x1, x0 = np.meshgrid(x, x, x, indexing="ij")
    coords = np.stack((x0, x1, x2), -1)
    r = (coords**2).sum(-1) ** 0.5

    assert r[D // 2, D // 2, D // 2] == 0.0
    vol1 = fft.fftn_center(vol1)
    vol2 = fft.fftn_center(vol2)

    # logger.info(r[D//2, D//2, D//2:])
    prev_mask = np.zeros((D, D, D), dtype=bool)
    fsc = [1.0]
    for i in range(1, D // 2):
        mask = r < i
        shell = np.where(mask & np.logical_not(prev_mask))
        v1 = vol1[shell]
        v2 = vol2[shell]
        p = np.vdot(v1, v2) / (np.vdot(v1, v1) * np.vdot(v2, v2)) ** 0.5
        fsc.append(float(p.real))
        prev_mask = mask

    return pd.DataFrame(dict(pixres=np.arange(D // 2) / D, fsc=fsc), dtype=float)


def main(args):
    vol1 = ImageSource.from_file(args.volumes[0])
    vol2 = ImageSource.from_file(args.volumes[1])
    fsc_vals = calculate_fsc(vol1.images(), vol2.images(), args.mask)

    if args.outtxt:
        logger.info(f"Saving FSC values to {args.outtxt}")
        fsc_vals.to_csv(args.outtxt, sep=" ", index=False)
    else:
        fsc_str = fsc_vals.round(4).to_csv(sep="\t", index=False)
        logger.info(f"\n{fsc_str}")

    if fsc_vals.shape[0] > 1 and (fsc_vals.fsc >= 0.5).any():
        res = fsc_vals.pixres[fsc_vals.fsc >= 0.5].max()
        logger.info("0.5: {:.7g} ang".format((1 / res) * args.Apix))
    else:
        logger.warning("0.5: N/A")

    if fsc_vals.shape[0] > 1 and (fsc_vals.fsc >= 0.143).any():
        res = fsc_vals.pixres[fsc_vals.fsc >= 0.143].max()
        logger.info("0.143: {:.7g} ang".format((1 / res) * args.Apix))
    else:
        logger.warning("0.143: N/A")

    if args.plot:
        if isinstance(args.plot, bool):
            if args.outtxt:
                plot_file = "".join([os.path.splitext(args.outtxt)[0], ".png"])
            else:
                plot_file = "fsc-plot.png"
        else:
            plot_file = str(args.plot)

        logger.info(f"Saving plot to {plot_file}")
        create_fsc_plot(fsc_vals=fsc_vals, outfile=plot_file, Apix=args.Apix)
