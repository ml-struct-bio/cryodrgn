"""Compute Fourier shell correlation between two volumes, applying an optional mask.

Example usage
-------------
$ cryodrgn_utils fsc volume1.mrc volume2.mrc

# Save FSC values to file and produce an FSC plot
$ cryodrgn_utils fsc vol1.mrc vol2.mrc -o fsc.txt -p

# Also apply a mask before computing FSCs
$ cryodrgn_utils fsc vol1.mrc vol2.mrc --mask test-mask.mrc -o fsc.txt -p fsc-plot.png

# Also apply phase randomization at Fourier shells for resolutions < 10 angstroms
$ cryodrgn_utils fsc vol1.mrc vol2.mrc --mask test-mask.mrc -o fsc.txt -p fsc-plot.png
                                       --corrected 10

"""
import os
import argparse
import logging
import numpy as np
import pandas as pd
import random
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
        "--corrected",
        type=float,
        help="use cryoSPARC-style high resolution phase randomization beyond this "
        "resolution to correct for possible effects of tight masking",
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
    vol1: np.ndarray,
    vol2: np.ndarray,
    mask: Optional[np.ndarray] = None,
    phase_randomization: Optional[float] = None,
) -> pd.DataFrame:
    if mask is not None:
        vol1 *= mask
        vol2 *= mask

    D = vol1.shape[0]
    x = np.arange(-D // 2, D // 2)
    x2, x1, x0 = np.meshgrid(x, x, x, indexing="ij")
    coords = np.stack((x0, x1, x2), -1)
    dists = (coords**2).sum(-1) ** 0.5
    assert dists[D // 2, D // 2, D // 2] == 0.0

    vol1 = fft.fftn_center(vol1)
    vol2 = fft.fftn_center(vol2)

    # logger.info(r[D//2, D//2, D//2:])
    prev_mask = np.zeros((D, D, D), dtype=bool)
    fsc = [1.0]
    for i in range(1, D // 2):
        mask = dists < i
        shell = np.where(mask & np.logical_not(prev_mask))
        v1 = vol1[shell]
        v2 = vol2[shell]
        p = np.vdot(v1, v2) / (np.vdot(v1, v1) * np.vdot(v2, v2)) ** 0.5
        fsc.append(float(p.real))
        prev_mask = mask

    if phase_randomization is not None:
        phase_D = int(phase_randomization * D)
        phase_mask = dists > phase_D
        phase_inds = list(range(phase_mask.sum()))
        random.shuffle(phase_inds)
        vol1[phase_mask] = vol1[phase_mask][phase_inds]
        random.shuffle(phase_inds)
        vol2[phase_mask] = vol2[phase_mask][phase_inds]

        for i in range(1, D // 2):
            mask = dists < i
            shell = np.where(mask & np.logical_not(prev_mask))

            if i > phase_D:
                v1 = vol1[shell]
                v2 = vol2[shell]
                p = float(
                    (np.vdot(v1, v2) / (np.vdot(v1, v1) * np.vdot(v2, v2)) ** 0.5).real
                )
                if not np.isnan(p) and p < 1.0:
                    fsc[i - 1] = (fsc[i - 1] - p) / (1 - p)

            prev_mask = mask

    return pd.DataFrame(dict(pixres=np.arange(D // 2) / D, fsc=fsc), dtype=float)


def get_fsc_thresholds(
    fsc_vals: pd.DataFrame, apix: float, verbose: bool = True
) -> tuple[float, float]:
    if ((fsc_vals.pixres > 0) & (fsc_vals.fsc >= 0.5)).any():
        res_05 = fsc_vals.pixres[fsc_vals.fsc >= 0.5].max()
        if verbose:
            logger.info("res @ FSC=0.5: {:.4g} ang".format((1 / res_05) * apix))
    else:
        res_05 = None
        if verbose:
            logger.warning("res @ FSC=0.5: N/A")

    if ((fsc_vals.pixres > 0) & (fsc_vals.fsc >= 0.143)).any():
        res_143 = fsc_vals.pixres[fsc_vals.fsc >= 0.143].max()
        if verbose:
            logger.info("res @ FSC=0.143: {:.4g} ang".format((1 / res_143) * apix))
    else:
        res_143 = None
        if verbose:
            logger.warning("res @ FSC=0.143: N/A")

    return res_05, res_143


def main(args: argparse.Namespace) -> None:
    vol1 = ImageSource.from_file(args.volumes[0])
    vol2 = ImageSource.from_file(args.volumes[1])
    mask = ImageSource.from_file(args.mask).images() if args.mask is not None else None

    if args.corrected is not None:
        if args.corrected >= 1:
            args.corrected = (args.corrected / args.Apix) ** -1

    fsc_vals = calculate_fsc(vol1.images(), vol2.images(), mask, args.corrected)

    if args.outtxt:
        logger.info(f"Saving FSC values to {args.outtxt}")
        fsc_vals.to_csv(args.outtxt, sep=" ", index=False)
    else:
        fsc_str = fsc_vals.round(4).to_csv(sep="\t", index=False)
        logger.info(f"\n{fsc_str}")

    _ = get_fsc_thresholds(fsc_vals, args.Apix)

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
