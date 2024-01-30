"""Compute Fourier shell correlation between two volumes.

Example usages
--------------
$ cryodrgn_utils fsc volume1.mrc volume2.mrc -o fsc.txt -p
$ cryodrgn_utils fsc volume1.mrc volume2.mrc --ind 5k-ind.pkl -o fsc.txt

"""
import argparse
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from typing import Optional
from cryodrgn import fft
from cryodrgn.source import ImageSource

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument("volumes", nargs=2, help="volumes to compare")

    parser.add_argument("--mask")
    parser.add_argument(
        "--plot", "-p", action="store_true", help="create plot of FSC curve?"
    )
    parser.add_argument("--Apix", type=float, default=1)
    parser.add_argument("-o", help="Output")


def calculate_fsc(
    vol1: torch.Tensor, vol2: torch.Tensor, mask_file: Optional[str] = None
) -> pd.Series:

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

    return pd.Series(fsc, index=np.arange(D // 2) / D)


def main(args):
    vol1 = ImageSource.from_file(args.volumes[0])
    vol2 = ImageSource.from_file(args.volumes[1])
    fsc_vals = calculate_fsc(vol1.images(), vol2.images(), args.mask)

    if args.o:
        np.savetxt(args.o, fsc_vals.values)
    else:
        logger.info(fsc_vals)

    if (fsc_vals >= 0.5).any():
        res = fsc_vals[fsc_vals >= 0.5].index.max()
        logger.info("0.5: {}".format(1 / res * args.Apix))
    if (fsc_vals >= 0.143).any():
        res = fsc_vals[fsc_vals >= 0.143].index.max()
        logger.info("0.143: {}".format(1 / res * args.Apix))

    if args.plot:
        plt.plot(fsc_vals.index, fsc_vals.values)
        plt.ylim((0, 1))
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
