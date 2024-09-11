"""Translate a particle stack by applying 2D shifts.

Example usage
-------------
$ cryodrgn_utils translate_mrcs projections.mrcs trans.pkl -o projections.trans.mrcs

"""
import argparse
import os
import logging
import matplotlib.pyplot as plt
import numpy as np
from cryodrgn import fft, utils
from cryodrgn.source import ImageSource, write_mrc

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("mrcs", help="Input particles (.mrcs, .cs, .star, or .txt)")
    parser.add_argument("trans", help="Pose or translations pickle (.pkl)")
    parser.add_argument(
        "--tscale",
        type=float,
        default=1.0,
        help="Scale translations by this amount (default: %(default)s)",
    )
    parser.add_argument(
        "--datadir",
        help="Optionally overwrite path to starfile .mrcs if loading from a starfile",
    )
    parser.add_argument(
        "-o", type=os.path.abspath, required=True, help="Output particle stack (.mrcs)"
    )
    parser.add_argument("--out-png")


def plot_projections(out_png: str, imgs: np.ndarray) -> None:
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(10, 10))
    axes = axes.ravel()
    for i in range(min(len(imgs), 9)):
        axes[i].imshow(imgs[i])

    plt.savefig(out_png)


def main(args: argparse.Namespace) -> None:
    # load particles
    particles = ImageSource.from_file(args.mrcs, datadir=args.datadir)
    logger.info(particles.shape)

    trans = utils.load_pkl(args.trans)
    if type(trans) is tuple:
        trans = trans[1]
    trans *= args.tscale
    if np.any(trans > 1):
        raise ValueError(
            "Old pose format detected "
            "— translations must now be in units of fraction of box!"
        )
    trans *= particles.D  # convert to pixels
    assert len(trans) == particles.n

    new_dim = -particles.D / 2, particles.D / 2
    xx, yy = np.meshgrid(np.arange(*new_dim), np.arange(*new_dim))
    TCOORD = np.stack([xx, yy], axis=2) / particles.D  # DxDx2

    imgs = np.empty(particles.shape, dtype=np.float32)
    for ii in range(particles.n):
        if ii % 1000 == 0:
            logger.info(f"Processing image {ii}")

        ff = fft.fft2_center(particles[ii])
        tfilt = np.dot(TCOORD, trans[ii]) * -2 * np.pi
        tfilt = np.cos(tfilt) + np.sin(tfilt) * 1j
        imgs[ii] = fft.ifftn_center(ff * tfilt)

    logger.info(f"Writing {args.o}")
    write_mrc(args.o, imgs)

    if args.out_png:
        plot_projections(args.out_png, imgs[:9])
