"""
Preprocess a dataset for more streamlined cryoDRGN training
"""

import argparse
import numpy as np
import math
import os
from typing import List
import logging
import torch
from cryodrgn import fft, utils
from cryodrgn.utils import window_mask
from cryodrgn.source import ImageSource
from cryodrgn.mrc import MRCFile

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument(
        "mrcs", help="Input particles or volume (.mrc, .mrcs, .star, or .txt)"
    )
    parser.add_argument(
        "-o", metavar="MRCS", type=os.path.abspath, required=True, help="Output .mrcs"
    )
    parser.add_argument(
        "--datadir",
        help="Optionally provide path to input .mrcs if loading from a .star or .cs file",
    )

    group = parser.add_argument_group("Image preprocessing settings")
    group.add_argument(
        "--ind",
        type=os.path.abspath,
        metavar="PKL",
        help="Filter particle stack by these indices",
    )
    group.add_argument(
        "-D", type=int, help="New box size in pixels (if downsampling), must be even"
    )
    group.add_argument(
        "--uninvert-data",
        dest="invert_data",
        action="store_false",
        help="Do not invert data sign",
    )
    group.add_argument(
        "--window-r",
        default=0.85,
        type=float,
        help="Circular windowing mask inner radius (default: %(default)s)",
    )
    group.add_argument(
        "--no-window",
        dest="window",
        action="store_false",
        help="Turn off real space windowing of dataset",
    )

    group = parser.add_argument_group("Extra arguments for volume generation")
    group.add_argument(
        "-b",
        type=int,
        default=5000,
        help="Batch size for processing images (default: %(default)s)",
    )
    group.add_argument(
        "--chunk",
        type=int,
        default=100000,
        help="Chunksize (in # of images) to split particle stack when saving",
    )
    group.add_argument(
        "--no-lazy",
        dest="lazy",
        action="store_false",
        help="Load whole dataset (faster than loading in batches)",
    )
    group.add_argument(
        "--max-threads",
        type=int,
        default=16,
        help="Maximum number of CPU cores for parallelization (default: %(default)s)",
    )

    group = parser.add_argument_group("GPU acceleratation")
    group.add_argument(
        "--use-cupy", action="store_true", help="Use cupy to replace numpy"
    )
    return parser


def mkbasedir(out):
    if not os.path.exists(os.path.dirname(out)):
        os.makedirs(os.path.dirname(out))


def warnexists(out):
    if os.path.exists(out):
        logger.warning(f"Warning: {out} already exists. Overwriting.")


def main(args):
    mkbasedir(args.o)
    warnexists(args.o)
    assert args.o.endswith(".mrcs") or args.o.endswith(
        ".txt"
    ), "Must specify output in .mrcs file format"

    # load images
    lazy = args.lazy

    # filter images
    ind = None
    if args.ind is not None:
        logger.info(f"Filtering image dataset with {args.ind}")
        ind = utils.load_pkl(args.ind).astype(int)

    images = ImageSource.from_file(
        args.mrcs, lazy=lazy, datadir=args.datadir, indices=ind
    )
    original_D = images.L

    logger.info(f"Loading {len(images)} {original_D}x{original_D} images")
    window = args.window
    invert_data = args.invert_data
    downsample = args.D and args.D < original_D
    if downsample:
        assert (
            args.D <= original_D
        ), f"New box size {args.D} cannot be larger than the original box size {original_D}"
        assert args.D % 2 == 0, "New box size must be even"
        start = int(original_D / 2 - args.D / 2)
        stop = int(original_D / 2 + args.D / 2)
        D = args.D
        logger.info(f"Downsampling images to {D}x{D}")
    else:
        D = original_D

    def preprocess(imgs: torch.Tensor):
        if window:
            imgs *= window_mask(original_D, args.window_r, 0.99)

        ret = fft.ht2_center(imgs)
        if invert_data:
            ret *= -1
        if downsample:
            ret = ret[:, start:stop, start:stop]
        ret = fft.symmetrize_ht(ret)
        return ret

    def preprocess_in_batches(imgs, b):
        ret = np.empty((len(imgs), D + 1, D + 1), dtype=np.float32)
        Nbatches = math.ceil(len(imgs) / b)
        for ii in range(Nbatches):
            logger.info(f"Processing batch of {b} images ({ii+1} of {Nbatches})")
            ret[ii * b : (ii + 1) * b, :, :] = torch.Tensor(  # type: ignore
                preprocess(imgs[ii * b : (ii + 1) * b])
            )
        return ret

    nchunks = math.ceil(len(images) / args.chunk)
    out_mrcs = [f".{i}.ft".join(os.path.splitext(args.o)) for i in range(nchunks)]
    chunk_names = [os.path.basename(x) for x in out_mrcs]
    for i in range(nchunks):
        logger.info(f"Processing chunk {i+1} of {nchunks}")
        chunk = images[i * args.chunk : (i + 1) * args.chunk]
        new = preprocess_in_batches(chunk, args.b)
        logger.info(f"New shape: {new.shape}")
        logger.info(f"Saving {out_mrcs[i]}")
        MRCFile.write(out_mrcs[i], new, is_vol=False)

    out_txt = f"{os.path.splitext(args.o)[0]}.ft.txt"
    logger.info(f"Saving summary txt file {out_txt}")
    with open(out_txt, "w") as f:
        f.write("\n".join(chunk_names))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    main(add_args(parser).parse_args())
