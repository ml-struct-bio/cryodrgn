"""
Downsample an image stack or volume by clipping fourier frequencies
"""

import argparse
import math
import os
import logging
import numpy as np
from cryodrgn import fft, utils
from cryodrgn.mrc import MRCHeader, MRCFile
from cryodrgn.source import ImageSource

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument(
        "mrcs", help="Input particles or volume (.mrc, .mrcs, .star, .cs, or .txt)"
    )
    parser.add_argument(
        "-D", type=int, required=True, help="New box size in pixels, must be even"
    )
    parser.add_argument(
        "-o",
        metavar="MRCS",
        type=os.path.abspath,
        required=True,
        help="Output projection stack (.mrcs)",
    )
    parser.add_argument(
        "-b",
        type=int,
        default=5000,
        help="Batch size for processing images (default: %(default)s)",
    )
    parser.add_argument(
        "--is-vol", action="store_true", help="Flag if input .mrc is a volume"
    )
    parser.add_argument(
        "--chunk",
        type=int,
        help="Chunksize (in # of images) to split particle stack when saving",
    )
    parser.add_argument(
        "--datadir",
        help="Optionally provide path to input .mrcs if loading from a .star or .cs file",
    )
    parser.add_argument(
        "--max-threads",
        type=int,
        default=16,
        help="Maximum number of CPU cores for parallelization (default: %(default)s)",
    )
    parser.add_argument(
        "--ind",
        type=os.path.abspath,
        metavar="PKL",
        help="Filter particle stack by these indices",
    )
    return parser


def mkbasedir(out):
    if not os.path.exists(os.path.dirname(out)):
        os.makedirs(os.path.dirname(out))


def warnexists(out):
    if os.path.exists(out):
        logger.warning("Warning: {} already exists. Overwriting.".format(out))


def main(args):
    mkbasedir(args.o)
    warnexists(args.o)
    assert args.o.endswith(".mrcs") or args.o.endswith(
        "mrc"
    ), "Must specify output in .mrc(s) file format"

    lazy = not args.is_vol

    ind = None
    if args.ind is not None:
        assert not args.is_vol
        logger.info(f"Filtering image dataset with {args.ind}")
        ind = utils.load_pkl(args.ind).astype(int)

    old = ImageSource.from_file(args.mrcs, lazy=lazy, indices=ind, datadir=args.datadir)

    oldD = old.D
    assert (
        args.D <= oldD
    ), f"New box size {args.D} cannot be larger than the original box size {oldD}"
    assert args.D % 2 == 0, "New box size must be even"

    D = args.D
    start = int(oldD / 2 - D / 2)
    stop = int(oldD / 2 + D / 2)

    # Downsample volume
    if args.is_vol:
        old = old.images()
        oldft = fft.htn_center(old)
        logger.info(oldft.shape)
        newft = oldft[start:stop, start:stop, start:stop]
        logger.info(newft.shape)
        new = np.array(fft.ihtn_center(newft)).astype(np.float32)
        logger.info(f"Saving {args.o}")
        MRCFile.write(args.o, array=new, is_vol=True)

    # Downsample images
    else:

        def transform_fn(chunk, indices):
            oldft = fft.ht2_center(chunk)
            newft = oldft[:, start:stop, start:stop]
            new = fft.iht2_center(newft)
            return new

        if args.chunk is None:
            logger.info("Saving {}".format(args.o))
            header = MRCHeader.make_default_header(
                nz=old.n, ny=D, nx=D, data=None, is_vol=args.is_vol
            )
            MRCFile.write(
                filename=args.o,
                array=old,
                header=header,
                is_vol=args.is_vol,
                transform_fn=transform_fn,
                chunksize=args.b,
            )

        else:
            # Downsample images, saving chunks of N images
            nchunks = math.ceil(len(old) / args.chunk)
            out_mrcs = [
                ".{}".format(i).join(os.path.splitext(args.o)) for i in range(nchunks)
            ]
            chunk_names = [os.path.basename(x) for x in out_mrcs]
            for i in range(nchunks):
                logger.info("Processing chunk {}".format(i))
                chunk = old[i * args.chunk : (i + 1) * args.chunk]

                header = MRCHeader.make_default_header(
                    nz=len(chunk), ny=D, nx=D, data=None, is_vol=args.is_vol
                )
                logger.info(f"Saving {out_mrcs[i]}")
                MRCFile.write(
                    filename=out_mrcs[i],
                    array=chunk,
                    header=header,
                    is_vol=args.is_vol,
                    transform_fn=transform_fn,
                    chunksize=args.b,
                )

            # Write a text file with all chunks
            out_txt = "{}.txt".format(os.path.splitext(args.o)[0])
            logger.info(f"Saving {out_txt}")
            with open(out_txt, "w") as f:
                f.write("\n".join(chunk_names))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    main(add_args(parser).parse_args())
