"""Downsample an image stack or volume by clipping fourier frequencies.

Example usage
-------------
$ cryodrgn downsample my_particle_stack.mrcs -D 128 -o particles.128.mrcs
$ cryodrgn downsample my_particle_stack.mrcs -D 164 -o particles.164.mrcs \
                                                    --ind chosen_particles.pkl
$ cryodrgn downsample my_particle_stack.star -D 128 -o particles.128.mrcs \
                                             --datadir folder_with_subtilts/

# Try a smaller processing batch size if you are running into memory issues, or a
# larger size for faster processing
$ cryodrgn downsample my_particle_stack.txt -D 256 -o particles.256.mrcs -b 2000
$ cryodrgn downsample my_particle_stack.txt -D 256 -o particles.256.mrcs -b 20000

# This will create files particles.256.0.mrcs, particles.256.1.mrcs, ...,
# particles.256.i.mrcs, where i is equal to particle count // 10000, in addition to
# output file particles.256.txt that indexes all of them
$ cryodrgn downsample my_particle_stack.mrcs -D 256 -o particles.256.mrcs --chunk 10000

"""
import argparse
import math
import os
import logging
import numpy as np
import pandas as pd
from typing import Optional
from cryodrgn import fft, utils
from cryodrgn.source import ImageSource, MRCHeader, write_mrc

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "input", help="Input particles or volume (.mrc, .mrcs, .star, .cs, or .txt)"
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
        help="Size of chunks (in # of images, each in its own file) to split particle "
        "stack when saving",
    )
    parser.add_argument(
        "--datadir",
        help="Optionally provide folder containing input .mrcs files "
        "if loading from a .star or .cs file",
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


def mkbasedir(out):
    if not os.path.exists(os.path.dirname(out)):
        os.makedirs(os.path.dirname(out))


def warnexists(out):
    if os.path.exists(out):
        logger.warning("Warning: {} already exists. Overwriting.".format(out))


def downsample_mrc_images(
    src: ImageSource,
    new_D: int,
    out_fl: str,
    batch_size: int,
    chunk_size: Optional[int] = None,
):
    if new_D > src.D:
        raise ValueError(
            f"New box size {new_D} cannot be larger than the original box size {src.D}!"
        )

    old_apix = src.apix
    if old_apix is None:
        old_apix = 1.0

    new_apix = round(old_apix * src.D / new_D, 6)
    start = int(src.D / 2 - new_D / 2)
    stop = int(src.D / 2 + new_D / 2)

    if isinstance(new_apix, pd.Series):
        new_apix = new_apix.unique()

        if len(new_apix) > 1:
            logger.warning(
                f"Found multiple A/px values in {src.filenames}, using the first one "
                f"found {new_apix[0]} for the output .mrcs!"
            )

        new_apix = new_apix[0]

    def transform_fn(chunk, indices):
        oldft = fft.ht2_center(chunk)
        newft = oldft[:, start:stop, start:stop]
        new = fft.iht2_center(newft)
        return new

    if chunk_size is None:
        logger.info(f"Saving {out_fl}")

        header = MRCHeader.make_default_header(
            nz=src.n, ny=new_D, nx=new_D, Apix=new_apix, data=None, is_vol=False
        )

        write_mrc(
            filename=out_fl,
            array=src,
            header=header,
            is_vol=False,
            transform_fn=transform_fn,
            chunksize=batch_size,
        )

    # downsample images, saving one chunk of N images at a time
    else:
        nchunks = math.ceil(len(src) / chunk_size)
        out_mrcs = [
            ".{}".format(i).join(os.path.splitext(out_fl)) for i in range(nchunks)
        ]
        chunk_names = [os.path.basename(x) for x in out_mrcs]

        for i in range(nchunks):
            logger.info("Processing chunk {}".format(i))
            chunk = src[i * chunk_size : (i + 1) * chunk_size]

            header = MRCHeader.make_default_header(
                nz=len(chunk),
                ny=new_D,
                nx=new_D,
                Apix=new_apix,
                data=None,
                is_vol=False,
            )

            logger.info(f"Saving {out_mrcs[i]}")
            write_mrc(
                filename=out_mrcs[i],
                array=chunk,
                header=header,
                is_vol=False,
                transform_fn=transform_fn,
                chunksize=batch_size,
            )

        # Write a text file with all chunks
        out_txt = "{}.txt".format(os.path.splitext(out_fl)[0])
        logger.info(f"Saving {out_txt}")
        with open(out_txt, "w") as f:
            f.write("\n".join(chunk_names))


def main(args):
    mkbasedir(args.o)
    warnexists(args.o)
    assert args.o.endswith(".mrcs") or args.o.endswith(
        "mrc"
    ), "Must specify output in .mrc(s) file format"

    ind = None
    if args.ind is not None:
        assert not args.is_vol
        logger.info(f"Filtering image dataset with {args.ind}")
        ind = utils.load_pkl(args.ind).astype(int)

    old_src = ImageSource.from_file(
        args.input, lazy=not args.is_vol, indices=ind, datadir=args.datadir
    )
    assert (
        args.D <= old_src.D
    ), f"New box size {args.D} cannot be larger than the original box size {old_src.D}"
    assert args.D % 2 == 0, "New box size must be even"

    # downsample volume
    if args.is_vol:
        start = int(old_src.D / 2 - args.D / 2)
        stop = int(old_src.D / 2 + args.D / 2)

        old = old_src.images()
        oldft = fft.htn_center(old)
        logger.info(oldft.shape)
        newft = oldft[start:stop, start:stop, start:stop]
        logger.info(newft.shape)

        new = np.array(fft.ihtn_center(newft)).astype(np.float32)
        logger.info(f"Saving {args.o}")
        write_mrc(args.o, array=new, is_vol=True)

    # downsample images
    else:
        downsample_mrc_images(old_src, args.D, args.o, args.b, args.chunk)
