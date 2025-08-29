"""Downsample an image stack or volume to a lower resolution by clipping Fourier freqs.

Example usage
-------------
# Downsample an image stack file to a new stack file
$ cryodrgn downsample my_particle_stack.mrcs -D 128 -o particles.128.mrcs

# Downsample an image stack and also apply a filtering index
$ cryodrgn downsample my_particle_stack.mrcs -D 164 -o particles.164.mrcs \
                                                    --ind chosen_particles.pkl

# Downsample an image stack saved across many files referenced by a .star file to a
# single new image stack file
$ cryodrgn downsample my_particle_stack.star -D 128 -o particles.128.mrcs \
                                             --datadir folder_with_subtilts/

# Downsample a .star image stack, preserving the original image stack file structure
# and creating a new .star file with image stack saved alongside it (i.e. no
# --datadir necessary for future use)
$ cryodrgn downsample my_particle_stack.star -D 128 -o particles.128.star \
                                             --datadir folder_with_subtilts/

# Same case as above except specifying a new --datadir using --outdir:
$ cryodrgn downsample my_particle_stack.star -D 128 -o particles.128.star \
                                             --datadir folder_with_subtilts/ \
                                             --outdir my_new_datadir/

# Downsample an image stack saved across many files referenced by a .txt file to a
# single new image stack file
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
import time
import numpy as np
import torch
from collections.abc import Iterable
from typing import Optional
from cryodrgn import fft, utils
from cryodrgn.source import ImageSource, StarfileSource, TxtFileSource
from cryodrgn.mrcfile import write_mrc, MRCHeader

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    """The command-line arguments for use with `cryodrgn downsample`."""

    parser.add_argument(
        "input", help="Input particles or volume (.mrc, .mrcs, .star, .cs, or .txt)"
    )
    parser.add_argument(
        "-D", type=int, required=True, help="New box size in pixels, must be even"
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=os.path.abspath,
        required=True,
        help="Output projection stack (.mrc, .mrcs, .star, or .txt)",
    )
    parser.add_argument(
        "--outdir",
        type=os.path.abspath,
        help="Output image stack directory, (default: placed alongside --outfile)",
    )
    parser.add_argument(
        "-b",
        "--batch-size",
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


def downsample_mrc_images(
    src: ImageSource,
    new_box_size: int,
    out_fl: str,
    batch_size: int,
    chunk_size: Optional[int] = None,
) -> None:
    """Downsample the images in a single particle stack into a new .mrcs file.

    This utility function also simplifies handling each of the individual .mrcs files
    listed in a .txt or .star file into a new collection of downsampled .mrcs files.

    Arguments
    ---------
    src         A loaded particle image stack.
    new_box_size       The new downsampled box size in pixels.
    out_fl      The output .mrcs file name.
    batch_size  The batch size for processing images;
                useful for avoiding out-of-memory issues.
    chunk_size  If given, divide output into files each containing this many images.

    """
    if new_box_size > src.D:
        raise ValueError(
            f"New box size {new_box_size} cannot be larger "
            f"than the original box size {src.D}!"
        )

    old_apix = src.apix
    if old_apix is None:
        old_apix = 1.0

    new_apix = np.round(old_apix * src.D / new_box_size, 6)
    if isinstance(new_apix, Iterable):
        new_apix = tuple(set(new_apix))

        if len(new_apix) > 1:
            logger.warning(
                f"Found multiple A/px values in {src.filenames}, using the first one "
                f"found {new_apix[0]} for the output .mrcs!"
            )
        new_apix = new_apix[0]

    def downsample_transform(chunk: torch.Tensor, indices: np.ndarray) -> torch.Tensor:
        """Downsample an image array by clipping Fourier frequencies."""

        start = int(src.D / 2 - new_box_size / 2)
        stop = int(src.D / 2 + new_box_size / 2)
        oldft = fft.ht2_center(chunk)
        newft = oldft[:, start:stop, start:stop]
        new_chunk = fft.iht2_center(newft)

        return new_chunk

    if chunk_size is None:
        logger.info(f"Saving {out_fl}")

        header = MRCHeader.make_default_header(
            nz=src.n,
            ny=new_box_size,
            nx=new_box_size,
            Apix=new_apix,
            dtype=src.dtype,
            is_vol=False,
        )
        src.write_mrc(
            output_file=out_fl,
            header=header,
            transform_fn=downsample_transform,
            chunksize=batch_size,
        )

    # Downsample images, saving one chunk of N images at a time
    else:
        nchunks = math.ceil(len(src) / chunk_size)
        out_mrcs = [f".{i}".join(os.path.splitext(out_fl)) for i in range(nchunks)]
        chunk_names = [os.path.basename(x) for x in out_mrcs]

        for i, out_mrc in enumerate(out_mrcs):
            logger.info(f"Processing chunk {i}")
            chunk_indices = np.arange(i * chunk_size, min((i + 1) * chunk_size, src.n))

            header = MRCHeader.make_default_header(
                nz=len(chunk_indices),
                ny=new_box_size,
                nx=new_box_size,
                Apix=new_apix,
                dtype=src.dtype,
                is_vol=False,
            )

            logger.info(f"Saving {out_mrcs[i]}")
            src.write_mrc(
                output_file=out_mrc,
                header=header,
                indices=chunk_indices,
                transform_fn=downsample_transform,
                chunksize=batch_size,
            )

        # Write a text file with all chunks
        out_txt = f"{os.path.splitext(out_fl)[0]}.txt"
        logger.info(f"Saving {out_txt}")
        with open(out_txt, "w") as f:
            f.write("\n".join(chunk_names))


def main(args: argparse.Namespace) -> None:
    """Downsampling the given particle stack into a new stack (see `add_args` above)."""

    utils.create_basedir(args.outfile)
    utils.warn_file_exists(args.outfile)
    out_ext = os.path.splitext(args.outfile)[1]
    t0 = time.time()

    # Load filtering indices if provided
    ind = None
    if args.ind is not None:
        assert not args.is_vol
        logger.info(f"Filtering image dataset with {args.ind}")
        ind = utils.load_pkl(args.ind).astype(int)

    # Load image data (using pointers instead of raw data unless loading a volume)
    src = ImageSource.from_file(
        args.input,
        lazy=not args.is_vol,
        indices=ind,
        datadir=args.datadir,
        max_threads=args.max_threads,
    )
    if args.D > src.D:
        raise ValueError(
            f"New box size {args.D=} can't be larger "
            f"than the original box size {src.D}!"
        )
    if args.D % 2 != 0:
        raise ValueError(f"New box size {args.D=} is not even!")

    # Downsample a single file containing a volume into another volume
    if args.is_vol:
        start = int(src.D / 2 - args.D / 2)
        stop = int(src.D / 2 + args.D / 2)

        old_imgs = src.images()
        oldft = fft.htn_center(old_imgs)
        logger.info(oldft.shape)
        newft = oldft[start:stop, start:stop, start:stop]
        logger.info(newft.shape)

        new = np.array(fft.ihtn_center(newft)).astype(np.float32)
        logger.info(f"Saving {args.outfile}")
        write_mrc(args.outfile, array=new, is_vol=True)

    # Downsample images into a .mrcs image stack file, no matter the input format was
    elif out_ext in {".mrcs", ".mrc"}:
        downsample_mrc_images(src, args.D, args.outfile, args.batch_size, args.chunk)

    # Downsample images referenced in a .star file, using the original image stack file
    # structure where possible
    elif out_ext in {".star", ".txt"}:
        if args.chunk is not None:
            raise ValueError("Chunked output only available for .mrcs output!")

        if out_ext == ".star" and not isinstance(src, StarfileSource):
            raise ValueError(
                f"To write output to a .star file you must use a .star file as input, "
                f"instead was given {args.input=} !"
            )
        if out_ext == ".txt" and not isinstance(src, TxtFileSource):
            raise ValueError(
                f"To specify a .txt file as output you must use a .txt file as input, "
                f"instead was given {args.input=} !"
            )

        outdir = args.outdir or os.path.dirname(args.outfile)
        if out_ext == ".txt" and not os.path.isabs(outdir):
            outdir = os.path.join(
                os.path.dirname(os.path.abspath(args.outfile)),
                outdir,
            )
        outdir = os.path.abspath(outdir)
        os.makedirs(outdir, exist_ok=True)
        logger.info(f"Storing downsampled stacks in new --datadir `{outdir}`...")

        basepath = os.path.dirname(os.path.commonprefix([p for p, _ in src.sources]))
        newpaths = {
            oldpath: os.path.join(outdir, os.path.relpath(oldpath, basepath))
            for oldpath, _ in src.sources
        }
        for fl, fl_src in src.sources:
            os.makedirs(os.path.dirname(newpaths[fl]), exist_ok=True)
            downsample_mrc_images(
                fl_src, args.D, newpaths[fl], args.batch_size, chunk_size=None
            )

        src.df["__mrc_filepath"] = src.df["__mrc_filepath"].map(newpaths)
        if out_ext == ".star":
            src.df["__mrc_filename"] = [
                os.path.relpath(newpath, outdir) for newpath in src.df["__mrc_filepath"]
            ]
        else:
            src.df["__mrc_filename"] = [
                os.path.relpath(newpath, os.path.dirname(args.outfile))
                for newpath in src.df["__mrc_filepath"]
            ]

        if "_rlnImageName" in src.df.columns:
            src.df["_rlnImageName"] = [
                "@".join([old_name.split("@")[0], os.path.relpath(newpath, outdir)])
                for old_name, newpath in zip(
                    src.df["_rlnImageName"].values, src.df["__mrc_filepath"].values
                )
            ]
        if isinstance(src, StarfileSource) and src.resolution is not None:
            apix = src.apix or 1.0
            src.set_optics_values(
                "_rlnImagePixelSize", np.round(apix * src.resolution / args.D, 6)
            )

        src.write(args.outfile)

    else:
        raise ValueError(
            f"Unrecognized output extension `{out_ext}` "
            f"not in {{.mrc,.mrcs,.star,.txt}}!"
        )

    logger.info(f"Downsampling completed in {time.time() - t0:.1f} seconds")
