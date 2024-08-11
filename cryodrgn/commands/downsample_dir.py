"""Downsample image stack distributed across files referenced in a .star, .cs, or .txt.

Example usage
-------------
# will save downsampled .mrcs files in my_stack_256/ to downsampled.128/
$ cryodrgn downsample_dir my_particle_stack.star --datadir my_stack_256/ \
                                                 -D 128 -o particles.128.star

"""
import argparse
import os
import logging
from cryodrgn import utils
from cryodrgn.source import ImageSource, MRCDataFrameSource, StarfileSource
from cryodrgn.commands.downsample import mkbasedir, warnexists, downsample_mrc_images

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "input",
        help="Input particles (.star, .cs, or .txt)",
    )
    parser.add_argument(
        "--datadir",
        help="Optionally provide folder containing input .mrcs files "
        "if loading from a .star or .cs file",
    )

    parser.add_argument(
        "-D", type=int, required=True, help="New box size in pixels, must be even"
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=os.path.abspath,
        required=True,
        help="Output projection stack (.star, .cs, or .txt)",
    )
    parser.add_argument(
        "--outdir",
        type=os.path.abspath,
        help="Output image stack directory, (default: `downsampled.<-D>/`)",
    )
    parser.add_argument(
        "-b",
        type=int,
        default=5000,
        help="Batch size for processing images (default: %(default)s)",
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


def main(args):
    assert args.D % 2 == 0, "New box size must be even"
    mkbasedir(args.outfile)
    warnexists(args.outfile)

    ind = None
    if args.ind is not None:
        logger.info(f"Filtering image dataset with {args.ind}")
        ind = utils.load_pkl(args.ind).astype(int)

    src = ImageSource.from_file(
        filepath=args.input, lazy=True, indices=ind, datadir=args.datadir
    )
    if not isinstance(src, MRCDataFrameSource):
        raise TypeError(
            f"Input file {args.input} is not an index of a collection of .mrc/.mrcs "
            f"files such as .star, .cs, or .txt, "
            f"see `cryodrgn downsample for other input formats!"
        )

    outdir = args.outdir or f"downsampled.{args.D}"
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)
    logger.info(f"Storing downsampled stacks in new --datadir `{outdir}`...")

    new_fls = dict()
    for fl, fl_src in src.sources:
        new_fls[fl] = os.path.join(outdir, os.path.basename(fl))
        downsample_mrc_images(fl_src, args.D, new_fls[fl], args.b, chunk_size=None)

    src.df["__mrc_filepath"] = src.df["__mrc_filepath"].map(new_fls)

    if isinstance(src, StarfileSource):
        if src.relion31:
            if "_rlnImagePixelSize" in src.data_optics:
                src.data_optics["_rlnImagePixelSize"] = round(
                    src.data_optics["_rlnImagePixelSize"] * src.resolution / args.D, 6
                )

    src.write(args.outfile)
