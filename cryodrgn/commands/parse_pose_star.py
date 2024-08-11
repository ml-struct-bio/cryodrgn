"""Parse image poses from RELION .star file into a separate file for use in cryoDRGN.

This command is often used as a part of preparing inputs for training commands such as
`train_vae` and `abinit_homo` when particles are coming from a .star file.

Example usage
-------------
$ cryodrgn parse_pose_star particles_from_M.star -o pose.pkl

# override image parameters even if given in file
$ cryodrgn parse_pose_star particles_from_M.star -o pose.pkl -D 294 --Apix 1.7

"""
import argparse
import os
import pickle
import logging
import numpy as np
from cryodrgn import utils
from cryodrgn.starfile import Starfile

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input", help="RELION .star file")
    parser.add_argument(
        "--outpkl",
        "-o",
        metavar="PKL",
        type=os.path.abspath,
        required=True,
        help="Output pose.pkl",
    )

    group = parser.add_argument_group("Optionally provide missing image parameters")
    group.add_argument(
        "-D", type=int, help="override box size of reconstruction (pixels)"
    )
    group.add_argument(
        "--Apix",
        type=float,
        help="Pixel size (A); override if translations are specified in Angstroms",
    )


def main(args: argparse.Namespace) -> None:
    if not args.input.endswith(".star"):
        raise ValueError("Input file must be a .star file!")
    if not args.outpkl.endswith(".pkl"):
        raise ValueError("Output file must be a .pkl file (pickled Python format)!")

    starfile = Starfile(args.input)
    logger.info(f"{len(starfile)} particles")
    apix, resolution = starfile.apix, starfile.resolution

    if args.D is not None:
        resolution = np.array([args.D for _ in range(len(starfile))])
    if args.Apix is not None:
        apix = np.array([args.Apix for _ in range(len(starfile))])
    if resolution is None:
        raise ValueError(
            f"Must provide image size with -D as none found in `{args.input}`!"
        )

    # parse rotations
    euler = np.zeros((len(starfile), 3))
    euler[:, 0] = starfile.df["_rlnAngleRot"]
    euler[:, 1] = starfile.df["_rlnAngleTilt"]
    euler[:, 2] = starfile.df["_rlnAnglePsi"]
    logger.info("Euler angles (Rot, Tilt, Psi):")
    logger.info(euler[0])
    logger.info("Converting to rotation matrix:")
    rot = np.asarray([utils.R_from_relion(*x) for x in euler])

    logger.info(rot[0])

    # parse translations
    trans = np.zeros((len(starfile), 2))
    if "_rlnOriginX" in starfile.df.columns and "_rlnOriginY" in starfile.df.columns:
        # translations in pixels
        trans[:, 0] = starfile.df["_rlnOriginX"]
        trans[:, 1] = starfile.df["_rlnOriginY"]
    elif (
        "_rlnOriginXAngst" in starfile.df.columns
        and "_rlnOriginYAngst" in starfile.df.columns
    ):
        # translation in Angstroms (Relion 3.1)
        if apix is None:
            raise ValueError(
                f"Must provide --Apix argument to convert _rlnOriginXAngst and "
                f"_rlnOriginYAngst translation units as A/px not "
                f"found in `{args.input}`!"
            )
        trans[:, 0] = starfile.df["_rlnOriginXAngst"]
        trans[:, 1] = starfile.df["_rlnOriginYAngst"]
        trans /= apix.reshape(-1, 1)
    else:
        logger.warning(
            "Warning: Neither _rlnOriginX/Y nor _rlnOriginX/YAngst found. "
            "Defaulting to 0s."
        )

    logger.info("Translations (pixels):")
    logger.info(trans[0])

    # convert translations from pixels to fraction
    trans /= resolution.reshape(-1, 1)

    # write output
    logger.info(f"Writing {args.outpkl}")
    with open(args.outpkl, "wb") as f:
        pickle.dump((rot, trans), f)
