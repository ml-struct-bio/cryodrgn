"""Parse image poses from RELION .star file into a separate file for use in cryoDRGN.

This command is often used as a part of preparing inputs for training commands such as
`train_vae` and `abinit_homo` when particles are coming from a .star file.

Example usages
--------------
$ cryodrgn parse_pose_star particles_from_M.star -o pose.pkl -D 294 --Apix 1.7

"""
import argparse
import os
import pickle
import logging
import numpy as np
from cryodrgn import starfile, utils

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument("input", help="RELION .star file")
    parser.add_argument(
        "--outpkl",
        "-o",
        metavar="PKL",
        type=os.path.abspath,
        required=True,
        help="Output pose.pkl",
    )

    parser.add_argument(
        "--ignore-optics",
        action="store_true",
        help="ignore presence of optics tables from RELION 3.1 format",
    )

    group = parser.add_argument_group("Optionally provide missing image parameters")
    group.add_argument("-D", type=int, help="Box size of reconstruction (pixels)")
    group.add_argument(
        "--Apix",
        type=float,
        help="Pixel size (A); Required if translations are specified in Angstroms",
    )
    return parser


def main(args):
    if not args.input.endswith(".star"):
        raise ValueError("Input file must be a .star file!")
    if not args.outpkl.endswith(".pkl"):
        raise ValueError("Output file must be a .pkl file (pickled Python format)!")

    s = starfile.Starfile.load(args.input)
    if s.relion31:  # Get image stats from data_optics table
        assert s.data_optics is not None

        if len(s.data_optics.df) > 1:
            if not args.ignore_optics:
                raise ValueError(
                    "Datasets with only one optics group are supported "
                    "— use --ignore-optics flag if necessary!"
                )

        args.Apix = float(s.data_optics.df["_rlnImagePixelSize"][0])
        args.D = int(float(s.data_optics.df["_rlnImageSize"][0]))

    if args.D is None and "_rlnImageSize" in s.headers:
        args.D = int(float(s.df["_rlnImageSize"][0]))
    assert args.D is not None, "Must provide image size with -D"

    N = len(s.df)
    logger.info(f"{N} particles")

    # parse rotations
    euler = np.zeros((N, 3))
    euler[:, 0] = s.df["_rlnAngleRot"]
    euler[:, 1] = s.df["_rlnAngleTilt"]
    euler[:, 2] = s.df["_rlnAnglePsi"]
    logger.info("Euler angles (Rot, Tilt, Psi):")
    logger.info(euler[0])
    logger.info("Converting to rotation matrix:")
    rot = np.asarray([utils.R_from_relion(*x) for x in euler])

    logger.info(rot[0])

    # parse translations
    trans = np.zeros((N, 2))
    if "_rlnOriginX" in s.headers and "_rlnOriginY" in s.headers:
        # translations in pixels
        trans[:, 0] = s.df["_rlnOriginX"]
        trans[:, 1] = s.df["_rlnOriginY"]

    elif "_rlnOriginXAngst" in s.headers and "_rlnOriginYAngst" in s.headers:
        # translation in Angstroms (Relion 3.1)
        assert (
            args.Apix is not None
        ), "Must provide --Apix argument to convert _rlnOriginXAngst and _rlnOriginYAngst translation units"
        trans[:, 0] = s.df["_rlnOriginXAngst"]
        trans[:, 1] = s.df["_rlnOriginYAngst"]
        trans /= args.Apix

    else:
        logger.warning(
            "Warning: Neither _rlnOriginX/Y nor _rlnOriginX/YAngst found. Defaulting to 0s."
        )

    logger.info("Translations (pixels):")
    logger.info(trans[0])

    # convert translations from pixels to fraction
    trans /= args.D

    # write output
    logger.info(f"Writing {args.outpkl}")
    with open(args.outpkl, "wb") as f:
        pickle.dump((rot, trans), f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    main(add_args(parser).parse_args())
