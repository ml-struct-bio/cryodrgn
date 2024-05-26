"""Parse image poses from RELION .star file into a separate file for use in cryoDRGN.

This command is often used as a part of preparing inputs for training commands such as
`train_vae` and `abinit_homo` when particles are coming from a .star file.

Example usages
--------------
$ cryodrgn parse_pose_star particles_from_M.star -o pose.pkl

# override image parameters even if given in file
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

    group = parser.add_argument_group("Optionally provide missing image parameters")
    group.add_argument(
        "-D", type=int, help="override box size of reconstruction (pixels)"
    )
    group.add_argument(
        "--Apix",
        type=float,
        help="Pixel size (A); override if translations are specified in Angstroms",
    )
    return parser


def main(args):
    if not args.input.endswith(".star"):
        raise ValueError("Input file must be a .star file!")
    if not args.outpkl.endswith(".pkl"):
        raise ValueError("Output file must be a .pkl file (pickled Python format)!")

    resolution = None
    apix = None
    stardata = starfile.Starfile.load(args.input)
    N = len(stardata.df)
    logger.info(f"{N} particles")

    if stardata.relion31:  # Get image stats from data_optics table
        assert stardata.data_optics is not None
        optics_df = stardata.data_optics.df.set_index("_rlnOpticsGroup")
        apix = np.array(
            [
                float(optics_df.loc[g, "_rlnImagePixelSize"])
                for g in stardata.df["_rlnOpticsGroup"].values
            ]
        )
        resolution = np.array(
            [
                int(float(optics_df.loc[g, "_rlnImageSize"]))
                for g in stardata.df["_rlnOpticsGroup"].values
            ]
        )

    if resolution is None and "_rlnImageSize" in stardata.headers:
        resolution = np.array(
            [int(float(stardata.df["_rlnImageSize"][0])) for _ in range(N)]
        )

    if args.D is not None:
        resolution = np.array([args.D for _ in range(N)])
    if args.Apix is not None:
        apix = np.array([args.Apix for _ in range(N)])

    assert resolution is not None, "Must provide image size with -D"

    # parse rotations
    euler = np.zeros((N, 3))
    euler[:, 0] = stardata.df["_rlnAngleRot"]
    euler[:, 1] = stardata.df["_rlnAngleTilt"]
    euler[:, 2] = stardata.df["_rlnAnglePsi"]
    logger.info("Euler angles (Rot, Tilt, Psi):")
    logger.info(euler[0])
    logger.info("Converting to rotation matrix:")
    rot = np.asarray([utils.R_from_relion(*x) for x in euler])

    logger.info(rot[0])

    # parse translations
    trans = np.zeros((N, 2))
    if "_rlnOriginX" in stardata.headers and "_rlnOriginY" in stardata.headers:
        # translations in pixels
        trans[:, 0] = stardata.df["_rlnOriginX"]
        trans[:, 1] = stardata.df["_rlnOriginY"]

    elif (
        "_rlnOriginXAngst" in stardata.headers
        and "_rlnOriginYAngst" in stardata.headers
    ):
        # translation in Angstroms (Relion 3.1)
        assert apix is not None, (
            "Must provide --Apix argument to convert _rlnOriginXAngst "
            "and _rlnOriginYAngst translation units"
        )
        trans[:, 0] = stardata.df["_rlnOriginXAngst"]
        trans[:, 1] = stardata.df["_rlnOriginYAngst"]
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    main(add_args(parser).parse_args())
