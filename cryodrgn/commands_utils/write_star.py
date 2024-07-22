"""Create a Relion .star file from a given particle stack and CTF parameters.

Example usage
-------------
$ cryodrgn_utils write_star particles.128.mrcs -o particles.128.star --ctf ctf.pkl
$ cryodrgn_utils write_star particles.128.mrcs -o particles.128.star --ctf ctf.pkl \
                                               --ind good-ind.pkl

"""
import argparse
import os
import numpy as np
import pandas as pd
import logging
from cryodrgn import utils
from cryodrgn.source import ImageSource, StarfileSource
from cryodrgn.starfile import Starfile
from cryodrgn.mrc import MRCHeader

logger = logging.getLogger(__name__)


CTF_HEADERS = [
    "_rlnImageSize",
    "_rlnImagePixelSize",
    "_rlnDefocusU",
    "_rlnDefocusV",
    "_rlnDefocusAngle",
    "_rlnVoltage",
    "_rlnSphericalAberration",
    "_rlnAmplitudeContrast",
    "_rlnPhaseShift",
]

POSE_HDRS = [
    "_rlnAngleRot",
    "_rlnAngleTilt",
    "_rlnAnglePsi",
    "_rlnOriginX",
    "_rlnOriginY",
]


def add_args(parser):
    parser.add_argument("particles", help="Input particles (.mrcs, .txt, .star)")

    parser.add_argument(
        "-o", type=os.path.abspath, required=True, help="Output .star file"
    )

    parser.add_argument("--ctf", help="Input ctf.pkl")
    parser.add_argument("--poses", help="Optionally include pose.pkl")
    parser.add_argument(
        "--ind", help="Optionally filter by array of selected indices (.pkl)"
    )

    # TODO: is this needed any more?
    parser.add_argument(
        "--full-path",
        action="store_true",
        help="Write the full path to particles (default: relative paths)",
    )

    parser.add_argument(
        "--datadir",
        type=str,
        # TODO Check if .cs files work with this (or remove .cs here!)
        help="Path prefix to particle stack if loading relative paths from a .star or .cs file",
        default="",
    )

    parser.add_argument(
        "--relion30",
        action="store_true",
        help="Write output in RELION 3.0 format instead of the default 3.1 format.",
    )


def main(args):
    assert args.o.endswith(".star"), "Output file must be .star file"

    input_ext = os.path.splitext(args.particles)[-1]
    assert input_ext in (
        ".mrcs",
        ".txt",
        ".star",
    ), "Input file must be .mrcs/.txt/.star"

    # Either accept an input star file, or an input .mrcs/.txt with optional ctf/pose pkl file(s)
    ctf = poses = eulers = trans = None
    if input_ext == ".star":
        assert (
            args.poses is None
        ), "--poses cannot be specified when input is a starfile (poses are obtained from starfile)"
        assert (
            args.ctf is None
        ), "--ctf cannot be specified when input is a starfile (ctf information are obtained from starfile)"
    else:
        if not args.ctf:
            raise ValueError("--ctf must be specified when input is not a starfile!")

    # TODO REMOVE:ImageSource->StarFileSource->_MRCDataFrameSource
    # TODO REMOVE:datadir is passed through and used to set __mrc_filename values with realtive directories
    # TODO REMOVE:indices is used via super. If lazy=True, not utilised at all. If lazy=False, it is used to create an subselected stack of mrcs (in np fom)
    particles = ImageSource.from_file(args.particles, datadir=args.datadir, lazy=True) # indices=args.ind)
    if args.ctf:
        ctf = utils.load_pkl(args.ctf)
        assert ctf.shape[1] == 9, "Incorrect CTF pkl format"
        assert len(particles) == len(
            ctf
        ), f"{len(particles)} != {len(ctf)}, Number of particles != number of CTF parameters"

    if args.poses:
        poses = utils.load_pkl(args.poses)
        assert len(particles) == len(
            poses[0]
        ), f"{len(particles)} != {len(poses)}, Number of particles != number of poses"
    logger.info(f"{len(particles)} particles in {args.particles}")

    ind = np.arange(particles.n)
    if args.ind:
        ind = utils.load_pkl(args.ind)
        logger.info(f"Filtering to {len(ind)} particles")
        if ctf is not None:
            ctf = ctf[ind]
        if poses is not None:
            poses = (poses[0][ind], poses[1][ind])

    if input_ext == ".star":
        assert isinstance(particles, StarfileSource)
        df = particles.df.loc[ind]
        if ~args.relion30:
            optics = Starfile.load(args.particles).data_optics
        else:
            optics = None
    else:
        if input_ext == ".txt":
            with open(args.particles, "r") as f:
                mrcs_files = f.read().splitlines()

            base = os.path.dirname(os.path.abspath(args.particles))
            counts = [
                MRCHeader.parse(os.path.join(base, f)).fields["nz"] for f in mrcs_files
            ]
            ind_lbls = np.concatenate([np.arange(count) for count in counts])[ind] + 1
            image_names = np.repeat(mrcs_files, counts)[ind]

        else:
            ind_lbls = [str(i + 1) for i in ind]
            image_names = particles.filenames[ind]

        if args.full_path:
            image_names = [os.path.abspath(image_name) for image_name in image_names]
        names = [f"{lbl}@{name}" for lbl, name in zip(ind_lbls, image_names)]

        # convert poses
        if poses is not None:
            eulers = utils.R_to_relion_scipy(poses[0])
            D = particles[0].shape[-1]
            trans = poses[1] * D  # convert from fraction to pixels

        # Create a new dataframe with required star file headers
        data = {"_rlnImageName": names}
        ctf_cols = {2, 3, 4, 8}
        if ctf is not None:
            for ctf_col in ctf_cols:
                data[CTF_HEADERS[ctf_col]] = ctf[:, ctf_col]

        # figure out what the optics groups are using voltage, spherical aberration,
        # and amplitude contrast to assign a unique group to each particle
        optics_cols = list(set(range(9)) - ctf_cols)
        optics_headers = [CTF_HEADERS[optics_col] for optics_col in optics_cols]
        optics_groups, optics_indx = np.unique(
            ctf[:, optics_cols], return_inverse=True, axis=0
        )
        data["_rlnOpticsGroup"] = optics_indx + 1
        optics_groups = pd.DataFrame(optics_groups, columns=optics_headers)
        optics_groups["_rlnOpticsGroup"] = np.array(
            [str(i + 1) for i in range(optics_groups.shape[0])],
        )
        optics = Starfile(df=optics_groups, relion31=False, headers=None)

        if eulers is not None and trans is not None:
            for i in range(3):
                data[POSE_HDRS[i]] = eulers[:, i]  # type: ignore
            for i in range(2):
                data[POSE_HDRS[3 + i]] = trans[:, i]

        df = pd.DataFrame(data=data)
    
    # TODO check if setting __mrc_index __mrc_filename or 
    # __mrc_filepath for starfile input to starfile output
    # is relion compatible
    s = Starfile(headers=None, df=df, relion31=not args.relion30, data_optics=optics)
    s.write(args.o)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
