"""Create a Relion .star file from a given particle stack and CTF parameters.

Example usage
-------------
# If using a .mrcs or .txt stack, you must provide an additional CTF file
$ cryodrgn_utils write_star particles.128.mrcs -o particles.128.star --ctf ctf.pkl

# You can also add particle pose data, and filter final output by given particle index
$ cryodrgn_utils write_star particles.128.mrcs -o particles.128.star --ctf ctf.pkl \
                            --poses pose.pkl --ind good-ind.pkl

# Can be used to filter a particle stack that is already a .star file, no CTF needed
$ cryodrgn_utils write_star particles.128.star -o particles.128_good.star \
                            --ind good-ind.pkl

# If using a .txt stack, you can use `--full-path` to avoid needing to specify
# any `--datadir` when using the new .star stack
$ cryodrgn_utils write_star particles.256.txt -o particles.256.star --ctf ctf.pkl \
                            --full-path

"""
import argparse
import os
import numpy as np
import pandas as pd
import logging
from cryodrgn import utils
from cryodrgn.source import ImageSource, StarfileSource, MRCHeader
from cryodrgn.starfile import write_star

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


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("particles", help="Input particles (.mrcs, .txt, .star)")
    parser.add_argument(
        "-o", type=os.path.abspath, required=True, help="Output .star file"
    )

    parser.add_argument("--ctf", help="Input ctf.pkl")
    parser.add_argument("--poses", help="Optionally include pose.pkl")
    parser.add_argument(
        "--ind", help="Optionally filter by array of selected indices (.pkl)"
    )
    parser.add_argument(
        "--datadir",
        type=str,
        help="Path prefix to particle stack if input contains relative paths "
        "or absolute paths that will be replaced",
    )
    parser.add_argument(
        "--full-path",
        action="store_true",
        help="Update relative particle paths to absolute paths (default: leave as is)",
    )
    parser.add_argument(
        "--relion30",
        action="store_true",
        help="Write output in RELION 3.0 format instead of the default 3.1 format.",
    )


def main(args: argparse.Namespace) -> None:
    assert args.o.endswith(".star"), "Output file must be .star file"

    input_ext = os.path.splitext(args.particles)[-1]
    assert input_ext in (
        ".mrcs",
        ".txt",
        ".star",
    ), "Input file must be .mrcs/.txt/.star"

    # Either accept an input star file, or an input .mrcs/.txt with CTF .pkl
    # and an optional pose .pkl file(s)
    ctf = poses = eulers = trans = None
    if input_ext == ".star":
        if args.poses is not None:
            raise ValueError(
                "--poses cannot be specified when input is a starfile "
                "(poses are obtained from starfile)"
            )
        if args.ctf is not None:
            raise ValueError(
                "--ctf cannot be specified when input is a starfile "
                "(ctf information are obtained from starfile)"
            )
    else:
        if not args.ctf:
            raise ValueError("--ctf must be specified when input is not a starfile!")

    particles = ImageSource.from_file(args.particles, lazy=True, datadir=args.datadir)
    logger.info(f"{len(particles)} particles in {args.particles}")

    if args.ctf:
        ctf = utils.load_pkl(args.ctf)
        if ctf.shape[1] != 9:
            raise ValueError(
                f"Incorrect CTF .pkl format "
                f"— expected 9 columns but found {ctf.shape[1]}!"
            )
        if len(particles) != len(ctf):
            raise ValueError(
                f"{len(particles)} != {len(ctf)}, "
                f"Number of particles != number of CTF parameters"
            )
    if args.poses:
        poses = utils.load_pkl(args.poses)
        if len(particles) != len(poses[0]):
            raise ValueError(
                f"{len(particles)} != {len(poses)}, "
                f"Number of particles != number of poses"
            )

    # load the particle filter if given and apply it to the CTF and poses data
    ind = np.arange(particles.n)
    if args.ind:
        ind = utils.load_pkl(args.ind)
        logger.info(f"Filtering to {len(ind)} particles")
        if ctf is not None:
            ctf = ctf[ind]
        if poses is not None:
            poses = (poses[0][ind], poses[1][ind])

    # When the input is already a .star file, we just filter the data table directly
    if input_ext == ".star":
        assert isinstance(particles, StarfileSource)
        df = particles.df.loc[ind]
        optics = None if args.relion30 else particles.data_optics

    # When the input is not a .star file, we have to create the data table fields
    else:
        base_dir = os.path.dirname(os.path.abspath(args.particles))

        if input_ext == ".txt":
            with open(args.particles, "r") as f:
                mrcs_files = f.read().splitlines()

            counts = [
                MRCHeader.parse(os.path.join(base_dir, f)).fields["nz"]
                for f in mrcs_files
            ]
            ind_lbls = np.concatenate([np.arange(count) for count in counts])[ind] + 1
            image_names = np.repeat(mrcs_files, counts)[ind]
        else:
            ind_lbls = [str(i + 1) for i in ind]
            image_names = particles.filenames[ind]

        if args.full_path:
            image_names = [
                os.path.abspath(
                    os.path.join(os.path.dirname(args.particles), image_name)
                )
                for image_name in image_names
            ]

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

        # Figure out what the optics groups are using voltage, spherical aberration,
        # and amplitude contrast to assign a unique group to each particle
        if not args.relion30:
            optics_cols = list(set(range(9)) - ctf_cols)
            optics_headers = [CTF_HEADERS[optics_col] for optics_col in optics_cols]
            optics_groups, optics_indx = np.unique(
                ctf[:, optics_cols], return_inverse=True, axis=0
            )

            data["_rlnOpticsGroup"] = optics_indx + 1
            optics = pd.DataFrame(optics_groups, columns=optics_headers)
            optics["_rlnOpticsGroup"] = np.array(
                [str(i + 1) for i in range(optics.shape[0])],
            )
        else:
            optics = None

        if eulers is not None and trans is not None:
            for i in range(3):
                data[POSE_HDRS[i]] = eulers[:, i]  # type: ignore
            for i in range(2):
                data[POSE_HDRS[3 + i]] = trans[:, i]

        df = pd.DataFrame(data=data)

    write_star(args.o, data=df, data_optics=optics)
