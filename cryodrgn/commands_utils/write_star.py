"""
Create a Relion 3.0 star file from a particle stack and ctf parameters
"""

import argparse
import os

import numpy as np
import pandas as pd

from cryodrgn import dataset, mrc, utils
from cryodrgn.starfile import Starfile

log = utils.log

CTF_HEADERS = [
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
    parser.add_argument("--ctf", help="Input ctf.pkl")
    parser.add_argument("--poses", help="Optionally include pose.pkl")
    parser.add_argument(
        "--ind", help="Optionally filter by array of selected indices (.pkl)"
    )
    parser.add_argument(
        "--full-path",
        action="store_true",
        help="Write the full path to particles (default: relative paths)",
    )
    parser.add_argument(
        "-o", type=os.path.abspath, required=True, help="Output .star file"
    )

    return parser


def parse_chunk_size(txtfile):
    lines = open(txtfile, "r").readlines()

    def abspath(f):
        if os.path.isabs(f):
            return f
        base = os.path.dirname(os.path.abspath(txtfile))
        return os.path.join(base, f)

    lines = [abspath(x).strip() for x in lines]
    return [mrc.parse_header(f).fields["nz"] for f in lines]


def main(args):
    assert args.o.endswith(".star"), "Output file must be .star file"
    input_ext = os.path.splitext(args.particles)[-1]
    assert input_ext in (
        ".mrcs",
        ".txt",
        ".star",
    ), "Input file must be .mrcs/.txt/.star"

    # Either accept an input star file, or an input .mrcs/.txt with optional ctf/pose pkl file(s)
    starfile = None
    if input_ext == ".star":
        assert (
            args.poses is None
        ), "--poses cannot be specified when input is a starfile (poses are obtained from starfile)"
        assert (
            args.ctf is None
        ), "--ctf cannot be specified when input is a starfile (ctf information are obtained from starfile)"

        starfile = Starfile.load(args.particles)
        particles = starfile.get_particles(
            datadir=os.path.dirname(args.particles), lazy=True
        )

    else:
        particles = dataset.load_particles(args.particles, lazy=True)
        if args.ctf:
            ctf = utils.load_pkl(args.ctf)
            assert ctf.shape[1] == 9, "Incorrect CTF pkl format"
            assert len(particles) == len(
                ctf
            ), f"{len(particles)} != {len(ctf)}, Number of particles != number of CTF paraameters"
        if args.poses:
            poses = utils.load_pkl(args.poses)
            assert len(particles) == len(
                poses[0]
            ), f"{len(particles)} != {len(poses)}, Number of particles != number of poses"
    log(f"{len(particles)} particles in {args.particles}")

    if input_ext == ".star":
        particle_ind = np.arange(len(particles))
    # Get index for particles in each .mrcs file
    elif input_ext == ".txt":
        N_per_chunk = parse_chunk_size(args.particles)
        particle_ind = np.concatenate([np.arange(nn) for nn in N_per_chunk])
        assert len(particle_ind) == len(particles)
    else:  # single .mrcs file
        particle_ind = np.arange(len(particles))

    if args.ind:
        ind = utils.load_pkl(args.ind)
        log(f"Filtering to {len(ind)} particles")
        particles = [particles[ii] for ii in ind]
        if args.ctf:
            ctf = ctf[ind]
        if args.poses:
            poses = (poses[0][ind], poses[1][ind])
        particle_ind = particle_ind[ind]

    if input_ext == ".star":
        df = starfile.df.loc[particle_ind]
    else:
        particle_ind += 1  # CHANGE TO 1-BASED INDEXING
        image_names = [img.fname for img in particles]
        if args.full_path:
            image_names = [os.path.abspath(img.fname) for img in particles]
        names = [f"{i}@{name}" for i, name in zip(particle_ind, image_names)]

        if args.ctf:
            ctf = ctf[:, 2:]

        # convert poses
        if args.poses:
            eulers = utils.R_to_relion_scipy(poses[0])
            D = particles[0].get().shape[0]
            trans = poses[1] * D  # convert from fraction to pixels

        # Create a new dataframe with required star file headers
        data = {"_rlnImageName": names}
        if args.ctf:
            for i in range(7):
                data[CTF_HEADERS[i]] = ctf[:, i]

        if args.poses:
            for i in range(3):
                data[POSE_HDRS[i]] = eulers[:, i]
            for i in range(2):
                data[POSE_HDRS[3 + i]] = trans[:, i]
        df = pd.DataFrame(data=data)

    s = Starfile(headers=None, df=df)
    s.write(args.o)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
