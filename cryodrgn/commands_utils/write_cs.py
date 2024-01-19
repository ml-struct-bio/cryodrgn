"""
Create a CryoSparc .cs file from a particle stack and ctf parameters, or an input .cs file
"""

import argparse
import os
import numpy as np
import logging
from cryodrgn import utils
from cryodrgn.source import ImageSource

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument("particles", help="Input particles (.cs)")
    parser.add_argument("--datadir", help="Data/Project directory for cryoSPARC")
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


def main(args):
    assert args.o.endswith(".cs"), "Output file must be .cs file"
    input_ext = os.path.splitext(args.particles)[-1]
    assert input_ext in (".mrcs", ".txt", ".cs"), "Input file must be .mrcs/.txt/.cs"

    # Either accept an input cs file, or an input .mrcs/.txt with optional ctf/pose pkl file(s)
    particles = ImageSource.from_file(args.particles, lazy=True, datadir=args.datadir)
    if input_ext == ".cs":
        assert (
            args.poses is None
        ), "--poses cannot be specified when input is a cs file (poses are obtained from cs file)"
        assert (
            args.ctf is None
        ), "--ctf cannot be specified when input is a cs file (ctf information are obtained from cs file)"
    else:
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
    logger.info(f"{len(particles)} particles in {args.particles}")

    ind = None
    if args.ind:
        ind = utils.load_pkl(args.ind)
        logger.info(f"Filtering to {len(ind)} particles")
    particles = np.array(particles.images(ind))

    if input_ext == ".cs":
        pass  # Nothing to be done - we've already sub-setted the .cs data
    else:
        # Carefully construct a np array with all fields (ctf/pose/others) filled in
        raise NotImplementedError(
            "Generation of a .cs file from input .mrcs/.txt file is coming soon"
        )

    assert isinstance(particles, np.ndarray)
    with open(args.o, "wb") as f:
        np.save(f, particles)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
