"""Create a CryoSparc .cs file from a particle stack, using poses and CTF if necessary.

Example usage
-------------
$ cryodrgn_utils filter_cs particles.mrcs --poses pose.pkl --ctf ctf.pkl -o particles.cs
$ cryodrgn_utils write_cs particles.star --datadir=/scratch/empiar_10345/Micrographs \
                          -o particles.cs

"""
import argparse
import os
import numpy as np
import logging
from cryodrgn import utils
from cryodrgn.source import ImageSource

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
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
        "-o", "--output", type=os.path.abspath, required=True, help="Output .star file"
    )


def main(args: argparse.Namespace) -> None:
    input_ext = os.path.splitext(args.particles)[-1]
    if input_ext not in (".mrcs", ".txt", ".cs", ".star"):
        raise ValueError(
            f"Input must be .mrcs/.txt/.star/.cs, given {input_ext} file instead!"
        )
    if not args.output.endswith(".cs"):
        raise ValueError(f"Output file {args.output} is not a .cs file!")

    # We will filter particles using the set of indices if they are given
    ind = utils.load_pkl(args.ind) if args.ind is not None else None

    # Either accept an input cs file...
    if input_ext == ".cs":
        if args.poses is not None:
            raise ValueError(
                "--poses cannot be specified when input is a cs file (poses are "
                "obtained from cs file)!"
            )
        if args.ctf is not None:
            raise ValueError(
                "--ctf cannot be specified when input is a cs file (ctf "
                "information is obtained from cs file)!"
            )

        particles = np.load(args.particles)
        logger.info(f"{len(particles)} particles in {args.particles}")
        if args.ind:
            logger.info(f"Filtering to {len(ind)} particles")
            particles = particles[ind]

    # ...or an input .mrcs/.txt/.star with optional ctf/pose pkl file(s)
    else:
        particles = ImageSource.from_file(
            args.particles, lazy=True, datadir=args.datadir
        )
        logger.info(f"{particles.orig_n} particles in {args.particles}")

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

        if args.ind:
            logger.info(f"Filtering to {len(ind)} particles")
            particles = np.array(particles.images(ind))
        else:
            particles = np.array(particles.images())

    if input_ext == ".cs":
        pass  # Nothing to be done - we've already sub-setted the .cs data
    else:
        # Carefully construct a np array with all fields (ctf/pose/others) filled in
        raise NotImplementedError(
            "Generation of a .cs file from input .mrcs/.txt file is coming soon"
        )

    assert isinstance(particles, np.ndarray)
    with open(args.output, "wb") as f:
        np.save(f, particles)
