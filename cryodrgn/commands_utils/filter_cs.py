"""Create a CryoSparc .cs file from a particle stack, using poses and CTF if necessary.

Example usage
-------------
$ cryodrgn_utils filter_cs particles_exported.cs --ind good-particles.pkl \
                           -o filtered-particles.cs

"""
import argparse
import os
import numpy as np
import logging
from cryodrgn import utils

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("particles", help="Input particles (.cs)")
    parser.add_argument(
        "--ind",
        help="Optionally filter by array of selected indices (.pkl)",
        required=True,
    )
    parser.add_argument(
        "-o", "--output", type=os.path.abspath, required=True, help="Output .star file"
    )


def main(args: argparse.Namespace) -> None:
    input_ext = os.path.splitext(args.particles)[-1]
    if input_ext != ".cs":
        raise ValueError(f"Input must be .cs, given {input_ext} file instead!")
    if not args.output.endswith(".cs"):
        raise ValueError(f"Output file {args.output} is not a .cs file!")
    ind = utils.load_pkl(args.ind)

    particles = np.load(args.particles)
    logger.info(f"{len(particles)} particles in {args.particles}")
    logger.info(f"Filtering to {len(ind)} particles")
    particles = particles[ind]

    assert isinstance(particles, np.ndarray)
    with open(args.output, "wb") as f:
        np.save(f, particles)
