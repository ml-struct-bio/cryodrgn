"""
Filter a .star file
"""

import argparse
import os
import logging
from cryodrgn import starfile
from cryodrgn import utils
from cryodrgn.dataset import TiltSeriesData

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument("input", help="Input .star file")

    parser.add_argument("--ind", required=True, help="Array of selected indices (.pkl)")
    parser.add_argument(
        "--et", action="store_true", help="Set if .star file includes tilts"
    )
    parser.add_argument("-o", type=os.path.abspath, help="Output .star file")


def main(args):
    s = starfile.Starfile.load(args.input)
    ind = utils.load_pkl(args.ind)

    if args.et:
        particles_to_tilts, tilts_to_particles = TiltSeriesData.parse_particle_tilt(
            args.input
        )

        ind = utils.load_pkl(args.ind)
        tilt_indices = TiltSeriesData.particles_to_tilts(particles_to_tilts, ind)
        df = s.df.loc[tilt_indices]

    else:
        df = s.df.loc[ind]

    s = starfile.Starfile(headers=None, df=df)
    s.write(args.o)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
