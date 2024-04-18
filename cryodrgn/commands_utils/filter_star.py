"""Filter a .star file using a saved set of particle indices.

Example usages
--------------
$ cryodrgn_utils filter_star particles.star --ind good-particles.pkl \
                             -o filtered-particles.star
$ cryodrgn_utils filter_star tilts.star --et --ind good-particles.pkl \
                             -o filtered-tilts.star

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
    parser.add_argument("--tomo", action="store_true", help="Split output by micrograph name into separate .star files")
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

    if args.tomo:
        df = s.df.loc[ind]
        os.makedirs(args.o, exist_ok=True)  # Ensure the output directory exists
        for micrograph_name, group_df in df.groupby('_rlnMicrographName'):
            filename_without_extension = os.path.splitext(micrograph_name)[0]  # Remove extension
            output_path = os.path.join(args.o, f"{filename_without_extension}.star")
            new_star = starfile.Starfile(headers=None, df=group_df)
            new_star.write(output_path)
            logger.info(f"Written .star file for {filename_without_extension} to {output_path}")
        exit()
    else:
        df = s.df.loc[ind]

    s = starfile.Starfile(headers=None, df=df)
    s.write(args.o)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())

