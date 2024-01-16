"""Filter a .star file using a saved set of particle indices.

Example usage
-------------
$ cryodrgn_utils filter_star particles.star --ind good-particles.pkl \
                             -o filtered-particles.star
$ cryodrgn_utils filter_star tilts.star --et --ind good-particles.pkl \
                             -o filtered-tilts.star

"""
import argparse
import os
import logging
from cryodrgn import utils
from cryodrgn.dataset import TiltSeriesData
from cryodrgn.starfile import parse_star, write_star

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser):
    parser.add_argument("input", help="Input .star file")

    parser.add_argument("--ind", required=True, help="Array of selected indices (.pkl)")
    parser.add_argument(
        "--et", action="store_true", help="Set if .star file includes tilts"
    )
    parser.add_argument(
        "--micrograph-files",
        "-m",
        action="store_true",
        help="Split output by micrograph name into separate .star files",
    )
    parser.add_argument(
        "-o",
        type=os.path.abspath,
        help="Output .star file (or directory if --micrograph-files)",
    )


def main(args: argparse.Namespace):
    stardf, data_optics = parse_star(args.input)
    ind = utils.load_pkl(args.ind)

    if args.et:
        particles_to_tilts, _ = TiltSeriesData.parse_particle_tilt(args.input)
        tilt_indices = TiltSeriesData.particles_to_tilts(particles_to_tilts, ind)
        filtered_df = stardf.loc[tilt_indices]
    else:
        filtered_df = stardf.loc[ind]

    # filter data optics table by what optics groups are left in the particle table
    if data_optics is not None and "_rlnOpticsGroup" in filtered_df.columns:
        new_grps = set(filtered_df["_rlnOpticsGroup"])
        new_optics = data_optics.loc[data_optics._rlnOpticsGroup.isin(new_grps)]
    else:
        new_optics = None

    if args.micrograph_files:
        if "_rlnMicrographName" not in filtered_df.columns:
            raise ValueError(
                "Cannot write micrograph files for a .star file "
                "without a `_rlnMicrographName` field!"
            )

        os.makedirs(args.o, exist_ok=True)
        for micrograph_name, group_df in filtered_df.groupby("_rlnMicrographName"):
            filename_without_extension = os.path.splitext(micrograph_name)[0]
            output_path = os.path.join(args.o, f"{filename_without_extension}.star")

            # filter data optics table by what optics groups are left
            # in this micrograph's particle table
            if data_optics is not None and "_rlnOpticsGroup" in filtered_df.columns:
                micro_grps = set(group_df["_rlnOpticsGroup"])
                micro_optics = new_optics.loc[
                    new_optics._rlnOpticsGroup.isin(micro_grps)
                ]
            else:
                micro_optics = None

            write_star(output_path, data=group_df, data_optics=micro_optics)
            logger.info(
                f"Wrote .star file for {filename_without_extension} to {output_path}"
            )
    else:
        write_star(args.o, data=filtered_df, data_optics=new_optics)
