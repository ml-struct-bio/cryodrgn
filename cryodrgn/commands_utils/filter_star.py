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


def main(args):
    s = starfile.Starfile.load(args.input)
    ind = utils.load_pkl(args.ind)

    if args.et:
        particles_to_tilts, _ = TiltSeriesData.parse_particle_tilt(args.input)
        tilt_indices = TiltSeriesData.particles_to_tilts(particles_to_tilts, ind)
        df = s.df.loc[tilt_indices]
    else:
        df = s.df.loc[ind]

    # filter data optics table by what optics groups are left in the particle table
    if s.data_optics is not None:
        if "_rlnOpticsGroup" in df.columns:
            new_grps = set(df["_rlnOpticsGroup"])
            new_optics_df = s.data_optics.df.copy()
            new_optics_df = new_optics_df.loc[
                new_optics_df._rlnOpticsGroup.isin(new_grps)
            ]
            new_optics = starfile.Starfile(headers=None, df=new_optics_df)
        else:
            new_optics = s.data_optics
    else:
        new_optics = None

    if args.micrograph_files:
        if "_rlnMicrographName" not in df.columns:
            raise ValueError(
                "Cannot write micrograph files for a .star file "
                "without a `_rlnMicrographName` field!"
            )

        os.makedirs(args.o, exist_ok=True)
        for micrograph_name, group_df in df.groupby("_rlnMicrographName"):
            filename_without_extension = os.path.splitext(micrograph_name)[0]
            output_path = os.path.join(args.o, f"{filename_without_extension}.star")

            # filter data optics table by what optics groups are left
            # in this micrograph's particle table
            if s.data_optics is not None:
                if "_rlnOpticsGroup" in df.columns:
                    micro_grps = set(group_df["_rlnOpticsGroup"])
                    micro_optics_df = new_optics.df.copy()
                    micro_optics_df = micro_optics_df.loc[
                        micro_optics_df._rlnOpticsGroup.isin(micro_grps)
                    ]
                    micro_optics = starfile.Starfile(headers=None, df=micro_optics_df)
                else:
                    micro_optics = new_optics
            else:
                micro_optics = None

            micro_star = starfile.Starfile(
                headers=None,
                df=group_df,
                relion31=s.relion31,
                data_optics=micro_optics,
            )
            micro_star.write(output_path)
            logger.info(
                f"Wrote .star file for {filename_without_extension} to {output_path}"
            )
    else:
        s = starfile.Starfile(
            headers=None, df=df, relion31=s.relion31, data_optics=new_optics
        )
        s.write(args.o)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
