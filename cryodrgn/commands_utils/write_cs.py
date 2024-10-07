"""Create a CryoSparc .cs file from a particle stack, using poses and CTF if necessary.

Example usage
-------------
$ cryodrgn_utils write_cs particles.mrcs --poses pose.pkl --ctf ctf.pkl -o particles.cs
$ cryodrgn_utils write_cs particles.star --datadir=/scratch/empiar_10345/Micrographs \
                          -o particles.cs

"""
import argparse
import os
import logging
from cryodrgn.commands_utils.filter_cs import main as filter_cs_main

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
    logger.warning(
        "`cryodrgn write_cs` is deprecated as of cryoDRGN v3.4.1 and will be removed "
        "in a future version; use `filter_cs` instead to filter .cs files!",
    )
    filter_cs_main(args)
