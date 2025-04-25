"""Parse image CTF and poses from a cryoSPARC .cs metafile"""

import os
import argparse
from cryodrgn.commands.parse_ctf_csparc import main as ctf_main
from cryodrgn.commands.parse_pose_csparc import main as pose_main


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input", help="Cryosparc .cs file")
    parser.add_argument("--ctf", help="Output pkl of CTF parameters")
    parser.add_argument("--poses", help="Output pkl of poses")
    parser.add_argument(
        "--png", metavar="PNG", type=os.path.abspath, help="Optionally plot the CTF"
    )

    parser.add_argument(
        "--abinit",
        action="store_true",
        help="Flag if results are from ab-initio reconstruction",
    )
    parser.add_argument(
        "--hetrefine",
        action="store_true",
        help="Flag if results are from a heterogeneous refinements (default: homogeneous refinement)",
    )

    group = parser.add_argument_group("Optionally provide missing image parameters")
    group.add_argument("-D", type=int, help="Image size in pixels")
    group.add_argument("--Apix", type=float, help="Angstroms per pixel")


def main(args: argparse.Namespace) -> None:
    setattr(args, "cs", args.input)
    setattr(args, "o", args.ctf)
    ctf_main(args)

    setattr(args, "o", args.poses)
    pose_main(args)
