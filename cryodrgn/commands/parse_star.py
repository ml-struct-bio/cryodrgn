"""Parse image CTF and poses from RELION .star file into separate files for cryoDRGN.

This command is often used as a part of preparing inputs for training commands such as
`train_vae` and `abinit_homo` when particles are coming from a .star file.

Example usage
-------------
$ cryodrgn parse_star particles_from_M.star --ctf ctf.pkl --poses pose.pkl

# Override image parameters even if given in file
$ cryodrgn parse_star particles_from_M.star --ctf ctf.pkl --poses pose.pkl \
                                            -D 294 --Apix 1.7

"""
import os
import argparse
from cryodrgn.commands.parse_ctf_star import main as ctf_main
from cryodrgn.commands.parse_pose_star import main as pose_main


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input", help="Input .star file")
    parser.add_argument("--ctf", help="Output pkl of CTF parameters")
    parser.add_argument("--poses", help="Output pkl of poses")

    parser.add_argument(
        "--png", metavar="PNG", type=os.path.abspath, help="Optionally plot the CTF"
    )

    group = parser.add_argument_group("Optionally provide missing image parameters")
    group.add_argument("-D", type=int, help="Image size in pixels")
    group.add_argument("--Apix", type=float, help="Angstroms per pixel")
    group.add_argument("--kv", type=float, help="Accelerating voltage (kV)")
    group.add_argument("--cs", type=float, help="Spherical abberation (mm)")
    group.add_argument("-w", type=float, help="Amplitude contrast ratio")
    group.add_argument("--ps", type=float, help="Phase shift (deg)")

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files if they already exist",
    )


def main(args: argparse.Namespace) -> None:
    if not args.overwrite:
        if os.path.exists(args.ctf):
            raise RuntimeError(
                f"{args.ctf} already exists. Use --overwrite to overwrite."
            )
        if os.path.exists(args.poses):
            raise RuntimeError(
                f"{args.poses} already exists. Use --overwrite to overwrite."
            )

    setattr(args, "star", args.input)
    setattr(args, "o", args.ctf)
    ctf_main(args)

    setattr(args, "outpkl", args.poses)
    pose_main(args)
