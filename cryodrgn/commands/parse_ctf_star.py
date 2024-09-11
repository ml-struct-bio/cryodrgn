"""Parse contrast transfer function values from a RELION .star file into separate file.

This command is often used as a part of preparing inputs for training commands such as
`train_vae` and `abinit_homo` when particles are coming from a .star file.

Example usage
-------------
$ cryodrgn parse_ctf_star particles_from_M.star -o ctf.pkl -D 294 --Apix 1.7

"""
import argparse
import os
import pickle
import logging
import numpy as np
from cryodrgn import ctf
from cryodrgn.starfile import Starfile

logger = logging.getLogger(__name__)

HEADERS = [
    "_rlnDefocusU",
    "_rlnDefocusV",
    "_rlnDefocusAngle",
    "_rlnVoltage",
    "_rlnSphericalAberration",
    "_rlnAmplitudeContrast",
    "_rlnPhaseShift",
]


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("star", help="Input")

    parser.add_argument(
        "-o", type=os.path.abspath, required=True, help="Output pkl of CTF parameters"
    )
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


def main(args: argparse.Namespace) -> None:
    assert args.o.endswith(".pkl"), "Output CTF parameters must be .pkl file"
    stardata = Starfile(args.star)
    logger.info(f"{len(stardata)} particles")

    apix = stardata.apix
    if apix is None:
        if args.Apix is None:
            raise ValueError(
                f"Cannot find A/px values in {args.star} "
                f"— must be given manually with --Apix <val> !"
            )
    if args.Apix is not None:
        apix = args.Apix

    resolution = stardata.resolution
    if resolution is None:
        if args.D is None:
            raise ValueError(
                f"Cannot find image size values in {args.star} "
                f"— must be given manually with -D <val> !"
            )
    if args.D is not None:
        resolution = args.D

    # Sometimes CTF parameters are missing from the star file
    overrides = dict()
    if args.kv is not None:
        logger.info(f"Overriding accerlating voltage with {args.kv} kV")
        overrides[HEADERS[3]] = args.kv
    if args.cs is not None:
        logger.info(f"Overriding spherical abberation with {args.cs} mm")
        overrides[HEADERS[4]] = args.cs
    if args.w is not None:
        logger.info(f"Overriding amplitude contrast ratio with {args.w}")
        overrides[HEADERS[5]] = args.w
    if args.ps is not None:
        logger.info(f"Overriding phase shift with {args.ps}")
        overrides[HEADERS[6]] = args.ps

    ctf_params = np.zeros((len(stardata), 9))
    ctf_params[:, 0] = resolution
    ctf_params[:, 1] = apix
    for i, header in enumerate(
        [
            "_rlnDefocusU",
            "_rlnDefocusV",
            "_rlnDefocusAngle",
            "_rlnVoltage",
            "_rlnSphericalAberration",
            "_rlnAmplitudeContrast",
            "_rlnPhaseShift",
        ]
    ):
        ctf_params[:, i + 2] = (
            stardata.get_optics_values(header)
            if header not in overrides
            else overrides[header]
        )

    logger.info("CTF parameters for first particle:")
    ctf.print_ctf_params(ctf_params[0])
    logger.info("Saving {}".format(args.o))

    with open(args.o, "wb") as f:
        pickle.dump(ctf_params.astype(np.float32), f)

    if args.png:
        import matplotlib.pyplot as plt

        assert args.D, "Need image size to plot CTF"
        ctf.plot_ctf(args.D, args.Apix, ctf_params[0, 2:])
        plt.savefig(args.png)
        logger.info(args.png)
