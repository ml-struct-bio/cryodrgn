"""Parse CTF parameters from a RELION .star file"""

import argparse
import os
import pickle
import logging
import numpy as np
from cryodrgn import ctf, starfile

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


def add_args(parser):
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
    return parser


def main(args):
    assert args.star.endswith(".star"), "Input file must be .star file"
    assert args.o.endswith(".pkl"), "Output CTF parameters must be .pkl file"
    s = starfile.Starfile.load(args.star)
    N = len(s.df)
    logger.info(f"{N} particles")

    overrides = {}
    if s.relion31:
        assert s.data_optics is not None
        df = s.data_optics.df
        assert len(df) == 1, "Only one optics group supported"
        args.D = int(df["_rlnImageSize"][0])
        args.Apix = float(df["_rlnImagePixelSize"][0])
        overrides[HEADERS[3]] = float(df[HEADERS[3]][0])
        overrides[HEADERS[4]] = float(df[HEADERS[4]][0])
        overrides[HEADERS[5]] = float(df[HEADERS[5]][0])
    else:
        assert args.D is not None, "Must provide image size with -D"
        assert args.Apix is not None, "Must provide pixel size with --Apix"

    # Sometimes CTF parameters are missing from the star file
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

    ctf_params = np.zeros((N, 9))
    ctf_params[:, 0] = args.D
    ctf_params[:, 1] = args.Apix
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
            s.df[header] if header not in overrides else overrides[header]
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    main(add_args(parser).parse_args())
