"""Parse image poses from a cryoSPARC .cs metafile"""

import argparse
import os
import pickle
import logging
import numpy as np
import torch
from cryodrgn.models import lie_tools

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input", help="Cryosparc .cs file")
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
    parser.add_argument(
        "-D", type=int, required=True, help="Box size of reconstruction (pixels)"
    )
    parser.add_argument(
        "-o", metavar="PKL", type=os.path.abspath, required=True, help="Output pose.pkl"
    )
    return parser


def main(args: argparse.Namespace) -> None:
    assert args.input.endswith(".cs"), "Input format must be .cs file"
    assert args.o.endswith(".pkl"), "Output format must be .pkl"

    data = np.load(args.input)
    # view the first row
    for i in range(len(data.dtype)):
        print(i, data.dtype.names[i], data[0][i])

    if args.abinit:
        RKEY = "alignments_class_0/pose"
        TKEY = "alignments_class_0/shift"
    else:
        RKEY = "alignments3D/pose"
        TKEY = "alignments3D/shift"

    # parse rotations
    logger.info(f"Extracting rotations from {RKEY}")
    rot = np.array([x[RKEY] for x in data])
    rot = torch.tensor(rot)
    rot = lie_tools.expmap(rot)
    rot = rot.cpu().numpy()
    logger.info("Transposing rotation matrix")
    rot = np.array([x.T for x in rot])
    logger.info(rot.shape)

    # parse translations
    logger.info(f"Extracting translations from {TKEY}")
    trans = np.array([x[TKEY] for x in data])
    if args.hetrefine:
        logger.info("Scaling shifts by 2x")
        trans *= 2
    logger.info(trans.shape)

    # convert translations from pixels to fraction
    trans /= args.D

    # write output
    logger.info(f"Writing {args.o}")
    with open(args.o, "wb") as f:
        pickle.dump((rot, trans), f)
