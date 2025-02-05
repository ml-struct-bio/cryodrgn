"""View the first row of a cryosparc .cs file"""

import argparse
import numpy as np


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input", help="Input")
    return parser


def main(args: argparse.Namespace) -> None:
    x = np.load(args.input)
    print(f"{len(x)} particles")
    w = np.max([len(n) for n in x.dtype.names])
    for a, b in zip(x.dtype.names, x[0]):
        print(f"{a:{w}}:    {b}")
