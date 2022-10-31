"""View the first row of a cryosparc .cs file"""

import argparse

import numpy as np


def add_args(parser):
    parser.add_argument("input", help="Input")
    return parser


def main(args):
    x = np.load(args.input)
    print(f"{len(x)} particles")
    w = np.max([len(n) for n in x.dtype.names])
    for a, b in zip(x.dtype.names, x[0]):
        print(f"{a:{w}}:    {b}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
