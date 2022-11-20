"""
Generate trajectory along PCs
"""

import argparse
import os
import pickle

import numpy as np
from scipy.spatial.distance import cdist

from cryodrgn import analysis


def add_args(parser):
    parser.add_argument("z", help="Input z.pkl")
    parser.add_argument(
        "--dim", type=int, help="Choose PC (1-based indexing) (default: all)"
    )
    parser.add_argument(
        "-n",
        type=int,
        default=10,
        help="Number of samples along PC (default: %(default)s)",
    )
    parser.add_argument(
        "--lim",
        nargs=2,
        type=float,
        help="Start and end point of trajectory (default: 5/95th percentile)",
    )
    parser.add_argument(
        "--use-percentile-spacing",
        action="store_true",
        help="Use equally spaced percentiles of the distribution instead of equally spaced points along the PC",
    )
    parser.add_argument(
        "-o",
        type=os.path.abspath,
        required=True,
        help="Output directory for pc.X.txt files",
    )
    return parser


def analyze_data_support(z, traj, cutoff=3):
    d = cdist(traj, z)
    count = (d < cutoff).sum(axis=1)
    return count


def main(args):
    if not os.path.exists(args.o):
        os.makedirs(args.o)

    z = pickle.load(open(args.z, "rb"))
    zdim = z.shape[1]
    pc, pca = analysis.run_pca(z)

    # Use 1-based indexing
    dims = [args.dim] if args.dim else list(range(1, zdim + 1))
    lim = args.lim if args.lim else (5, 95)

    for dim in dims:
        print("PC{}".format(dim))
        if args.use_percentile_spacing:
            pc_values = np.percentile(
                pc[:, dim - 1], np.linspace(lim[0], lim[1], args.n)
            )
            print("Limits: {}, {}".format(pc_values[0], pc_values[-1]))
            traj = analysis.get_pc_traj(pca, zdim, args.n, dim, None, None, pc_values)
        else:
            start = np.percentile(pc[:, dim - 1], lim[0])
            stop = np.percentile(pc[:, dim - 1], lim[1])
            print("Limits: {}, {}".format(start, stop))
            traj = analysis.get_pc_traj(
                pca, zdim, args.n, dim, float(start), float(stop)
            )

        print("Neighbor count along trajectory:")
        print(analyze_data_support(z, traj))

        out = f"{args.o}/pc{dim}.txt"
        print(out)
        np.savetxt(out, traj)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
