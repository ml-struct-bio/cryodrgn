"""Construct a path of embeddings in latent space along principal components.

Example usages
--------------
$ cryodrgn pc_traversal zvals.pkl
$ cryodrgn pc_traversal zvals.pkl --pc 3
$ cryodrgn pc_traversal zvals.pkl --pc 4 -n 12 ---lim 0.10 0.90 -o z-path-new.txt

"""
import os
import argparse
import pickle
import numpy as np
from scipy.spatial.distance import cdist
from cryodrgn import analysis


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("zfile", help="Input .pkl file containing z-embeddings")
    parser.add_argument(
        "--pc", type=int, nargs="+", help="Choose PCs (1-based indexing) (default: all)"
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
        "--outdir",
        "-o",
        type=os.path.abspath,
        nargs="?",
        const="zpaths",
        metavar="Z-DIR",
        help="output folder for pc<i>.txt files path z-values; "
        "choose name automatically if flag given with no name",
    )


def analyze_data_support(z, traj, cutoff=3):
    d = cdist(traj, z)
    count = (d < cutoff).sum(axis=1)
    return count


def main(args):
    if args.outdir:
        os.makedirs(args.outdir)

    z = pickle.load(open(args.zfile, "rb"))
    zdim = z.shape[1]
    pc, pca = analysis.run_pca(z)
    dims = args.pc if args.pc is not None else list(range(1, zdim + 1))
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

        if args.outdir:
            np.savetxt(os.path.join(args.outdir, f"pc{dim}.txt"), traj)
        else:
            print(traj)
