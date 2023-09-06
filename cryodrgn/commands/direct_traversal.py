"""
Construct z path interpolating between anchor points
"""

import argparse
import numpy as np
import os

from cryodrgn import utils


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("z", help="Input z.pkl embeddings")
    parser.add_argument(
        "--ind",
        metavar=".TXT",
        required=True,
        help="Text file containing indices of anchor points",
    )
    parser.add_argument(
        "-n",
        type=int,
        default=6,
        help="Number of points in between anchors, inclusive (default: %(default)s)",
    )
    parser.add_argument("--loop", action="store_true", help="Loop to first point")
    parser.add_argument(
        "-o",
        metavar="Z.PATH.TXT",
        type=os.path.abspath,
        required=True,
        help="Output .txt file for z-values",
    )
    return parser


def main(args):
    assert args.n > 2
    if not os.path.exists(os.path.dirname(args.o)):
        os.makedirs(os.path.dirname(args.o))
    z_all = utils.load_pkl(args.z)
    zdim = z_all.shape[1]
    ind = np.loadtxt(args.ind).astype(int)
    if args.loop:
        ind = np.append(ind, ind[0])
    z_anchors = z_all[ind]
    z_path = []
    for i in range(len(ind) - 1):
        z_start = z_anchors[i]
        z_end = z_anchors[i + 1]
        z = np.repeat(np.arange(args.n, dtype=np.float32), zdim).reshape(args.n, zdim)
        z *= (z_end - z_start) / (args.n - 1)
        z += z_start
        z_path.append(z[:-1])
    z_path.append(z_end.reshape(1, -1))
    z_path = np.concatenate(z_path)
    print(z_path.shape)
    print(args.o)
    np.savetxt(args.o, z_path)


if __name__ == "__main__":
    main(parse_args().parse_args())
