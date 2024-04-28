"""Construct a path in z-latent-space interpolating directly between anchor points.

Example usages
--------------
$ cryodrgn direct_traversal zvals.pkl --anchors anchors.txt
$ cryodrgn direct_traversal zvals.pkl --anchors anchors.txt -n 20 -o z-path-new.txt
$ cryodrgn direct_traversal zvals.pkl --anchors anchors.txt -n 3 --loop -o

"""
import os
import argparse
import numpy as np
from cryodrgn import utils


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("zfile", help="Input .pkl file containing z-embeddings")
    parser.add_argument(
        "--anchors",
        required=True,
        nargs="+",
        help="List of anchor point indices in the given file, either given "
        "directly as integers, or in a .txt file(s).",
    )
    parser.add_argument(
        "-n",
        type=int,
        default=6,
        help="Number of points in between anchors, inclusive (default: %(default)s)",
    )
    parser.add_argument("--loop", action="store_true", help="Loop to first point")
    parser.add_argument(
        "--outtxt",
        "-o",
        type=os.path.abspath,
        nargs="?",
        const="z-path.txt",
        metavar="Z-PATH.TXT",
        help="output .txt file for path z-values; "
        "choose name automatically if flag given with no name",
    )


def parse_anchors(
    given_anchors: list[str], zvals: np.ndarray, zfile: str, loop: bool = False
) -> list[int]:
    anchors = list()

    for anchor_txt in given_anchors:
        if os.path.exists(anchor_txt):
            new_anchors = np.loadtxt(anchor_txt).astype(int).tolist()
        elif anchor_txt.isnumeric():
            new_anchors = [int(anchor_txt)]
        else:
            raise ValueError(
                f"Unrecognized anchor value `{anchor_txt}` which is neither an integer "
                f"nor a .txt file containing a list of integers!"
            )

        for new_anchor in new_anchors:
            if new_anchor < 0:
                raise ValueError(
                    f"Invalid anchor index {new_anchor} is not a positive integer!"
                )
            if new_anchor >= zvals.shape[0]:
                raise ValueError(
                    f"Invalid anchor index {new_anchor} too big for "
                    f"`{zfile}` containing {zvals.shape[0]} points!"
                )
        anchors += new_anchors

    if loop:
        anchors.append(anchors[0])

    if len(anchors) < 2:
        raise ValueError(
            f"Need at least two anchors for graph traversal; given {len(anchors)}!"
        )

    return anchors


def main(args: argparse.Namespace) -> None:
    z_all = utils.load_pkl(args.zfile)
    zdim = z_all.shape[1]
    ind = parse_anchors(args.anchors, z_all, args.zfile, args.loop)
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

    if args.outtxt:
        if not os.path.exists(os.path.dirname(args.outtxt)):
            os.makedirs(os.path.dirname(args.outtxt))
        np.savetxt(args.outtxt, z_path)
    else:
        print(z_path)
