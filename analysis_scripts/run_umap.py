"""
UMAP dimensionality reduction
"""

import argparse
import pickle
import warnings

import matplotlib.pyplot as plt
import umap

warnings.filterwarnings("ignore")  # ignore numba warnings from umap


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Input z.pkl")
    parser.add_argument("--stride", type=int, help="Stride the dataset")
    parser.add_argument("-o", help="Output UMAP embeddings (.pkl)")
    parser.add_argument("--show", action="store_true", help="Show UMAP plot")
    return parser


def main(args):
    z = pickle.load(open(args.input, "rb"))
    if args.stride:
        z = z[:: args.stride]
    print(z.shape)
    reducer = umap.UMAP()
    z_embedded = reducer.fit_transform(z)
    if args.o:
        pickle.dump(z_embedded, open(args.o, "wb"))
    if args.show:
        plt.scatter(z_embedded[:, 0], z_embedded[:, 1], s=2, alpha=0.05)
        plt.show()


if __name__ == "__main__":
    main(parse_args().parse_args())
