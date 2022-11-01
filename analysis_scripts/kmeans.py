"""K-means clustering"""

import argparse
import pickle

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA

from cryodrgn import analysis


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Input z.pkl")
    parser.add_argument("-k", type=int, required=True, help="# clusters")
    parser.add_argument("--stride", type=int, help="Stride the dataset")
    parser.add_argument("-o", help="Output labels (.pkl)")
    parser.add_argument("--out-png", help="Output image (.png)")
    parser.add_argument("--out-k", help="Output cluster centers z values (.txt)")
    parser.add_argument(
        "--on-data",
        action="store_true",
        help="Use nearest data point instead of cluster center",
    )
    parser.add_argument("--out-k-ind", help="Output cluster center indices (.txt)")
    parser.add_argument(
        "--reorder", action="store_true", help="Reorder cluster centers"
    )
    return parser


def main(args):
    fig, ax = plt.subplots()
    print(args)
    z = pickle.load(open(args.input, "rb"))
    if args.stride:
        z = z[:: args.stride]
    print("{} points".format(len(z)))

    # k-means clustering
    labels, centers = analysis.cluster_kmeans(
        z, args.k, on_data=args.on_data, reorder=args.reorder
    )

    # use the nearest data point instead of cluster centroid
    if args.on_data:
        centers_zi = cdist(centers, z).argmin(axis=1)
        print(centers_zi)
        centers_z = z[centers_zi]
        centers = centers_z
        if args.out_k_ind:
            np.savetxt(args.out_k_ind, centers_zi, fmt="%d")

    if args.o:
        with open(args.o, "wb") as f:
            pickle.dump(labels, f)

    if args.out_k:
        np.savetxt(args.out_k, centers)

    # dimensionality reduction for viz
    pca = PCA(z.shape[1])
    pca.fit(z)
    print("PCA explained variance ratio:")
    print(pca.explained_variance_ratio_)
    pc = pca.transform(z)

    for i in range(args.k):
        ii = np.where(labels == i)
        pc_sub = pc[ii]
        plt.scatter(
            pc_sub[:, 0], pc_sub[:, 1], s=2, alpha=0.1, label="cluster {}".format(i)
        )

    c = pca.transform(centers)
    plt.scatter(c[:, 0], c[:, 1], c="k")
    for i in range(args.k):
        ax.annotate(str(i), c[i, 0:2])

    xx, yy = 0, 1
    plt.xlabel("PC{} ({:3f})".format(xx + 1, pca.explained_variance_ratio_[xx]))
    plt.ylabel("PC{} ({:3f})".format(yy + 1, pca.explained_variance_ratio_[yy]))

    if args.out_png:
        plt.savefig(args.out_png)
    else:
        plt.show()


if __name__ == "__main__":
    main(parse_args().parse_args())
