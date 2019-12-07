'''Select closet point in original latent space based on 2D PCA projection or provided UMAP embedding'''

import argparse
import numpy as np
import sys, os
import pickle

from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('z', help='Input z.pkl')
    parser.add_argument('--umap', help='Select from UMAP embedding (or tSNE)')
    parser.add_argument('--use-pca', action='store_true', help='Select from 2D PCA projection')
    parser.add_argument('--xy', nargs=2, type=float, help='Selected points from 2D space')
    return parser

def pca_transform(z):
    pca = PCA(z.shape[1])
    pca.fit(z)
    print('Explained variance ratio:')
    print(pca.explained_variance_ratio_)
    pc = pca.transform(z)
    return pc, pca

def main(args):
    z = pickle.load(open(args.z,'rb'))
    # Either perform PCA or load the umap embedding
    if args.use_pca:
        umap, _ = pca_transform(z)
    else:
        umap = pickle.load(open(args.umap,'rb'))

    # Compute distance
    tmp = np.linalg.norm(umap - np.array(args.xy), axis=1)
    ind = tmp.argsort()

    # Print out the closest point
    i = ind[0]
    print(i)
    print(tmp[i])
    print(' '.join([str(ii) for ii in z[i].round(5)]))

    # Get 10 closest points -- how much variation between these points?
    zclose = z[ind[:10]]
    print(zclose.std(axis=0))

    plt.scatter(umap[:,0],umap[:,1],s=1,alpha=.05)
    plt.scatter(umap[i,0],umap[i,1],color='k')
    plt.show()

if __name__ == '__main__':
    main(parse_args().parse_args())

