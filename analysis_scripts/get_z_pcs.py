'''
Generate trrajectory along PCs
'''

import argparse
import numpy as np
import sys, os
import pickle
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('z', help='Input z.pkl')
    parser.add_argument('--dim', type=int, help='Choose PC (default: all)')
    parser.add_argument('-n', type=int, default=10, help='Number of samples along PC')
    parser.add_argument('--lim', nargs=2, type=float, help='Start and end point of trajectory (default: 5/95th percentile)')
    parser.add_argument('-o', help='Output txt file of z values or prefix to output .pc.txt files')
    return parser

def pca_transform(z):
    pca = PCA(z.shape[1])
    pca.fit(z)
    print('Explained variance ratio:')
    print(pca.explained_variance_ratio_)
    pc = pca.transform(z)
    return pc, pca

def get_pc_traj(pca, D, numpoints, dim, start, end):
    traj_pca = np.zeros((numpoints,D))
    traj_pca[:,dim] = np.linspace(start, end, numpoints)
    ztraj_pca = pca.inverse_transform(traj_pca)
    return ztraj_pca

def analyze_data_support(z, traj, cutoff=3):
    d = cdist(traj,z)
    count = (d < cutoff).sum(axis=1)
    return count

def main(args):
    z = pickle.load(open(args.z,'rb'))
    D = z.shape[1]
    pc, pca = pca_transform(z)

    if args.dim:
        dim = args.dim
        print('PC{}'.format(dim+1))
        if args.lim:
            start, stop = args.lim
        else:
            start = np.percentile(pc[:,dim], 5)
            stop = np.percentile(pc[:,dim], 95)
        print('Limits: {}, {}'.format(start, stop))
        traj = get_pc_traj(pca, D, args.n, args.dim, start, stop)
        print('Neighbor count along trajectory:')
        print(analyze_data_support(z, traj))
        print(args.o)
        np.savetxt(args.o, traj)
    else:
        for dim in range(D):
            print('PC{}'.format(dim+1))
            start = np.percentile(pc[:,dim], 5)
            stop = np.percentile(pc[:,dim], 95)
            print('Limits: {}, {}'.format(start, stop))
            traj = get_pc_traj(pca, D, args.n, dim, start, stop)
            print('Neighbor count along trajectory:')
            print(analyze_data_support(z, traj))
            out = '{}.pc{}.txt'.format(args.o,dim+1)
            print(out)
            np.savetxt(out, traj)

if __name__ == '__main__':
    main(parse_args().parse_args())

