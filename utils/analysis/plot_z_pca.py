'''
'''

import argparse
import numpy as np
import sys, os
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
from sklearn.decomposition import PCA

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',  help='Input z pkl')
    parser.add_argument('-o', '--out-png', help='Output PNG')
    parser.add_argument('--ms', default=2, type=float, help='Marker size')
    parser.add_argument('--alpha', default=.1, type=float, help='Marker size')
    parser.add_argument('--sample1', type=int, help='Plot z value for N points')
    parser.add_argument('--sample2', type=int, help='Plot median z after chunking into N chunks')
    parser.add_argument('--out-s', help='Save sampled z values')
    parser.add_argument('--color', action='store_true')
    parser.add_argument('--seed', default=0, type=int)
    parser.add_argument('--annotate', action='store_true', help='Annotate selected points')
    parser.add_argument('--stride', action='store_true', help='Use strided points instead of random samples')
    return parser


def main(args):
    np.random.seed(args.seed)
    fig, ax = plt.subplots()
    print(args.input)
    x = pickle.load(open(args.input,'rb'))
    
    # PCA
    pca = PCA(x.shape[1])
    pca.fit(x)
    print('Explained variance ratio:')
    print(pca.explained_variance_ratio_)
    pc = pca.transform(x)
    
    if args.color:
        plt.scatter(pc[:,0], pc[:,1], c=np.arange(len(x)), label=args.input, alpha=args.alpha, s=args.ms, cmap='hsv')
    else:
        plt.scatter(pc[:,0], pc[:,1], label=args.input, alpha=args.alpha, s=args.ms)
    plt.xlabel('PC1 ({:3f})'.format(pca.explained_variance_ratio_[0]))
    plt.ylabel('PC2 ({:3f})'.format(pca.explained_variance_ratio_[1]))

    if args.sample1:
        ii = np.random.choice(len(x), args.sample1)
        if args.stride:
            ii = np.arange(len(x))[::len(x)//args.sample1]
        xd = x[ii]
        xd_pc = pca.transform(xd)
        plt.scatter(xd_pc[:,0],xd_pc[:,1],c=np.arange(len(xd)),cmap='hsv')
        if args.annotate:
            for i in range(args.sample1):
                ax.annotate(str(i), xd_pc[i,0:2])
    if args.sample2:
        xsplit = np.array_split(x,args.sample2)
        print([len(k) for k in xsplit])
        xd = np.array([np.median(xs,axis=0) for xs in xsplit])
        #xd = np.array([np.mean(xs,axis=0) for xs in xsplit])
        print(len(xd))
        xd_pc = pca.transform(xd)
        plt.scatter(xd_pc[:,0],xd_pc[:,1],c=np.arange(len(xd)),cmap='hsv')

    if args.out_s:
        np.savetxt(args.out_s, xd)
    plt.legend(loc='best')
    if args.out_png: 
        plt.savefig(args.out_png)
    else:
        plt.show()

if __name__ == '__main__':
    main(parse_args().parse_args())
