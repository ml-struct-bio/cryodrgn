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
    parser.add_argument('input', nargs='+', help='Input')
    parser.add_argument('-o', help='Output')
    parser.add_argument('--sample1', type=int, help='Plot z value for N points')
    parser.add_argument('--sample2', type=int, help='Plot median z after chunking into N chunks')
    parser.add_argument('--out-s', help='Save sampled z values')
    return parser

def main(args):
    for f in args.input:
        print(f)
        x = pickle.load(open(f,'rb'))
        assert x.shape[1] > 2
        
        # PCA
        pca = PCA(2)
        pca.fit(x)
        print('Explained variance ratio:')
        print(pca.explained_variance_ratio_)
        pc = pca.transform(x)

        plt.scatter(pc[:,0], pc[:,1], c=np.arange(len(x)), label=f, alpha=.1, s=2, cmap='hsv')

    if args.sample1:
        d = len(x) // args.sample1
        xd = x[::d]
        print(len(xd))
        xd_pc = pca.transform(xd)
        plt.scatter(xd_pc[:,0],xd_pc[:,1],c=np.arange(len(xd)),cmap='hsv')
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
    if args.o: 
        plt.savefig(args.o)
    else:
        plt.show()

if __name__ == '__main__':
    main(parse_args().parse_args())
