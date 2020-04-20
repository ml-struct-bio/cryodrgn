'''
Plot PCA projection of latent space
'''

import argparse
import numpy as np
import sys, os
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',  help='Input z pkl')
    parser.add_argument('-o', '--out-png', help='Output PNG')
    parser.add_argument('--axis', type=int, nargs=2, default=[0,1], help='Dimensions to plot (default: %(default)s)')
    parser.add_argument('--ms', default=2, type=float, help='Marker size for plotting (default: %(default)s)')
    parser.add_argument('--alpha', default=.1, type=float, help='Alpha value for plotting (default: %(default)s)')
    parser.add_argument('--sample1', type=int, help='Optionally plot z value for N randomly sampled points')
    parser.add_argument('--sample2', type=int, help='Optionally plot median z after chunking into N chunks')
    parser.add_argument('--out-s', help='Save sampled z values (.txt)')
    parser.add_argument('--color', action='store_true', help='Color points by image index')
    parser.add_argument('--seed', default=0, type=int, help='Random seed (default: %(default)s)')
    parser.add_argument('--annotate', action='store_true', help='Annotate sampled points in plot')
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
    
    ii,jj = args.axis
    if args.color:
        plt.scatter(pc[:,ii], pc[:,jj], c=np.arange(len(x)), label=args.input, alpha=args.alpha, s=args.ms, cmap='hsv')
    else:
        plt.scatter(pc[:,ii], pc[:,jj], label=args.input, alpha=args.alpha, s=args.ms)
    plt.xlabel('PC{} ({:3f})'.format(ii+1, pca.explained_variance_ratio_[ii]))
    plt.ylabel('PC{} ({:3f})'.format(jj+1, pca.explained_variance_ratio_[jj]))

    if args.sample1:
        s = np.random.choice(len(x), args.sample1)
        print(s)
        xd = x[s]
        xd_pc = pca.transform(xd)
        plt.scatter(xd_pc[:,ii],xd_pc[:,jj],c=np.arange(len(xd)),cmap='hsv')
        if args.annotate:
            for i in range(args.sample1):
                ax.annotate(str(i), xd_pc[i,args.axis])
    if args.sample2:
        xsplit = np.array_split(x,args.sample2)
        print([len(k) for k in xsplit])
        xd = np.array([np.median(xs,axis=0) for xs in xsplit])
        #xd = np.array([np.mean(xs,axis=0) for xs in xsplit])
        print(len(xd))
        xd_pc = pca.transform(xd)
        plt.scatter(xd_pc[:,ii],xd_pc[:,jj],c='k')#np.arange(len(xd)),cmap='hsv')

    if args.out_s:
        np.savetxt(args.out_s, xd)
    plt.legend(loc='best')
    if args.out_png: 
        plt.savefig(args.out_png)
    else:
        plt.show()

if __name__ == '__main__':
    main(parse_args().parse_args())
