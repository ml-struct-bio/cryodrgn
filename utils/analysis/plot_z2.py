'''
Plot 2D latent space
'''

import argparse
import numpy as np
import sys, os
import pickle
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',  help='Input z pkl')
    parser.add_argument('-o', '--out-png', help='Output PNG')
    parser.add_argument('--ms', default=2, type=float, help='Marker size for plotting (default: %(default)s)')
    parser.add_argument('--alpha', default=.1, type=float, help='Alpha value for plotting (default: %(default)s)')
    parser.add_argument('--sample1', type=int, help='Optionally plot z value for N randomly sampled points')
    parser.add_argument('--sample2', type=int, help='Optionally plot median z after chunking into N chunks')
    parser.add_argument('--out-s', help='Save sampled z values (.txt)')
    parser.add_argument('--color', action='store_true', help='Color points by image index')
    parser.add_argument('--seed', default=0, type=int, help='Random seed (default: %(default)s)')
    parser.add_argument('--annotate', action='store_true', help='Annotate sampled points in plot')
    parser.add_argument('--kde', action='store_true', help='KDE plot instead of scatter')
    parser.add_argument('--stride', type=int, help='Stride dataset')
    return parser

def main(args):
    np.random.seed(args.seed)
    f = args.input
    print(f)
    x = pickle.load(open(f,'rb'))
    if args.stride:
        x = x[::args.stride]
    print(x.shape)

    # seaborn jointpoint
    if args.kde:
        g = sns.jointplot(x[:,0],x[:,1], kind='kde')
        ax = g.ax_joint 

    # scatter plot
    else:
        fig,ax = plt.subplots()
        if args.color:
            plt.scatter(x[:,0], x[:,1], c=np.arange(len(x[:,0])), label=f, alpha=args.alpha, s=args.ms, cmap='hsv')
        else:
            plt.scatter(x[:,0], x[:,1], label=f, alpha=args.alpha, s=args.ms)
        plt.xlabel('z1')
        plt.ylabel('z2')
        plt.legend(loc='best')

    if args.sample1:
        ii = np.random.choice(len(x), args.sample1)
        print(ii)
        xd = x[ii]
        print(xd)
        plt.scatter(xd[:,0],xd[:,1],c=np.arange(len(xd)),cmap='hsv')
        if args.annotate:
            for i in range(args.sample1):
                ax.annotate(str(i), xd[i])
    if args.sample2:
        xsplit = np.array_split(x,args.sample2)
        xd = np.array([np.median(xs,axis=0) for xs in xsplit])
        print(len(xd))
        print(xd)
        plt.scatter(xd[:,0],xd[:,1],c='k')#np.arange(len(xd)),cmap='hsv')
    if args.out_s:
        np.savetxt(args.out_s, xd)
    if args.out_png: 
        plt.savefig(args.out_png)
    else:
        plt.show()

if __name__ == '__main__':
    main(parse_args().parse_args())
