'''
'''

import argparse
import numpy as np
import sys, os
import pickle
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',  help='Input z.pkl')
    parser.add_argument('-o', help='Output PNG')
    parser.add_argument('--ms', default=2, type=float, help='Marker size for plotting (default: %(default)s)')
    parser.add_argument('--alpha', default=.1, type=float, help='Alpha value for plotting (default: %(default)s)')
    parser.add_argument('--ylim', nargs=2, type=float)
    parser.add_argument('--sample1', type=int, help='Plot z value for N randomly sampled points')
    parser.add_argument('--sample2', type=int, help='Plot median z after chunking into N chunks')
    parser.add_argument('--seed', default=0, type=int, help='Random seed (default: %(default)s)')
    parser.add_argument('--out-s', help='Save sampled z values (.txt)')
    return parser

def main(args):
    np.random.seed(args.seed)
    f = args.input
    print(f)
    fi = open(f,'rb')
    x = pickle.load(fi)
    N = len(x)
    plt.scatter(np.arange(N),x, label=f, alpha=args.alpha, s=args.ms)
    #plt.scatter(np.arange(N), x, c=np.arange(len(x[:,0])), label=f, alpha=.1, s=2, cmap='hsv')
    plt.xlim((0,N))
    if args.sample1:
        s = np.random.choice(len(x), args.sample1)
        xd = x[s]
        print(xd)
        plt.plot(s,xd,'o')

    if args.sample2:
        t = np.array_split(np.arange(len(x)),args.sample2)
        t = np.array([np.median(tt,axis=0) for tt in t])
        xsplit = np.array_split(x,args.sample2)
        xd = np.array([np.median(xs,axis=0) for xs in xsplit])
        print(len(xd))
        print(xd)
        plt.plot(t,xd,'o',color='k')
    if args.out_s:
        np.savetxt(args.out_s, xd)
    if args.ylim:
        plt.ylim(args.ylim)
    plt.xlabel('image')
    plt.ylabel('latent encoding')
    plt.legend(loc='best')
    if args.o: 
        plt.savefig(args.o)

    # Plot histogram
    plt.figure()
    sns.distplot(x)
    plt.show()

if __name__ == '__main__':
    main(parse_args().parse_args())
