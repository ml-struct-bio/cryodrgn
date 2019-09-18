'''
'''

import argparse
import numpy as np
import sys, os
import pickle
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', nargs='+', help='Input')
    parser.add_argument('-o', help='Output')
    return parser

def pca(xin):
    xout = PCA(2).fit_transform(xin)
    return xout

def main(args):
    for f in args.input:
        print(f)
        x = pickle.load(open(f,'rb'))
        assert x.shape[1] > 2
        pc = pca(x)
        #plt.plot(pc[:,0], pc[:,1], 'o', label=f, alpha=.02, ms=2)
        plt.scatter(pc[:,0], pc[:,1], c=np.arange(len(x)), label=f, alpha=.02, s=2, cmap='hsv')
    plt.legend(loc='best')
    if args.o: 
        plt.savefig(args.o)
    else:
        plt.show()

if __name__ == '__main__':
    main(parse_args().parse_args())
