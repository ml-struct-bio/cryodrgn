'''Plot FSC txtfile'''

import argparse
import numpy as np
import sys, os
import pickle

import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', nargs='*', help='Input')
    parser.add_argument('-t', type=float, default=.143, help='Cutoff for resolution estimation (default: %(default)s)')
    parser.add_argument('--labels', nargs='*', help='Labels for plotting')
    parser.add_argument('-o')
    return parser

def main(args):
    labels = args.labels if args.labels else args.input
    assert len(labels) == len(args.input)
    for i, f in enumerate(args.input):
        print(f)
        x = np.loadtxt(f)
        plt.plot(x[:,0],x[:,1], label=labels[i])
        w = np.where(x[:,1]<args.t)
        print(w)
        print(x[:,0][w])
        print(1/x[:,0][w])
    plt.legend(loc='best')
    plt.ylim((0,1))
    plt.ylabel('FSC')
    plt.xlabel('frequency')
    if args.o:
        plt.savefig(args.o)
    else:
        plt.show()

if __name__ == '__main__':
    main(parse_args().parse_args())
