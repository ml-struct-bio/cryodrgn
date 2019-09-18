'''
'''

import argparse
import numpy as np
import sys, os
import pickle
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', nargs='+', help='Input')
    parser.add_argument('-o', help='Output')
    return parser

def main(args):
    for f in args.input:
        print(f)
        x = pickle.load(open(f,'rb'))
        #plt.plot(x[:,0],x[:,1], 'o', label=f, alpha=.01, ms=2)
        plt.scatter(x[:,0], x[:,1], c=np.arange(len(x[:,0])), label=f, alpha=.02, s=4, cmap='hsv')
    plt.xlabel('z1')
    plt.ylabel('z2')
    plt.legend(loc='best')
    if args.o: 
        plt.savefig(args.o)
    else:
        plt.show()

if __name__ == '__main__':
    main(parse_args().parse_args())
