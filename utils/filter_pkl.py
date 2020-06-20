'''Filter cryodrgn data stored in a .pkl file'''

import argparse
import numpy as np
import sys, os
import pickle

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='Input data .pkl')
    parser.add_argument('--ind', required=True, help='Selected indices array (.pkl)')
    parser.add_argument('-o', help='Output data .pkl')
    return parser

def load_pkl(x):
    return pickle.load(open(x,'rb'))

def main(args):
    x = load_pkl(args.input)
    ind = load_pkl(args.ind)

    # pose.pkl contains a tuple of rotations and translations
    if type(x) == tuple:
        print([xx.shape for xx in x])
        x = (xx[ind] for xx in x)
        x = tuple(x)
        print([xx.shape for xx in x])

    # all other cryodrgn pkls 
    else:
        print(x.shape)
        x = x[ind]
        print(x.shape)
    print(f'Saving {args.o}')
    pickle.dump(x, open(args.o,'wb'))

if __name__ == '__main__':
    main(parse_args().parse_args())
