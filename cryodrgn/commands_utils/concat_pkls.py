'''Concatenate arrays from multiple .pkl files'''

import argparse
import numpy as np
import sys, os
import pickle

log = print

def add_args(parser):
    parser.add_argument('input', nargs='+', help='Input .pkl files')
    parser.add_argument('-o', required=True, help='Output .pkl file')
    return parser

def main(args):
    x = [pickle.load(open(f,'rb')) for f in args.input]
    if type(x[0]) == tuple: # pose tuples
        r = [xx[0] for xx in x]
        t = [xx[1] for xx in x]
        r2 = np.concatenate(r)
        t2 = np.concatenate(t)
        log(r2.shape)
        log(t2.shape)
        x2 = (r2,t2)
    else:
        for i in x:
            log(i.shape)
        x2 = np.concatenate(x)
        log(x2.shape)
    pickle.dump(x2, open(args.o,'wb'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
