'''
Filter a .star file
'''

import argparse
import numpy as np
import sys, os
import pickle

from cryodrgn import utils
from cryodrgn import starfile
log = utils.log 

def add_args(parser):
    parser.add_argument('input', help='Input .star file')
    parser.add_argument('--ind', required=True, help='Array of selected indices (.pkl)')
    parser.add_argument('-o', type=os.path.abspath, help='Output .star file')
    return parser

def main(args):
    s = starfile.Starfile.load(args.input)
    ind = utils.load_pkl(args.ind)
    log('Loaded {} particles'.format(len(s.df)))
    log(f'Index array: {ind}')
    s.df = s.df.loc[ind]
    log('Filtered to {} particles'.format(len(s.df)))
    s.write(args.o)
    log(f'Saved {args.o}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
