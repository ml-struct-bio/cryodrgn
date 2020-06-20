'''
Filter a .star file

(Note: RELION 3.1 star files are not supported)
'''

import argparse
import numpy as np
import sys, os
import pickle

from cryodrgn import utils
from cryodrgn import starfile
log = utils.log 

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='Input .star file')
    parser.add_argument('--ind', required=True, help='Selected indices array (.pkl)')
    parser.add_argument('-o', help='Output .star file')
    return parser

def main(args):
    s = starfile.Starfile.load(args.input, relion31=False)
    ind = utils.load_pkl(args.ind)
    log('{} particles'.format(len(s.df)))
    s.df = s.df.loc[ind]
    log(len(s.df))
    log(len(s.df.index))
    s.write(args.o)

if __name__ == '__main__':
    main(parse_args().parse_args())
