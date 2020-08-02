'''Filter a particle stack'''

import argparse
import numpy as np
import sys, os
import pickle

from cryodrgn import utils
from cryodrgn import dataset
from cryodrgn import mrc
log = utils.log 

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='Input particles (.mrcs, .txt, .star, .cs)')
    parser.add_argument('--ind', required=True, help='Selected indices array (.pkl)')
    parser.add_argument('-o', help='Output .mrcs file')
    return parser

def main(args):
    x = dataset.load_particles(args.input, lazy=True)
    log(f'Loaded {len(x)} particles')
    ind = utils.load_pkl(args.ind)
    x = np.array([x[i].get() for i in ind])
    log(f'New stack dimensions: {x.shape}')
    mrc.write(args.o,x)

if __name__ == '__main__':
    main(parse_args().parse_args())
