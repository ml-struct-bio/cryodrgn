'''Parse image poses from a cryoSPARC .cs metafile'''

import argparse
import numpy as np
import sys, os
import pickle
import torch

from cryodrgn import lie_tools
from cryodrgn import utils

log = utils.log

def add_args(parser):
    parser.add_argument('input', help='Cryosparc .cs file')
    parser.add_argument('--abinit', action='store_true', help='Flag if results are from ab-initio reconstruction') 
    parser.add_argument('--hetrefine', action='store_true', help='Flag if results are from a heterogeneous refinements (default: homogeneous refinement)')
    parser.add_argument('-D', type=int, required=True, help='Box size of reconstruction (pixels)')
    parser.add_argument('-o', metavar='PKL', type=os.path.abspath, required=True, help='Output pose.pkl')
    return parser

def main(args):
    assert args.input.endswith('.cs'), "Input format must be .cs file"
    assert args.o.endswith('.pkl'), "Output format must be .pkl"

    data = np.load(args.input)
    # view the first row
    for i in range(len(data.dtype)):
        print(i, data.dtype.names[i], data[0][i])

    if args.abinit:
        RKEY = 'alignments_class_0/pose'
        TKEY = 'alignments_class_0/shift'
    else:
        RKEY = 'alignments3D/pose'
        TKEY = 'alignments3D/shift'

    # parse rotations
    log(f'Extracting rotations from {RKEY}')
    rot = np.array([x[RKEY] for x in data])
    rot = torch.tensor(rot)
    rot = lie_tools.expmap(rot)
    rot = rot.numpy()
    log('Transposing rotation matrix')
    rot = np.array([x.T for x in rot])
    log(rot.shape)

    # parse translations
    log(f'Extracting translations from {TKEY}')
    trans = np.array([x[TKEY] for x in data])
    if args.hetrefine:
        log('Scaling shifts by 2x')
        trans *= 2
    log(trans.shape)
    
    # convert translations from pixels to fraction 
    trans /= args.D

    # write output
    log(f'Writing {args.o}')
    with open(args.o,'wb') as f:
        pickle.dump((rot,trans),f)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    main(add_args(parser).parse_args())
