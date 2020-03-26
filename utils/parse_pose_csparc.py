'''Parse pose from a cryoSPARC .cs metafile'''

import argparse
import numpy as np
import sys, os
import pickle
import torch
sys.path.insert(0,'{}/../lib-python'.format(os.path.dirname(os.path.abspath(__file__))))
import utils
import lie_tools

log = utils.log

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='Cryosparc .cs file')
    parser.add_argument('--abinit', action='store_true', help='Flag if results are from ab-initio reconstruction') 
    parser.add_argument('--hetrefine', action='store_true', help='Flag if results are from a heterogeneous refinements (default: homogeneous refinement)')
    parser.add_argument('-D', type=int, required=True, help='Box size of reconstruction (pixels)')
    parser.add_argument('-o', required=True, help='Output prefix for appending .rot.pkl and .trans.pkl')
    return parser

def main(args):
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
    out_rot = '{}.rot.pkl'.format(args.o)
    log('Writing {}'.format(out_rot))
    with open(out_rot,'wb') as f:
        pickle.dump(rot,f)
    out_trans = '{}.trans.pkl'.format(args.o)
    log('Writing {}'.format(out_trans))
    with open(out_trans,'wb') as f:
        pickle.dump(trans,f)

if __name__ == '__main__':
    main(parse_args().parse_args())
