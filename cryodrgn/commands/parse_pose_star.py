'''Parse image poses from RELION .star file'''

import argparse
import numpy as np
import sys, os
import pickle

from cryodrgn import utils
from cryodrgn import starfile
log = utils.log

def add_args(parser):
    parser.add_argument('input', help='RELION .star file')
    parser.add_argument('-D', type=int, required=True, help='Box size of reconstruction (pixels)')
    parser.add_argument('-o', metavar='PKL', type=os.path.abspath, required=True, help='Output pose.pkl')
    return parser

def main(args):
    assert args.input.endswith('.star'), "Input file must be .star file"
    assert args.o.endswith('.pkl'), "Output format must be .pkl"

    s = starfile.Starfile.load(args.input)
    N = len(s.df)
    log('{} particles'.format(N))
    
    # parse rotations
    keys = ('_rlnAngleRot','_rlnAngleTilt','_rlnAnglePsi')
    euler = np.empty((N,3))
    euler[:,0] = s.df['_rlnAngleRot']
    euler[:,1] = s.df['_rlnAngleTilt']
    euler[:,2] = s.df['_rlnAnglePsi']
    log('Euler angles (Rot, Tilt, Psi):')
    log(euler[0])
    log('Converting to rotation matrix:')
    rot = np.asarray([utils.R_from_relion(*x) for x in euler])
    log(rot[0])

    # parse translations
    trans = np.empty((N,2))
    trans[:,0] = s.df['_rlnOriginX']
    trans[:,1] = s.df['_rlnOriginY']
    log('Translations:')
    log(trans[0])
    
    # convert translations from pixels to fraction 
    trans /= args.D

    # write output
    log(f'Writing {args.o}')
    with open(args.o,'wb') as f:
        pickle.dump((rot,trans),f)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    main(add_args(parser).parse_args())
