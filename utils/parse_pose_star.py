'''Parse image poses from RELION .star file'''

import argparse
import numpy as np
import sys, os
import pickle
sys.path.insert(0,'{}/../lib-python'.format(os.path.dirname(os.path.abspath(__file__))))
import utils
import starfile
log = utils.log

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='RELION .star file')
    parser.add_argument('-D', type=int, required=True, help='Box size of reconstruction (pixels)')
    parser.add_argument('-o', required=True, help='Output prefix for appending .rot.pkl and .trans.pkl')
    return parser

def main(args):
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
