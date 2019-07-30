'''Parse out 3D alignment poses from cryosparc .cs metafiles'''

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
    parser.add_argument('-o', help='Output prefix for appending .rot.pkl and .trans.pkl')
    return parser

def main(args):
    data = np.load(args.input)
    # view the first row
    for i in range(len(data.dtype)):
        print(i, data.dtype.names[i], data[0][i])
    rot = np.array([x['alignments_class_0/pose'] for x in data])
    rot = torch.tensor(rot)
    rot = lie_tools.expmap(rot)
    rot = rot.numpy()
    log('Transposing rotation matrix')
    rot = np.array([x.T for x in rot])
    log(rot.shape)
    trans = np.array([x['alignments_class_0/shift'] for x in data])
    log('Scaling shifts by 2x')
    trans *= 2
    log(trans.shape)
    
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
