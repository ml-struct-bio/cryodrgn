'''Save BNB rotations and translation from weights.pkl as separate files'''

import argparse
import numpy as np
import sys, os

import pickle
import torch

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('weights', help='Input')
    parser.add_argument('-o', help='Output prefix, to append .rot.pkl, .trans.pkl')
    return parser

def main(args):
    d = torch.load(args.weights)
    pose = d['bnb_pose']
    ind = [x[0] for x in pose]
    ind = np.concatenate(ind)
    rot = [x[1][0] for x in pose]
    rot = np.concatenate(rot)
    rot = rot[np.argsort(ind)]
    print(rot.shape)
    if len(pose[0][1]) == 2:
        trans = [x[1][1] for x in pose]
        trans = np.concatenate(trans)
        trans = trans[np.argsort(ind)]
        print(trans.shape)
    else:
        trans = None
    
    print('{}.rot.pkl'.format(args.o))
    with open('{}.rot.pkl'.format(args.o),'wb') as f:
        pickle.dump(rot,f)
    if trans is not None:
        print('{}.trans.pkl'.format(args.o))
        with open('{}.trans.pkl'.format(args.o),'wb') as f:
            pickle.dump(trans,f)

if __name__ == '__main__':
    main(parse_args().parse_args())
