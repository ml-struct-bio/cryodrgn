'''Align and compute distance between two series of rotation matrices'''

import argparse
import numpy as np
import sys, os
sys.path.insert(0,'{}/../lib-python'.format(os.path.dirname(os.path.abspath(__file__))))
import utils
log = utils.log
import matplotlib.pyplot as plt
from scipy.linalg import logm
import pickle


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('rot1', help='Input rotations')
    parser.add_argument('rot2', help='Input rotations')
    parser.add_argument('-i', type=int, default=0, help='Index to align on')
    parser.add_argument('-N', type=int, help='Compare first N images')
    parser.add_argument('--flip',action='store_true',help='Flip hand')
    parser.add_argument('--show', action='store_true', help='Show histogram of errors')
    parser.add_argument('-o', help='Output pickle of distances')
    return parser

def geodesic_so3(A,B):
    return np.sum(logm(np.dot(A.T,B))**2)**.5

def fast_dist(a,b):
    return np.sum((a-b)**2)

def align_rot(rotA,rotB,i,flip=False):
    if flip:
        x = -1*np.ones((3,3))
        x[2] = 1
        rotB = [x*r for r in rotB] # left multiplication of [-1,-1,1] diag matrix
                                    
    ref = rotA[i]
    rot = np.array([np.dot(x, ref.T) for x in rotA])
    ref = rotB[i]
    rot_hat = np.array([np.dot(x, ref.T) for x in rotB])
    return rot, rot_hat

def align_rot2(rotA,rotB,i,flip=False):
    if flip:
        x = -1*np.ones((3,3))
        x[2] = 1
        rotB = [x*r for r in rotB] # left multiplication of [-1,-1,1] diag matrix
    ref = np.dot(rotB[i].T, rotA[i])
    print(ref)
    rot_hat = np.array([np.dot(x, ref) for x in rotB])
    return rot_hat

def main(args):
    rot1 = utils.load_pkl(args.rot1)
    if type(rot1) == tuple:
        rot1 = rot1[0]
    #rot1 = np.array([utils.R_from_eman(*x) for x in rot1])
    rot2 = utils.load_pkl(args.rot2)
    if type(rot2) == tuple:
        rot2 = rot2[0]
    assert rot1.shape == rot2.shape

    #rot2 = align_rot2(rot1,rot2,args.i,args.flip)
    rot1, rot2 = align_rot(rot1,rot2,args.i,args.flip)
    if args.N:
        rot1 = rot1[:args.N]
        rot2 = rot2[:args.N]
    dists = np.sum((rot1-rot2)**2, axis=(1,2))
    #dists = np.asarray([geodesic_so3(a,b) for a,b in zip(rot1, rot2)])

    log('Mean error: {}'.format(np.mean(dists)))
    log('Median error: {}'.format(np.median(dists)))
    if args.show:
        plt.hist(dists)
        plt.show()

    if args.o:
        pickle.dump(dists,open(args.o,'wb'))

if __name__ == '__main__':
    main(parse_args().parse_args())
