'''Align and compute distance between two series of rotation matrices'''

import argparse
import numpy as np
import sys, os
sys.path.insert(0,'{}/../lib-python'.format(os.path.dirname(os.path.abspath(__file__))))
import utils
log = utils.log
import matplotlib.pyplot as plt
from scipy.linalg import logm


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('rot1', help='Input rotations')
    parser.add_argument('rot2', help='Input rotations')
    parser.add_argument('-N', type=int, default=30)
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
    rot2 = utils.load_pkl(args.rot2)
    assert rot1.shape == rot2.shape

    mean = []
    med = []
    for i in range(args.N):
        r1, r2 = align_rot(rot1,rot2,i,False)
        dists = np.asarray([fast_dist(a,b) for a,b in zip(r1, r2)])
        mean.append(np.mean(dists))
        med.append(np.median(dists))

    for i in range(args.N):
        r1, r2 = align_rot(rot1,rot2,i,True)
        dists = np.asarray([fast_dist(a,b) for a,b in zip(r1, r2)])
        mean.append(np.mean(dists))
        med.append(np.median(dists))

    log(np.argmin(mean))
    log('Mean error: {}'.format(np.min(mean)))
    log(np.argmin(med))
    log('Median error: {}'.format(np.min(med)))

if __name__ == '__main__':
    main(parse_args().parse_args())
