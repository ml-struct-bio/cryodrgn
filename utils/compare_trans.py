'''Compute error between two series of particle shifts'''

import argparse
import numpy as np
import sys, os
sys.path.insert(0,'{}/../lib-python'.format(os.path.dirname(os.path.abspath(__file__))))
import utils
log = utils.log
vlog = utils.vlog
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trans1', help='Input translations')
    parser.add_argument('trans2', help='Input translations')
    parser.add_argument('--s1',type=float,default=1.0,help='Scale for trans1')
    parser.add_argument('--s2',type=float,default=1.0,help='Scale for trans2')
    parser.add_argument('--show', action='store_true', help='Show histogram')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbosity')
    return parser

def main(args):
    trans1 = utils.load_pkl(args.trans1)
    if type(trans1) == tuple:
        trans1 = trans1[1]
    trans1 *= args.s1
    vlog(trans1.shape)
    vlog(trans1)
    trans2 = utils.load_pkl(args.trans2)
    if type(trans2) == tuple:
        trans2 = trans2[1]
    trans2 *= args.s2
    vlog(trans2.shape)
    vlog(trans2)
    assert trans1.shape == trans2.shape
    vlog(np.mean(trans1,axis=0))
    vlog(np.mean(trans2,axis=0))

    dists = np.sum((trans1-trans2)**2,axis=1)**.5
    vlog(dists.shape)

    log('Mean error: {}'.format(np.mean(dists)))
    log('Median error: {}'.format(np.median(dists)))
    if args.show:
        plt.figure(1)
        plt.hist(dists)
        plt.figure(2)
        plt.scatter(trans1[:,0],trans1[:,1], s=1, alpha=.1)
        plt.figure(3)
        plt.scatter(trans2[:,0],trans2[:,1], s=1, alpha=.1)
        plt.figure(4)
        d = trans1 - trans2
        plt.scatter(d[:,0],d[:,1], s=1, alpha=.1)
        plt.show()

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)
