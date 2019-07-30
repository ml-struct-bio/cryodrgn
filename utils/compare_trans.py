'''Compute error between two series of particle shifts'''

import argparse
import numpy as np
import sys, os
sys.path.insert(0,'{}/../lib-python'.format(os.path.dirname(os.path.abspath(__file__))))
import utils
log = utils.log
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trans1', help='Input translations')
    parser.add_argument('trans2', help='Input translations')
    parser.add_argument('--s1',type=float,default=1.0,help='Scale for trans1')
    parser.add_argument('--s2',type=float,default=1.0,help='Scale for trans2')
    parser.add_argument('--show', action='store_true', help='Show histogram')
    return parser

def main(args):
    trans1 = utils.load_pkl(args.trans1)*args.s1
    log(trans1.shape)
    log(trans1)
    trans2 = utils.load_pkl(args.trans2)*args.s2
    log(trans2.shape)
    log(trans2)
    assert trans1.shape == trans2.shape
    log(np.mean(trans1,axis=0))
    log(np.mean(trans2,axis=0))

    dists = np.sum((trans1-trans2)**2,axis=1)**.5
    log(dists.shape)

    log('Mean error: {}'.format(np.mean(dists)))
    log('Median error: {}'.format(np.median(dists)))
    if args.show:
        plt.hist(dists)
        plt.show()

if __name__ == '__main__':
    main(parse_args().parse_args())
