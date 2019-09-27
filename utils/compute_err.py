'''Compare two (aligned) volumes'''

import argparse
import numpy as np
import sys, os
import pickle

sys.path.insert(0, '{}/../lib-python'.format(os.path.dirname(os.path.realpath(__file__))))
import mrc
import utils
log = utils.log 

import matplotlib.pyplot as plt
from skimage.measure import compare_ssim

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('ref', help='Input')
    parser.add_argument('align', nargs='+', help='Input')
    parser.add_argument('--plot', action='store_true')
    return parser

def norm(x):
    return (x-x.mean())/x.std()

def NMSE(ref,a):
    return np.mean((ref-a)**2)/np.mean(ref**2)

def CC(ref,a):
    return np.mean(ref*a)/ref.std()/a.std()

def NCC(ref,a):
    return np.mean((ref-ref.mean())*(a-a.mean()))/ref.std()/a.std()

def main(args):
    ref,_,_ = mrc.parse_mrc(args.ref)
    log('Min   : {}'.format(np.min(ref)))
    log('Max   : {}'.format(np.max(ref)))
    log('Median: {}'.format(np.median(ref)))
    mask = np.where(ref)
    for x in args.align:
        log(x)
        a,_,_ = mrc.parse_mrc(x)
        log('Min   : {}'.format(np.min(a)))
        log('Max   : {}'.format(np.max(a)))
        log('Median: {}'.format(np.median(a)))
        assert ref.shape == a.shape
        log('MSE: {}'.format(np.mean((ref-a)**2)))
        log('NMSE: {}'.format(NMSE(ref,a)))
        #log('Normalize+MSE: {}'.format(np.mean((norm(ref)-norm(a))**2)))
        #log('Max Normalize+MSE: {}'.format(np.mean(ref/ref.max() - a/a.max())**2))
        #log('Masked MSE: {}'.format(np.mean((norm(ref[mask])-norm(a[mask]))**2)))
        log('SSIM: {}'.format(compare_ssim(ref,a)))
        #log('CC: {}'.format(CC(ref,a)))
        log('NCC: {}'.format(NCC(ref,a)))
        print('')
        if args.plot:
            plt.hist(a.ravel(), 20, alpha=0.1, label=x)
    if args.plot:
        plt.hist(ref.ravel(), 20, alpha=0.1, label='ref', range=(a.min(),a.max()))
        plt.yscale('log')
        plt.legend(loc='best')
        plt.show()

if __name__ == '__main__':
    main(parse_args().parse_args())
