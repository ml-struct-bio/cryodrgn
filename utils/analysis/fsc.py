'''Compute FSC between two volumes'''

import argparse
import numpy as np
import sys, os
import matplotlib.pyplot as plt

from cryodrgn import utils
from cryodrgn import mrc
from cryodrgn import ctf
from cryodrgn import fft

log = utils.log

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vol1', help='Input')
    parser.add_argument('vol2', help='Input')
    parser.add_argument('--mask')
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--Apix', type=float, default=1)
    parser.add_argument('-o', help='Output')
    return parser

def main(args):
    vol1 , _ = mrc.parse_mrc(args.vol1)
    vol2 , _ = mrc.parse_mrc(args.vol2)

    if args.mask:
        mask = mrc.parse_mrc(args.mask)[0]
        vol1 *= mask
        vol2 *= mask
    
    D = vol1.shape[0]
    x = np.arange(-D//2, D//2)
    x2, x1, x0 = np.meshgrid(x,x,x, indexing='ij')
    coords = np.stack((x0,x1,x2), -1)
    r = (coords**2).sum(-1)**.5

    assert r[D//2, D//2, D//2] == 0.0

    vol1 = fft.fftn_center(vol1)
    vol2 = fft.fftn_center(vol2)

    #log(r[D//2, D//2, D//2:])
    prev_mask = np.zeros((D,D,D), dtype=bool)
    fsc = [1.0]
    for i in range(1,D//2):
        mask = r < i
        shell = np.where(mask & np.logical_not(prev_mask))
        v1 = vol1[shell]
        v2 = vol2[shell]
        p = np.vdot(v1,v2) / (np.vdot(v1,v1)*np.vdot(v2,v2))**.5
        fsc.append(p.real)
        prev_mask = mask
    fsc = np.asarray(fsc)
    x = np.arange(D//2)/D

    res = np.stack((x,fsc),1)
    if args.o:
        np.savetxt(args.o, res)
    else:
        log(res)

    w = np.where(fsc < 0.5)
    if w:
        log('0.5: {}'.format(1/x[w[0]]*args.Apix))

    w = np.where(fsc < 0.143)
    if w:
        log('0.143: {}'.format(1/x[w[0]]*args.Apix))
    
    if args.plot:
        plt.plot(x,fsc)
        plt.ylim((0,1))
        plt.show()


    
   

if __name__ == '__main__':
    main(parse_args().parse_args())
