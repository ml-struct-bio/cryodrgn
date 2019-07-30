'''Skeleton script'''

import argparse
import numpy as np
import sys, os
import pickle

sys.path.insert(0,'/home/zhonge/dev/cryovae/master/lib-python')
import utils
log = utils.log 

import mrc
import fft


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='Input')
    parser.add_argument('-o', help='Output')
    return parser

def main(args):
    vol, _, _ = mrc.parse_mrc(args.input)
    volf = fft.htn_center(vol)
    D = vol.shape[0]
    xx = np.linspace(-1,1,D,endpoint=False)
    z,y,x = np.meshgrid(xx,xx,xx)
    coords = np.stack((x,y,z),-1)
    r = np.sum(coords**2,axis=-1)**.5
    log('Zeroing {} pixels'.format(len(np.where(r>1)[0])))
    volf[np.where(r>1)] = 0
    vol = fft.ihtn_center(volf)
    mrc.write(args.o,vol)

if __name__ == '__main__':
    main(parse_args().parse_args())
