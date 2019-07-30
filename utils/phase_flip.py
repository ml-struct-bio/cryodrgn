'''Phase flip an image stack by its CTF'''

import argparse
import numpy as np
import sys, os

sys.path.insert(0, '/home/zhonge/dev/cryovae/master/lib-python')

import utils
import mrc
import ctf
import fft

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('particles', help='Input')
    parser.add_argument('ctf_params', help='Input')
    parser.add_argument('-o', help='Output')
    return parser

def main(args):
    imgs, _, _ = mrc.parse_mrc(args.particles,lazy=True)
    ctf_params = utils.load_pkl(args.ctf_params)
    assert len(imgs) == len(ctf_params)
    
    D = imgs[0].get().shape[0]
    fx, fy= np.meshgrid(np.linspace(-.5,.5,D,endpoint=False),np.linspace(-.5,.5,D,endpoint=False))
    freqs = np.stack([fx.ravel(), fy.ravel()],1)
    
    imgs_flip = np.empty((len(imgs),D,D),dtype=np.float32)
    for i in range(len(imgs)):
        if i%1000 == 0: print(i)
        c = ctf.compute_ctf_np(freqs/ctf_params[i,0], *ctf_params[i,1:])
        c = c.reshape((D,D))
        ff = fft.fft2_center(imgs[i].get())
        ff *= np.sign(c)
        img = fft.ifftn_center(ff)
        imgs_flip[i] = img.astype(np.float32)

    mrc.write(args.o, imgs_flip)

    

if __name__ == '__main__':
    main(parse_args().parse_args())
