'''Phase flip an image stack by its CTF'''

import argparse
import numpy as np
import sys, os

from cryodrgn import utils
from cryodrgn import mrc
from cryodrgn import dataset
from cryodrgn import ctf
from cryodrgn import fft

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mrcs', help='Input particles (.mrcs, .txt, .star, or .cs)')
    parser.add_argument('ctf_params', help='Input CTF parameters (.pkl)')
    parser.add_argument('--datadir', help='Optionally overwrite path to starfile .mrcs if loading from a starfile')
    parser.add_argument('-o', help='Output')
    return parser

def main(args):
    imgs = dataset.load_particles(args.mrcs, lazy=True, datadir=args.datadir)
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
