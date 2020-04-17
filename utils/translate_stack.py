'''Translate a particle stack'''

import argparse
import numpy as np
import sys, os
import matplotlib.pyplot as plt

from cryodrgn import utils
from cryodrgn import dataset
from cryodrgn import mrc
from cryodrgn import fft

log = utils.log

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mrcs', help='Input particles (.mrcs, .cs, .star, or .txt)')
    parser.add_argument('trans', help='Pose or translations pickle (.pkl)')
    parser.add_argument('--tscale', type=float, default=1.0, help='Scale translations by this amount (default: %(default)s)')
    parser.add_argument('--datadir', help='Optionally overwrite path to starfile .mrcs if loading from a starfile')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output particle stack')
    parser.add_argument('--out-png')
    return parser

def plot_projections(out_png, imgs):
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(10,10))
    axes = axes.ravel()
    for i in range(min(len(imgs),9)):
        axes[i].imshow(imgs[i])
    plt.savefig(out_png)

def main(args):
    # load particles
    particles = dataset.load_particles(args.mrcs, datadir=args.datadir)
    log(particles.shape)
    Nimg, D, D = particles.shape

    trans = utils.load_pkl(args.trans)
    if type(trans) is tuple:
        trans = trans[1]
    trans *= args.tscale
    assert np.all(trans <= 1), "ERROR: Old pose format detected. Translations must be in units of fraction of box."
    trans *= D # convert to pixels
    assert len(trans) == Nimg

    xx,yy = np.meshgrid(np.arange(-D/2,D/2),np.arange(-D/2,D/2))
    TCOORD = np.stack([xx, yy],axis=2)/D # DxDx2
    
    imgs = []
    for ii in range(Nimg):
        ff = fft.fft2_center(particles[ii])
        tfilt = np.dot(TCOORD,trans[ii])*-2*np.pi
        tfilt = np.cos(tfilt) + np.sin(tfilt)*1j
        ff *= tfilt
        img = fft.ifftn_center(ff)
        imgs.append(img)

    imgs = np.asarray(imgs).astype(np.float32)
    mrc.write(args.o, imgs)

    if args.out_png:
        plot_projections(args.out_png, imgs[:9])

if __name__ == '__main__':
    main(parse_args().parse_args())
