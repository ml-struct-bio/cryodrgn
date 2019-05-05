'''Translate a particle stack'''

import argparse
import numpy as np
import sys, os
import matplotlib.pyplot as plt

sys.path.insert(0,'{}/../lib-python'.format(os.path.dirname(os.path.abspath(__file__))))
import utils
import mrc
import fft

log = utils.log

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=os.path.abspath, help='Input particle stack')
    parser.add_argument('trans', help='Pickle with image shifts')
    parser.add_argument('--tscale', type=float, 'Scale translations by this amount')
    parser.add_argument('-o', type=os.path.abspath, required=True, 
        help='Output particle stack')
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
    particles, _, _ = mrc.parse_mrc(args.input)
    log(particles.shape)
    Nimg, D, D = particles.shape

    trans = utils.load_pkl(args.trans)
    if args.tscale:
        trans *= args.tscale
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
