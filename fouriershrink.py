'''
Resize an image stack by clipping fourier frequencies
'''

import argparse
import numpy as np
import sys, os

sys.path.insert(0,'{}/lib-python'.format(os.path.dirname(os.path.abspath(__file__))))
import utils
import mrc
import fft

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

log = utils.log

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mrcs', help='Input projection stack (.mrcs)')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output projection stack (.mrcs)')
    parser.add_argument('--out-png', type=os.path.abspath, help='Montage of first 9 projections')
    parser.add_argument('-D', type=int, required=True, help='New image size, must be even')
    return parser

def mkbasedir(out):
    if not os.path.exists(os.path.dirname(out)):
        os.makedirs(os.path.dirname(out))

def warnexists(out):
    if os.path.exists(out):
        log('Warning: {} already exists. Overwriting.'.format(out))

def plot_projections(out_png, imgs):
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(10,10))
    axes = axes.ravel()
    for i in range(min(len(imgs),9)):
        axes[i].imshow(imgs[i])
    plt.savefig(out_png)


def main(args):
    mkbasedir(args.o)
    warnexists(args.o)

    old, _, _ = mrc.parse_mrc(args.mrcs)
    
    oldD = old.shape[1]
    assert args.D < oldD
    assert args.D % 2 == 0
    
    D = args.D

    start = int(oldD/2 - D/2)
    stop = int(oldD/2 + D/2)

    oldft = np.array([fft.ht2_center(img) for img in old])
    log(oldft.shape)
    newft = oldft[:,start:stop, start:stop]
    log(newft.shape)

    assert oldft[0,int(oldD/2),int(oldD/2)] == newft[0,int(D/2),int(D/2)]

    new = np.array([fft.ihtn_center(img) for img in newft])

    log('Saving {}'.format(args.o))
    mrc.write(args.o,new.astype(np.float32))

    if args.out_png:
        log('Saving {}'.format(args.out_png))
        plot_projections(args.out_png, new[:9])

if __name__ == '__main__':
    main(parse_args().parse_args())
