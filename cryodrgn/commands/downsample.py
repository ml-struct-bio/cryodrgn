'''
Downsample an image stack or volume by clipping fourier frequencies
'''

import argparse
import numpy as np
import sys, os
import math
import multiprocessing as mp
from multiprocessing import Pool

from cryodrgn import utils
from cryodrgn import mrc
from cryodrgn import fft
from cryodrgn import dataset

log = utils.log

def add_args(parser):
    parser.add_argument('mrcs', help='Input particles or volume (.mrc, .mrcs, .star, or .txt)')
    parser.add_argument('-D', type=int, required=True, help='New box size in pixels, must be even')
    parser.add_argument('-o', metavar='MRCS', type=os.path.abspath, required=True, help='Output projection stack (.mrcs)')
    parser.add_argument('-b', type=int, default=5000, help='Batch size for processing images (default: %(default)s)')
    parser.add_argument('--is-vol',action='store_true', help='Flag if input .mrc is a volume')
    parser.add_argument('--chunk', type=int, help='Chunksize (in # of images) to split particle stack when saving')
    parser.add_argument('--relion31', action='store_true', help='Flag for relion3.1 star format')
    parser.add_argument('--datadir', help='Optionally provide path to input .mrcs if loading from a .star or .cs file')
    parser.add_argument('--max-threads', type=int, default=16, help='Maximum number of CPU cores for parallelization (default: %(default)s)')
    return parser

def mkbasedir(out):
    if not os.path.exists(os.path.dirname(out)):
        os.makedirs(os.path.dirname(out))

def warnexists(out):
    if os.path.exists(out):
        log('Warning: {} already exists. Overwriting.'.format(out))

def main(args):
    mkbasedir(args.o)
    warnexists(args.o)
    assert (args.o.endswith('.mrcs') or args.o.endswith('mrc')), "Must specify output in .mrc(s) file format"

    lazy = not args.is_vol
    old = dataset.load_particles(args.mrcs, lazy=lazy, datadir=args.datadir, relion31=args.relion31)

    oldD = old[0].get().shape[0] if lazy else old.shape[-1]
    assert args.D <= oldD, f'New box size {args.D} cannot be larger than the original box size {oldD}'
    assert args.D % 2 == 0, 'New box size must be even'
    
    D = args.D
    start = int(oldD/2 - D/2)
    stop = int(oldD/2 + D/2)

    def _combine_imgs(imgs):
        ret = []
        for img in imgs:
            img.shape = (1,*img.shape) # (D,D) -> (1,D,D)
        cur = imgs[0]
        for img in imgs[1:]:
            if img.fname == cur.fname and img.offset == cur.offset + 4*np.product(cur.shape):
                cur.shape = (cur.shape[0] + 1, *cur.shape[1:])
            else:
                ret.append(cur)
                cur = img
        ret.append(cur)
        return ret

    def downsample_images(imgs):
        if lazy:
            imgs = _combine_imgs(imgs)
            imgs = np.concatenate([i.get() for i in imgs])
        with Pool(min(args.max_threads, mp.cpu_count())) as p:
            oldft = np.asarray(p.map(fft.ht2_center, imgs))
            newft = oldft[:, start:stop, start:stop]
            new = np.asarray(p.map(fft.iht2_center, newft))
        return new

    def downsample_in_batches(old, b):
        new = np.empty((len(old), D, D), dtype=np.float32)
        for ii in range(math.ceil(len(old)/b)):
            log(f'Processing batch {ii}')
            new[ii*b:(ii+1)*b,:,:] = downsample_images(old[ii*b:(ii+1)*b])
        return new

    ### Downsample volume ###
    if args.is_vol:
        oldft = fft.htn_center(old)
        log(oldft.shape)
        newft = oldft[start:stop,start:stop,start:stop]
        log(newft.shape)
        new = fft.ihtn_center(newft).astype(np.float32)
        log(f'Saving {args.o}')
        mrc.write(args.o, new, is_vol=True)

    ### Downsample images ###
    elif args.chunk is None:
        new = downsample_in_batches(old, args.b)
        log(new.shape)
        log('Saving {}'.format(args.o))
        mrc.write(args.o, new.astype(np.float32), is_vol=False)

    ### Downsample images, saving chunks of N images ###
    else:
        nchunks = math.ceil(len(old)/args.chunk)
        out_mrcs = ['.{}'.format(i).join(os.path.splitext(args.o)) for i in range(nchunks)]
        chunk_names = [os.path.basename(x) for x in out_mrcs]
        for i in range(nchunks):
            log('Processing chunk {}'.format(i))
            chunk = old[i*args.chunk:(i+1)*args.chunk]
            new = downsample_in_batches(chunk, args.b)
            log(new.shape)
            log(f'Saving {out_mrcs[i]}')
            mrc.write(out_mrcs[i], new, is_vol=False)
        # Write a text file with all chunks
        out_txt = '{}.txt'.format(os.path.splitext(args.o)[0])
        log(f'Saving {out_txt}')
        with open(out_txt,'w') as f:
            f.write('\n'.join(chunk_names))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    main(add_args(parser).parse_args())
