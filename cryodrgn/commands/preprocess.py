'''
Preprocess a dataset for more streamlined cryoDRGN training
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
    parser.add_argument('-o', metavar='MRCS', type=os.path.abspath, required=True, help='Output .mrcs')
    parser.add_argument('--datadir', help='Optionally provide path to input .mrcs if loading from a .star or .cs file')

    group = parser.add_argument_group('Image preprocessing settings')
    group.add_argument('--ind', type=os.path.abspath, metavar='PKL', help='Filter particle stack by these indices')
    group.add_argument('-D', type=int, help='New box size in pixels (if downsampling), must be even')
    group.add_argument('--uninvert-data', dest='invert_data', action='store_false', help='Do not invert data sign')
    group.add_argument('--window-r', default=0.85, type=float, help='Circular windowing mask inner radius (default: %(default)s)')
    group.add_argument('--no-window', dest='window', action='store_false', help='Turn off real space windowing of dataset')

    group = parser.add_argument_group('Extra arguments for volume generation')
    group.add_argument('-b', type=int, default=5000, help='Batch size for processing images (default: %(default)s)')
    group.add_argument('--chunk', type=int, default=100000, help='Chunksize (in # of images) to split particle stack when saving')
    group.add_argument('--no-lazy', dest='lazy', action='store_false', help='Load whole dataset (faster than loading in batches)')
    group.add_argument('--max-threads', type=int, default=16, help='Maximum number of CPU cores for parallelization (default: %(default)s)')
    return parser

def mkbasedir(out):
    if not os.path.exists(os.path.dirname(out)):
        os.makedirs(os.path.dirname(out))

def warnexists(out):
    if os.path.exists(out):
        log(f'Warning: {out} already exists. Overwriting.')

def main(args):
    mkbasedir(args.o)
    warnexists(args.o)
    assert (args.o.endswith('.mrcs') or args.o.endswith('.txt')), "Must specify output in .mrcs file format"

    # load images
    lazy = args.lazy
    images = dataset.load_particles(args.mrcs, lazy=lazy, datadir=args.datadir)
    
    # filter images
    if args.ind is not None: 
        log(f'Filtering image dataset with {args.ind}')
        ind = utils.load_pkl(args.ind).astype(int)
        images = [images[i] for i in ind] if lazy else images[ind]

    original_D = images[0].get().shape[0] if lazy else images.shape[-1]
    log(f'Loading {len(images)} {original_D}x{original_D} images')
    window = args.window
    invert_data = args.invert_data
    downsample = (args.D and args.D < original_D)
    if downsample:
        assert args.D <= original_D, f'New box size {args.D} cannot be larger than the original box size {D}'
        assert args.D % 2 == 0, 'New box size must be even'
        start = int(original_D/2 - args.D/2)
        stop = int(original_D/2 + args.D/2)
        D = args.D
        log(f'Downsampling images to {D}x{D}')
    else:
        D = original_D

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

    def preprocess(imgs):
        if lazy:
            imgs = _combine_imgs(imgs)
            imgs = np.concatenate([i.get() for i in imgs])
        with Pool(min(args.max_threads, mp.cpu_count())) as p:
            # todo: refactor as a routine in dataset.py

            # note: applying the window before downsampling is slightly 
            # different than in the original workflow
            if window:
                imgs *= dataset.window_mask(original_D, args.window_r, .99)
            ret = np.asarray(p.map(fft.ht2_center, imgs))
            if invert_data:
                ret *= -1
            if downsample:
                ret = ret[:, start:stop, start:stop]
            ret = fft.symmetrize_ht(ret)
        return ret

    def preprocess_in_batches(imgs, b):
        ret = np.empty((len(imgs), D+1, D+1), dtype=np.float32)
        Nbatches = math.ceil(len(imgs)/b)
        for ii in range(Nbatches):
            log(f'Processing batch of {b} images ({ii+1} of {Nbatches})')
            ret[ii*b:(ii+1)*b,:,:] = preprocess(imgs[ii*b:(ii+1)*b])
        return ret

    nchunks = math.ceil(len(images)/args.chunk)
    out_mrcs = [f'.{i}.ft'.join(os.path.splitext(args.o)) for i in range(nchunks)]
    chunk_names = [os.path.basename(x) for x in out_mrcs]
    for i in range(nchunks):
        log(f'Processing chunk {i+1} of {nchunks}')
        chunk = images[i*args.chunk:(i+1)*args.chunk]
        new = preprocess_in_batches(chunk, args.b)
        log(f'New shape: {new.shape}')
        log(f'Saving {out_mrcs[i]}')
        mrc.write(out_mrcs[i], new, is_vol=False)

    out_txt = f'{os.path.splitext(args.o)[0]}.ft.txt'
    log(f'Saving summary txt file {out_txt}')
    with open(out_txt,'w') as f:
        f.write('\n'.join(chunk_names))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    main(add_args(parser).parse_args())

