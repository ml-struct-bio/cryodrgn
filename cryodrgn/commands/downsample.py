'''
Downsample an image stack or volume by clipping fourier frequencies
'''

import argparse
import numpy as np
import sys, os
import math

from cryodrgn import utils
from cryodrgn import mrc
from cryodrgn import fft
from cryodrgn import dataset

log = utils.log

def add_args(parser):
    parser.add_argument('mrcs', help='Input particles or volume (.mrc, .mrcs, .star, or .txt)')
    parser.add_argument('-D', type=int, required=True, help='New box size in pixels, must be even')
    parser.add_argument('-o', metavar='MRCS', type=os.path.abspath, required=True, help='Output projection stack (.mrcs)')
    parser.add_argument('--is-vol',action='store_true', help='Flag if input .mrc is a volume')
    parser.add_argument('--chunk', type=int, help='Chunksize (in # of images) to split particle stack when saving')
    parser.add_argument('--datadir', help='Optionally provide path to input .mrcs if loading from a .star or .cs file')
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

    old = dataset.load_particles(args.mrcs, lazy=True, datadir=args.datadir)
    oldD = old[0].get().shape[0]
    assert args.D <= oldD, f'New box size {args.D} cannot be larger than the original box size {oldD}'
    assert args.D % 2 == 0, 'New box size must be even'
    
    D = args.D
    start = int(oldD/2 - D/2)
    stop = int(oldD/2 + D/2)

    ### Downsample volume ###
    if args.is_vol:
        oldft = fft.htn_center(np.array([x.get() for x in old]))
        log(oldft.shape)
        newft = oldft[start:stop,start:stop,start:stop]
        log(newft.shape)
        new = fft.ihtn_center(newft).astype(np.float32)
        log(f'Saving {args.o}')
        mrc.write(args.o, new, is_vol=True)

    ### Downsample images ###
    elif args.chunk is None:
        new = []
        for i in range(len(old)):
            if i % 1000 == 0:
                log(f'Processing image {i} of {len(old)}')
            img = old[i]
            oldft = fft.ht2_center(img.get()).astype(np.float32)
            newft = oldft[start:stop, start:stop]
            new.append(fft.ihtn_center(newft).astype(np.float32))
        assert oldft[int(oldD/2),int(oldD/2)] == newft[int(D/2),int(D/2)]
        new = np.asarray(new)
        log(new.shape)
        log('Saving {}'.format(args.o))
        mrc.write(args.o, new, is_vol=False)

    ### Downsample images, saving chunks of N images ###
    else:
        chunk_names = []
        nchunks = math.ceil(len(old)/args.chunk)
        for i in range(nchunks):
            log('Processing chunk {}'.format(i))
            out_mrcs = '.{}'.format(i).join(os.path.splitext(args.o))
            new = []
            for img in old[i*args.chunk:(i+1)*args.chunk]:
                oldft = fft.ht2_center(img.get()).astype(np.float32)
                newft = oldft[start:stop, start:stop]
                new.append(fft.ihtn_center(newft).astype(np.float32))
            assert oldft[int(oldD/2),int(oldD/2)] == newft[int(D/2),int(D/2)]
            new = np.asarray(new)
            log(new.shape)
            log(f'Saving {out_mrcs}'.format(out_mrcs))
            mrc.write(out_mrcs, new, is_vol=False)
            chunk_names.append(os.path.basename(out_mrcs))
        # Write a text file with all chunks
        out_txt = '{}.txt'.format(os.path.splitext(args.o)[0])
        log(f'Saving {out_txt}')
        with open(out_txt,'w') as f:
            f.write('\n'.join(chunk_names))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    main(add_args(parser).parse_args())
