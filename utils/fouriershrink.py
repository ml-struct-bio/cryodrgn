'''
Resize an image stack or volume by clipping fourier frequencies
'''

import argparse
import numpy as np
import sys, os

sys.path.insert(0,'{}/../lib-python'.format(os.path.dirname(os.path.abspath(__file__))))
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
    parser.add_argument('--is-vol',action='store_true')
    parser.add_argument('-D', type=int, required=True, help='New image size, must be even')
    parser.add_argument('--chunk', type=int, help='Chunksize (# images) to split projection stack')
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

    old, _, _ = mrc.parse_mrc(args.mrcs,lazy=True)
    oldD = old[0].get().shape[0]
    assert args.D < oldD
    assert args.D % 2 == 0
    
    D = args.D

    start = int(oldD/2 - D/2)
    stop = int(oldD/2 + D/2)

    if args.is_vol:
        oldft = fft.htn_center(np.array([x.get() for x in old]))
        log(oldft.shape)
        newft = oldft[start:stop,start:stop,start:stop]
        log(newft.shape)
        new = fft.ihtn_center(newft).astype(np.float32)

    if args.chunk is None:
        new = []
        for img in old:
            oldft = fft.ht2_center(img.get()).astype(np.float32)
            newft = oldft[start:stop, start:stop]
            new.append(fft.ihtn_center(newft).astype(np.float32))
        assert oldft[int(oldD/2),int(oldD/2)] == newft[int(D/2),int(D/2)]
        new = np.asarray(new)
        log(new.shape)
        log('Saving {}'.format(args.o))
        mrc.write(args.o,new)
    else:
        nchunks = len(old) // args.chunk + 1
        for i in range(nchunks):
            log('Processing chunk {}'.format(i))
            out = '.{}'.format(i).join(os.path.splitext(args.o))
            new = []
            for img in old[i*args.chunk:(i+1)*args.chunk]:
                oldft = fft.ht2_center(img.get()).astype(np.float32)
                newft = oldft[start:stop, start:stop]
                new.append(fft.ihtn_center(newft).astype(np.float32))
            assert oldft[int(oldD/2),int(oldD/2)] == newft[int(D/2),int(D/2)]
            new = np.asarray(new)
            log(new.shape)
            log('Saving {}'.format(out))
            mrc.write(out,new)

if __name__ == '__main__':
    main(parse_args().parse_args())
