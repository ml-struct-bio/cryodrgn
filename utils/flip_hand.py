'''Flip handedness of .mrc file'''

import argparse
import numpy as np
import sys, os
import pickle

from cryodrgn import utils
from cryodrgn import mrc
log = utils.log 

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='Input volume (.mrc)')
    parser.add_argument('-o', help='Output volume (.mrc)')
    return parser

def main(args):
    assert args.input.endswith('.mrc'), "Input volume must be .mrc file"
    assert args.o.endswith('.mrc'), "Output volume must be .mrc file"
    x, h = mrc.parse_mrc(args.input)
    x = x[::-1]
    mrc.write(args.o, x, header=h)
    log(f'Wrote {args.o}')

if __name__ == '__main__':
    main(parse_args().parse_args())
