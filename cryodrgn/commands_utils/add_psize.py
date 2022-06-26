'''Add pixel size to header of .mrc file'''

import argparse
import numpy as np
import sys, os
import pickle

from cryodrgn import utils
from cryodrgn import mrc
log = utils.log 

def add_args(parser):
    parser.add_argument('input', help='Input volume (.mrc)')
    parser.add_argument('--Apix', type=float, default=1, help='Angstrom/pixel (default: %(default)s)')
    parser.add_argument('-o', help='Output volume (.mrc)')
    return parser

def main(args):
    assert args.input.endswith('.mrc'), "Input volume must be .mrc file"
    assert args.o.endswith('.mrc'), "Output volume must be .mrc file"
    x, h = mrc.parse_mrc(args.input)
    h.update_apix(args.Apix)
    mrc.write(args.o, x, header=h)
    log(f'Wrote {args.o}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
