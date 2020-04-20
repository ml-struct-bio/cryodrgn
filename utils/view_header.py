'''View the header of a .mrc or .mrcs file'''

import argparse
import numpy as np
import sys, os

from cryodrgn import mrc
from cryodrgn import utils

log = utils.log

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='Particle stack')
    parser.add_argument('-o', help='Output PNG')
    return parser

def main(args):
    if not (args.input.endswith('.mrc') or args.input.endswith('.mrcs')):
        log(f'Warning: {args.input} does not appear to be a .mrc(s) file')
    header = mrc.parse_header(args.input)
    print(header)

if __name__ == '__main__':
    main(parse_args().parse_args())
