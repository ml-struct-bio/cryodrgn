'''View the header metadata of a .mrc or .mrcs file'''

import argparse
import numpy as np
import sys, os

from cryodrgn import mrc
from cryodrgn import utils

from pprint import pprint
log = utils.log

def add_args(parser):
    parser.add_argument('input', help='Particle stack (.mrcs) or density map (.mrc)')
    return parser

def main(args):
    if not (args.input.endswith('.mrc') or args.input.endswith('.mrcs')):
        log(f'Warning: {args.input} does not appear to be a .mrc(s) file')
    header = mrc.parse_header(args.input)
    pprint(header.fields)
    pprint(header.extended_header)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
