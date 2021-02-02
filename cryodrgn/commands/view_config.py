'''
Display config information of a cryoDRGN job
'''

import argparse
import numpy as np
import sys, os
import pickle
import subprocess
from datetime import datetime as dt

from cryodrgn import utils
from pprint import pprint

log = utils.log

def add_args(parser):
    parser.add_argument('workdir', type=os.path.abspath, help='Directory with cryoDRGN results')
    return parser

def main(args):
    f = open(f'{args.workdir}/config.pkl','rb')
    cfg = pickle.load(f)
    try:
        meta = pickle.load(f)
        log(f'Version: {meta["version"]}')
        log(f'Creation time: {meta["time"]}')
        log('Command:')
        print(' '.join(meta['cmd']))
    except:
        pass
    log('Config:')
    pprint(cfg)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
