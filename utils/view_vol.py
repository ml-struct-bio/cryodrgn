'''View 2D projection of a MRC volume'''

import argparse
import numpy as np
import sys, os

#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0,'{}/../lib-python'.format(os.path.dirname(os.path.abspath(__file__))))
import utils
import mrc

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='Input')
    parser.add_argument('-o', help='Output')
    return parser

def main(args):
    vol,_,_ = mrc.parse_mrc(args.input)
    plt.imshow(vol.sum(axis=0))
    plt.show()

if __name__ == '__main__':
    main(parse_args().parse_args())
