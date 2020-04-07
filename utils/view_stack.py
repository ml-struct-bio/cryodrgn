'''View the top 9 images in a particle stack'''

import argparse
import numpy as np
import sys, os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0,'{}/../lib-python'.format(os.path.dirname(os.path.abspath(__file__))))
import utils
import mrc
import analysis

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='Particle stack')
    parser.add_argument('-o', help='Output PNG')
    return parser

def main(args):
    stack, _ = mrc.parse_mrc(args.input,lazy=True)
    print('{} {}x{} images'.format(len(stack), *stack[0].get().shape))
    stack = [stack[x].get() for x in range(9)]
    analysis.plot_projections(stack)
    if args.o:
        plt.savefig(args.o)
    else:
        plt.show()

if __name__ == '__main__':
    main(parse_args().parse_args())
