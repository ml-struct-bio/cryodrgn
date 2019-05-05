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

def plot_projections(out_png, imgs):
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(10,10))
    axes = axes.ravel()
    for i in range(min(len(imgs),9)):
        axes[i].imshow(imgs[i])
    plt.savefig(out_png)

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='Particle stack')
    parser.add_argument('-o', help='Output PNG')
    return parser

def main(args):
    stack,_,_ = mrc.parse_mrc(args.input,lazy=True)
    stack = [stack[x].get() for x in range(9)]
    plot_projections(args.o, stack)

if __name__ == '__main__':
    main(parse_args().parse_args())
