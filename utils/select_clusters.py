'''Select particles based on kmeans cluster labels'''

import argparse
import numpy as np
import sys, os
import pickle

from cryodrgn import utils, analysis
log = utils.log 

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('labels', help='Input labels.pkl')
    parser.add_argument('--sel', nargs='+', type=int, help='Cluster ids')
    parser.add_argument('-o', type=os.path.abspath, help='Output particle index selection (.pkl)')

    group = parser.add_argument_group('Get original particle selection (if trained on a subset of the dataset with --ind)')
    group.add_argument('--parent-ind', type=os.path.abspath, help='Parent index .pkl')
    group.add_argument('--N-orig', type=int, help='Number of particles in original dataset')
    return parser

def main(args):
    labels = utils.load_pkl(args.labels)
    print(f'{len(labels)} particles')
    print(f'Selecting clusters {args.sel}')
    ind = analysis.get_ind_for_cluster(labels, args.sel)
    print(f'Selected {len(ind)} particles')
    print(ind)
    if args.parent_ind is not None:
        parent_ind = utils.load_pkl(args.parent_ind)
        assert args.N_orig
        ind = convert_original_indices(ind, N_orig, parent_ind)
        print('Converted to original indices')
        print(ind)
    print(f'Saved {args.o}')
    utils.save_pkl(ind, args.o)

if __name__ == '__main__':
    main(parse_args().parse_args())
