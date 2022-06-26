'''Select particle or volume data based on (kmeans) cluster labels'''

import argparse
import numpy as np
import sys, os
import pickle

from cryodrgn import utils, analysis
log = utils.log 

def add_args(parser):
    parser.add_argument('labels', help='Input labels.pkl')
    parser.add_argument('--sel', nargs='+', type=int, help='Ids of clusters to select')
    parser.add_argument('-o', type=os.path.abspath, help='Output particle index selection (.pkl)')

    group = parser.add_argument_group('Get original particle selection (if trained on a subset of the dataset with --ind)')
    group.add_argument('--parent-ind', type=os.path.abspath, help='Parent index .pkl')
    group.add_argument('--N-orig', type=int, help='Number of particles in original dataset')
    return parser

def main(args):
    labels = utils.load_pkl(args.labels)
    log(f'{len(labels)} particles')
    log(f'Selecting clusters {args.sel}')
    ind = analysis.get_ind_for_cluster(labels, args.sel)
    log(f'Selected {len(ind)} particles')
    log(ind)
    if args.parent_ind is not None:
        log('Converting to original indices')
        parent_ind = utils.load_pkl(args.parent_ind)
        assert args.N_orig
        ind = convert_original_indices(ind, N_orig, parent_ind)
        log(ind)
    utils.save_pkl(ind, args.o)
    log(f'Saved {args.o}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
