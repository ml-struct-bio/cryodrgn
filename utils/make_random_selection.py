'''Get an index selection (.pkl) for a random subset of particles'''

import argparse
import numpy as np
import sys, os
import pickle

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('N', type=int, help='Total number of particles')
    parser.add_argument('-o', required=True, help='Output selection (.pkl)')
    parser.add_argument('-s', help='Optionally save out inverted selection (.pkl)')
    parser.add_argument('--frac', type=float, help='Fraction of particles to select')
    parser.add_argument('-n', type=int, help='Number of particles to select')
    parser.add_argument('--seed', type=int, default=0, help='Random seed (default: %(default)s)')
    return parser

def main(args):
    print(f'{args.N} total particles')
    np.random.seed(args.seed)
    ind = np.arange(args.N)
    assert bool(args.frac) != bool(args.n), "Must specify --frac or -n"
    n = int(args.N*args.frac) if args.frac else args.n
    train = np.random.choice(ind, n, replace=False)
    train = np.array(sorted(train))
    test = set(ind)-set(train)
    test = np.array(sorted(test))
    
    print(f'{len(train)} particles in selection: {train}')
    print(f'Saving {args.o}')
    pickle.dump(train, open(args.o,'wb'))

    if args.s is not None:
        print(f'{len(test)} particles in inverted selection: {test}')
        print(f'Saving {args.s}')
        pickle.dump(test, open(args.s,'wb'))

if __name__ == '__main__':
    main(parse_args().parse_args())
