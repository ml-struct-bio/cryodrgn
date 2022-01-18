'''Get an index selection (.pkl) for a random subset of particles'''

import argparse
import numpy as np
import sys, os
import pickle

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('N', type=int, help='Total number of particles')
    parser.add_argument('-o', required=True, help='Output basename to append .train.pkl and .test.pkl')
    parser.add_argument('--frac', type=float, help='Fraction to include in train')
    parser.add_argument('-n', type=int, help='Number of particles to include in train')
    parser.add_argument('--seed', type=int, default=0, help='Random seed')
    return parser

def main(args):
    print(f'{args.N} total particles')
    np.random.seed(args.seed)
    ind = np.arange(args.N)
    assert bool(args.frac) != bool(args.n), "Must specify --frac or -n"
    n = int(args.N*args.frac) if args.frac else args.n
    train = np.random.choice(ind, n, replace=False)
    test = set(ind)-set(train)
    test = np.array(sorted(test))
    
    print(f'{len(train)} particles in training set: {train}')
    f = '{}.train.pkl'.format(args.o)
    print(f'Saving {f}')
    pickle.dump(train, open(f,'wb'))

    print(f'{len(test)} particles in test set: {test}')
    f = '{}.test.pkl'.format(args.o)
    print(f'Saving {f}')
    pickle.dump(test, open(f,'wb'))

if __name__ == '__main__':
    main(parse_args().parse_args())
