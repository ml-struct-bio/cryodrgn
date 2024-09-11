"""Create an index corresponding to the selection of a random subset of particles.

Example usage
-------------
# Sample a thousand indices from [0 ... 189042]
$ cryodrgn_utils select_random 189043 -o my-indices.pkl -n 100000

# Sample half of the indices from [0 ... 189042]
$ cryodrgn_utils select_random 189043 -o my-indices.pkl --frac 0.5

"""
import argparse
import pickle
import numpy as np


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("N", type=int, help="Total number of particles")
    parser.add_argument("-o", required=True, help="Output selection (.pkl)")
    parser.add_argument("-n", type=int, help="Number of particles to select")
    parser.add_argument("-s", help="Optionally save out inverted selection (.pkl)")
    parser.add_argument(
        "--frac", type=float, help="Optionally specify fraction of particles to select"
    )
    parser.add_argument(
        "--seed", type=int, default=0, help="Random seed (default: %(default)s)"
    )


def main(args: argparse.Namespace) -> None:
    print(f"{args.N} total particles")
    np.random.seed(args.seed)
    ind = np.arange(args.N)
    assert bool(args.frac) != bool(args.n), "Must specify --frac or -n"
    n = int(args.N * args.frac) if args.frac else args.n
    train = np.random.choice(ind, n, replace=False)
    train = np.array(sorted(train))
    test = set(ind) - set(train)
    test = np.array(sorted(test))

    print(f"{len(train)} particles in selection: {train}")
    print(f"Saving {args.o}")
    pickle.dump(train, open(args.o, "wb"))

    if args.s is not None:
        print(f"{len(test)} particles in inverted selection: {test}")
        print(f"Saving {args.s}")
        pickle.dump(test, open(args.s, "wb"))
