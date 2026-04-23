"""tSNE dimensionality reduction"""

import argparse
import pickle
import logging
from sklearn.manifold import TSNE

logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Input z.pkl")
    parser.add_argument("--stride", type=int, help="Stride the dataset")
    parser.add_argument(
        "-p", default=1000.0, type=float, help="Perplexity (default: %(default)s)"
    )
    parser.add_argument("-o", help="Output pickle")
    return parser


def main(args):
    with open(args.input, "rb") as f:
        z = pickle.load(f)

    if args.stride:
        z = z[:: args.stride]
        logger.info(
            f"Loaded zdim={z.shape[1]} latent space and used "
            f"striding to reduce to {z.shape[0]} datapoints"
        )
    else:
        logger.info(
            f"Loaded zdim={z.shape[1]} latent space with {z.shape[0]} datapoints"
        )

    logger.info("Fitting t-SNE...")
    z_embedded = TSNE(n_components=2, perplexity=args.p).fit_transform(z)
    with open(args.o, "wb") as f:
        pickle.dump(z_embedded, f)


if __name__ == "__main__":
    main(parse_args().parse_args())
