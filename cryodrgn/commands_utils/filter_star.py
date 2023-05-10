"""
Filter a .star file
"""

import argparse
import os
import logging
from cryodrgn import starfile
from cryodrgn import utils

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument("input", help="Input .star file")
    parser.add_argument("--ind", required=True, help="Array of selected indices (.pkl)")
    parser.add_argument("-o", type=os.path.abspath, help="Output .star file")
    return parser


def main(args):
    s = starfile.Starfile.load(args.input)
    ind = utils.load_pkl(args.ind)
    df = s.df.loc[ind]
    s = starfile.Starfile(headers=None, df=df)
    s.write(args.o)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
