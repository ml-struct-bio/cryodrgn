"""View the header metadata of a .mrc or .mrcs file"""

import argparse
from pprint import pprint
import logging
from cryodrgn.source import MRCHeader

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input", help="Particle stack (.mrcs) or density map (.mrc)")


def main(args: argparse.Namespace) -> None:
    if not (args.input.endswith(".mrc") or args.input.endswith(".mrcs")):
        logger.warning(f"Warning: {args.input} does not appear to be a .mrc(s) file")
    header = MRCHeader.parse(args.input)
    pprint(header.fields)
    pprint(header.extended_header)
