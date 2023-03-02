import os.path
import argparse
from cryodrgn.commands_utils import view_header

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


def test_view_header():
    args = view_header.add_args(argparse.ArgumentParser()).parse_args(
        [
            f"{DATA_FOLDER}/hand.mrcs",
        ]
    )
    view_header.main(args)
