import os.path
import argparse
from unittest.mock import patch
from cryodrgn.commands_utils import view_mrcs

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@patch("matplotlib.pyplot.show")
def test_invert_contrast(mock_pyplot_show):
    args = view_mrcs.add_args(argparse.ArgumentParser()).parse_args(
        [f"{DATA_FOLDER}/toy_projections.mrc", "--invert"]
    )
    view_mrcs.main(args)
