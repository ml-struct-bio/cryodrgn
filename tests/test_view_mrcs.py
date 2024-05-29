import pytest
import os.path
import argparse
from unittest.mock import patch
from cryodrgn.commands_utils import view_mrcs


@patch("matplotlib.pyplot.show")
def test_invert_contrast(mock_pyplot_show):
    args = view_mrcs.add_args(argparse.ArgumentParser()).parse_args(
        [os.path.join(pytest.DATADIR, "toy_projections.mrc"), "--invert"]
    )
    view_mrcs.main(args)
