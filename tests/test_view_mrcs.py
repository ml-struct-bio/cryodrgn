import pytest
import os.path
import argparse
from unittest.mock import patch
from cryodrgn.commands_utils import view_mrcs


@patch("matplotlib.pyplot.show")
def test_invert_contrast(mock_pyplot_show):
    parser = argparse.ArgumentParser()
    view_mrcs.add_args(parser)
    view_mrcs.main(
        parser.parse_args(
            [os.path.join(pytest.DATADIR, "toy_projections.mrc"), "--invert"]
        )
    )
