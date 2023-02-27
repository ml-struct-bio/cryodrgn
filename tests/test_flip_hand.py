import os.path
import argparse
import numpy as np
import pytest
from cryodrgn.source import ImageSource
from cryodrgn.commands_utils import flip_hand

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture
def mrcs_data():
    return ImageSource.from_file(f"{DATA_FOLDER}/toy_projections.mrc").images()


def test_flip_hand(mrcs_data):
    args = flip_hand.add_args(argparse.ArgumentParser()).parse_args(
        [
            f"{DATA_FOLDER}/toy_projections.mrc",
            "-o",
            "output/toy_projections_flipped.mrc",
        ]
    )
    flip_hand.main(args)

    flipped_data = ImageSource.from_file("output/toy_projections_flipped.mrc").images()
    # torch doesn't let us use a -ve stride, hence the conversion below
    assert np.allclose(np.array(flipped_data.cpu()), np.array(mrcs_data.cpu())[::-1])
