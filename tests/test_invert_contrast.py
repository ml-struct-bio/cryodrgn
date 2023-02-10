import os.path
import argparse
import numpy as np
import torch
import pytest
from cryodrgn import dataset, mrc
from cryodrgn.source import ImageSource
from cryodrgn.commands_utils import invert_contrast

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture
def mrcs_data():
    return ImageSource.from_mrcs(f"{DATA_FOLDER}/toy_projections.mrcs").images()


def test_invert_contrast(mrcs_data):
    args = invert_contrast.add_args(argparse.ArgumentParser()).parse_args(
        [
            f"{DATA_FOLDER}/toy_projections.mrc",
            "-o",
            "output/toy_projections_inverted.mrc",
        ]
    )
    invert_contrast.main(args)

    inverted_data = ImageSource.from_mrcs(
        "output/toy_projections_inverted.mrc"
    ).images()
    assert np.allclose(inverted_data, -mrcs_data)
