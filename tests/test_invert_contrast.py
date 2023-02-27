import os.path
import argparse
import torch
import pytest
from cryodrgn.source import ImageSource
from cryodrgn.commands_utils import invert_contrast

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture
def mrcs_data():
    return ImageSource.from_file(f"{DATA_FOLDER}/toy_projections.mrcs").images()


def test_invert_contrast(mrcs_data):
    args = invert_contrast.add_args(argparse.ArgumentParser()).parse_args(
        [
            f"{DATA_FOLDER}/toy_projections.mrc",
            "-o",
            "output/toy_projections_inverted.mrc",
        ]
    )
    invert_contrast.main(args)

    inverted_data = ImageSource.from_file(
        "output/toy_projections_inverted.mrc"
    ).images()
    assert torch.allclose(inverted_data, -mrcs_data)
