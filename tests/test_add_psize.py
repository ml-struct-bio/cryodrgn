import os.path
import argparse
import torch
import pytest
from cryodrgn.source import ImageSource
from cryodrgn.commands_utils import add_psize

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture
def mrcs_data():
    return ImageSource.from_file(f"{DATA_FOLDER}/toy_projections.mrcs").images()


def test_add_psize(mrcs_data):
    args = add_psize.add_args(argparse.ArgumentParser()).parse_args(
        [
            f"{DATA_FOLDER}/toy_projections.mrc",
            "-o",
            "output/toy_projections_added_psize.mrc",
        ]
    )
    add_psize.main(args)

    # Data is unchanged
    new_data = ImageSource.from_file("output/toy_projections_added_psize.mrc").images()
    assert torch.allclose(new_data, mrcs_data)
