import os.path
import argparse
import pytest
from cryodrgn.source import ImageSource
from cryodrgn.commands import pc_traversal

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture
def mrcs_data():
    return ImageSource.from_file(f"{DATA_FOLDER}/hand.mrcs").images()


def test_invert_contrast(mrcs_data):
    args = pc_traversal.add_args(argparse.ArgumentParser()).parse_args(
        [f"{DATA_FOLDER}/het_z.pkl", "-o", "output/pc_traversal"]
    )
    pc_traversal.main(args)
