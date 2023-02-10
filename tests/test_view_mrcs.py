import os.path
import argparse
import numpy as np
import torch
import pytest
from cryodrgn import dataset, mrc
from cryodrgn.source import ImageSource
from cryodrgn.commands_utils import view_mrcs

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture
def mrcs_data():
    return ImageSource.from_mrcs(f"{DATA_FOLDER}/toy_projections.mrcs").images()


def test_invert_contrast(mrcs_data):
    args = view_mrcs.add_args(argparse.ArgumentParser()).parse_args(
        [f"{DATA_FOLDER}/toy_projections.mrc", "--invert"]
    )
    view_mrcs.main(args)
