import argparse
import os.path
import pytest
import numpy as np
from cryodrgn import dataset
from cryodrgn.commands import preprocess

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture
def mrcs_file():
    return f"{DATA_FOLDER}/hand.mrcs"


def test_preprocess(mrcs_file):
    args = preprocess.add_args(argparse.ArgumentParser()).parse_args(
        [
            mrcs_file,
            "-o",
            "output/preprocessed.mrcs",
            "--uninvert-data",
            "--window-r",
            "0.5",
            "-D",
            "20",
            "--no-lazy",
            "-b",
            "10",
            "--chunk",
            "3",
        ]
    )
    preprocess.main(args)

    data = dataset.load_particles("output/preprocessed.ft.txt")
    assert isinstance(data, np.ndarray)
    assert data.shape == (100, 21, 21)
