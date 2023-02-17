import argparse
import os.path
import pytest
import torch
from cryodrgn.commands import preprocess
from cryodrgn.source import ImageSource

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

    data = ImageSource.from_file("output/preprocessed.ft.txt").images()
    assert isinstance(data, torch.Tensor)
    assert data.shape == (100, 21, 21)


def test_preprocess_lazy(mrcs_file):
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
            "-b",
            "10",
            "--chunk",
            "3",
        ]
    )
    preprocess.main(args)

    data = ImageSource.from_file("output/preprocessed.ft.txt").images()
    assert isinstance(data, torch.Tensor)
    assert data.shape == (100, 21, 21)
