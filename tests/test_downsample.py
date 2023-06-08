import argparse
import os
import os.path
import torch
import pytest
from cryodrgn.source import ImageSource
from cryodrgn.commands import downsample

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture
def input_star():
    # starfile refers to 13 random shuffled indices into mrcs between 0-1000
    return f"{DATA_FOLDER}/toy_projections_13.star"


def test_downsample(input_star):
    os.makedirs("output/downsampled", exist_ok=True)

    # Note - no filtering is possible in downsample currently
    args = downsample.add_args(argparse.ArgumentParser()).parse_args(
        [
            input_star,  # 13 particles
            "-D",
            "28",  # downsampled from 30x30
            "--datadir",
            DATA_FOLDER,  # If specified, prefixed to each _rlnImageName in starfile
            "-o",
            "output/downsampled/downsampled.mrcs",
        ]
    )
    downsample.main(args)

    output_data = ImageSource.from_file(
        "output/downsampled/downsampled.mrcs", lazy=False
    ).images()
    assert isinstance(output_data, torch.Tensor)
    assert output_data.shape == (13, 28, 28)


def test_downsample_in_chunks(input_star):
    os.makedirs("output", exist_ok=True)

    # TODO - no filtering is possible in downsample currently
    args = downsample.add_args(argparse.ArgumentParser()).parse_args(
        [
            input_star,  # 13 particles
            "-D",
            "28",
            "--datadir",
            DATA_FOLDER,  # If specified, prefixed to each _rlnImageName in starfile
            "--chunk",
            "5",
            "-o",
            "output/downsampled/downsampled.mrcs",
        ]
    )
    downsample.main(args)

    assert ImageSource.from_file(
        "output/downsampled/downsampled.txt", lazy=False
    ).shape == (13, 28, 28)

    assert ImageSource.from_file("output/downsampled/downsampled.0.mrcs").shape == (
        5,
        28,
        28,
    )
    assert ImageSource.from_file("output/downsampled/downsampled.1.mrcs").shape == (
        5,
        28,
        28,
    )
    assert ImageSource.from_file("output/downsampled/downsampled.2.mrcs").shape == (
        3,
        28,
        28,
    )
