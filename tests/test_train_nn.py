import os.path
import argparse
import numpy as np
import torch
import pytest
from cryodrgn import dataset, mrc
from cryodrgn.source import ImageSource
from cryodrgn.commands import train_nn

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture
def mrcs_data():
    return ImageSource.from_mrc(f"{DATA_FOLDER}/toy_projections.mrcs").images()


@pytest.fixture
def poses_file():
    return f"{DATA_FOLDER}/toy_angles.pkl"


def test_train_nn(mrcs_data, poses_file):
    args = train_nn.add_args(argparse.ArgumentParser()).parse_args(
        [
            f"{DATA_FOLDER}/toy_projections.mrcs",
            "--outdir",
            "output/train_nn",
            "--poses",
            poses_file,
            "--num-epochs",
            "3",
            "--no-amp",
        ]
    )
    train_nn.main(args)
