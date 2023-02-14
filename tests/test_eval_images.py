import os.path
import argparse
import torch
import pytest
from cryodrgn.source import ImageSource
from cryodrgn.commands import eval_images

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture
def mrcs_data():
    return ImageSource.from_mrcs(f"{DATA_FOLDER}/hand.mrcs").images()


def test_invert_contrast(mrcs_data):
    args = eval_images.add_args(argparse.ArgumentParser()).parse_args(
        [
            f"{DATA_FOLDER}/hand.mrcs",
            f"{DATA_FOLDER}/het_weights.pkl",
            "--config",
            f"{DATA_FOLDER}/het_config.pkl",
            "-o",
            f"output/out_eval_images_losses.pkl",
            "--out-z",
            f"output/out_eval_images_z.pkl",
            "--poses",
            f"{DATA_FOLDER}/hand_rot.pkl",
            "--log-interval",
            "1",
            "--verbose",
        ]
    )
    eval_images.main(args)
