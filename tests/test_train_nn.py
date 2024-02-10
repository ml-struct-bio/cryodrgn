"""Integration tests of neural network reconstruction with known poses.

The training done here has unrealistically low parameter values to allow the tests to
run in any environment in a reasonable amount of time with or without GPUs.

"""
import os
import pytest

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.mark.parametrize(
    "train_dir", [{"train_cmd": "train_nn", "epochs": 5}], indirect=True
)
@pytest.mark.parametrize("load_epoch", [0, 2])
@pytest.mark.parametrize("train_epochs", [4, 5, 6])
def test_reconstruct_and_from_checkpoint(trained_dir, load_epoch, train_epochs):
    trained_dir.train_load_epoch(load_epoch, train_epochs)
