"""Integration tests of ab initio reconstruction and downstream analyses.

The training done here has unrealistically low parameter values to allow the tests to
run in any environment in a reasonable amount of time with or without GPUs.

"""
import pytest


@pytest.mark.parametrize(
    "abinit_dir", [{"zdim": zdim} for zdim in [0, 4, 8]], indirect=True
)
def test_abinit_checkpoint_analysis_and_backproject(abinit_dir):
    abinit_dir.train()
    abinit_dir.train(load_epoch=0)
    abinit_dir.backproject()
    abinit_dir.view_config()


@pytest.mark.parametrize(
    "train_dir",
    [{"train_cmd": "train_nn", "epochs": 5}, {"train_cmd": "train_vae", "epochs": 5}],
    indirect=True,
)
@pytest.mark.parametrize("load_epoch", [0, 2])
@pytest.mark.parametrize("train_epochs", [4, 5, 6])
def test_frompose_train_and_from_checkpoint(trained_dir, load_epoch, train_epochs):
    trained_dir.train_load_epoch(load_epoch, train_epochs)
