"""Integration tests of ab initio reconstruction and downstream analyses.

Note that the training done here has unrealistically low parameter values to allow the
tests to run in any environment in a reasonable amount of time with or without GPUs.

"""
import pytest
import os
import argparse
from cryodrgn.commands import backproject_voxel, train_nn


@pytest.mark.parametrize(
    "particles", ["toy.mrcs", "toy.star", "toy.txt"], indirect=True
)
@pytest.mark.parametrize("poses", ["toy-poses", "toy-angles"], indirect=True)
@pytest.mark.parametrize("batch_size", ["24", "36"], ids=("batch24", "batch36"))
@pytest.mark.parametrize("use_amp", [False, True], ids=("no-amp", "use-amp"))
def test_homogeneous_with_poses(outdir, particles, poses, batch_size, use_amp):
    args = [
        particles,
        "-o",
        outdir,
        "-n",
        "10",
        "--poses",
        poses,
        "--dim",
        "12",
        "-b",
        batch_size,
    ]
    if not use_amp:
        args += ["--no-amp"]

    train_nn.main(train_nn.add_args(argparse.ArgumentParser()).parse_args(args))
    assert "weights.9.pkl" in os.listdir(outdir)


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


@pytest.mark.parametrize(
    "particles, poses, ctf",
    [
        ("toy.mrcs", "toy-poses", None),
        ("hand", "hand-rot", None),
        ("hand", "hand-poses", None),
        ("tilts.star", "tilt-poses", "CTF-Tilt"),
    ],
    indirect=True,
)
def test_backprojection(outdir, particles, poses, ctf):
    args = [particles, "--poses", poses, "-o", os.path.join(outdir, "vol.mrc")]

    if ctf is not None:
        args += ["--ctf", ctf]
    if "tilt" in particles:
        args += ["--tilt", "-d", "2.93"]

    parser = argparse.ArgumentParser()
    backproject_voxel.add_args(parser)
    backproject_voxel.main(parser.parse_args(args))
