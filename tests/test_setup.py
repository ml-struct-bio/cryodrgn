"""Volume reconstruction training with fixed poses, as well as downstream analyses."""

import pytest
import os.path
import argparse
from cryodrgn.commands import setup, train


def parse_setup_request(
    config_file,
    model,
    dataset,
    particles,
    ctf,
    poses,
    ind,
    datadir,
    reconstruction_type,
    z_dim,
    pose_estimation,
    tilt,
    cfgs,
    defaults,
):
    args = list()

    if config_file is not None:
        args += [config_file]
    if model is not None:
        args += ["--model", model]
    if dataset is not None:
        args += ["--dataset", dataset]
    if particles is not None:
        args += ["--particles", particles]
    if ctf is not None:
        args += ["--ctf", ctf]
    if poses is not None:
        args += ["--poses", poses]
    if ind is not None:
        args += ["--ind", ind]
    if datadir is not None:
        args += ["--datadir", datadir]
    if reconstruction_type is not None:
        args += ["--reconstruction-type", reconstruction_type]
    if z_dim is not None:
        args += ["--z-dim", z_dim]
    if pose_estimation is not None:
        args += ["--pose-estimation", pose_estimation]
    if tilt is not None and tilt is True:
        args += ["--tilt"]

    add_cfgs = ["num_epochs=2"]
    if cfgs:
        add_cfgs += cfgs
    args += ["--cfgs"] + add_cfgs
    if defaults:
        args += ["--defaults"]

    return args


@pytest.mark.parametrize(
    "model, ctf",
    [("hps", "CTF-Test"), ("hps", None), ("amort", "CTF-Test")],
    indirect=["ctf"],
)
@pytest.mark.parametrize("dataset", [None])
@pytest.mark.parametrize("particles", ["toy.mrcs", "toy.txt"], indirect=True)
@pytest.mark.parametrize("indices", [None, "first-100"], indirect=True)
@pytest.mark.parametrize("datadir", [None], indirect=True)
@pytest.mark.parametrize(
    "reconstruction_type, z_dim",
    [(None, None), (None, 0)],
)
@pytest.mark.parametrize(
    "pose_estimation, poses",
    [
        (None, None),
        (None, "toy_poses"),
        pytest.param(
            "fixed",
            None,
            marks=pytest.mark.xfail(
                raises=ValueError, reason="fixed estimation but no poses given"
            ),
        ),
        ("abinit", None),
        ("fixed", "toy_poses"),
        ("abinit", "toy_poses"),
    ],
    indirect=["poses"],
)
@pytest.mark.parametrize("tilt", [None])
@pytest.mark.parametrize("cfgs", [None])
@pytest.mark.parametrize("defaults", [False])
def test_setup_and_print_command(
    model,
    dataset,
    particles,
    ctf,
    poses,
    indices,
    datadir,
    reconstruction_type,
    z_dim,
    pose_estimation,
    tilt,
    cfgs,
    defaults,
):
    args = parse_setup_request(
        None,
        model,
        dataset,
        particles.path,
        ctf.path,
        poses.path,
        indices.path,
        datadir.path,
        reconstruction_type,
        z_dim,
        pose_estimation,
        tilt,
        cfgs,
        defaults,
    )

    parser = argparse.ArgumentParser()
    setup.add_args(parser)
    setup.main(parser.parse_args(args))


@pytest.mark.parametrize(
    "model, ctf",
    [("hps", "CTF-Test"), ("hps", None), ("amort", "CTF-Test")],
    indirect=["ctf"],
)
@pytest.mark.parametrize("dataset", [None])
@pytest.mark.parametrize("particles", ["toy.mrcs", "toy.txt"], indirect=True)
@pytest.mark.parametrize("indices", [None, "first-100"], indirect=True)
@pytest.mark.parametrize("datadir", [None], indirect=True)
@pytest.mark.parametrize(
    "reconstruction_type, z_dim",
    [(None, None), (None, 0)],
)
@pytest.mark.parametrize(
    "pose_estimation, poses",
    [
        (None, None),
        (None, "toy_poses"),
        pytest.param(
            "fixed",
            None,
            marks=pytest.mark.xfail(
                raises=ValueError, reason="fixed estimation but no poses given"
            ),
        ),
        ("abinit", None),
        ("fixed", "toy_poses"),
        ("abinit", "toy_poses"),
    ],
    indirect=["poses"],
)
@pytest.mark.parametrize("tilt", [None])
@pytest.mark.parametrize("cfgs", [None])
@pytest.mark.parametrize("defaults", [False])
def test_setup_and_run_command(
    tmpdir_factory,
    model,
    dataset,
    particles,
    ctf,
    poses,
    indices,
    datadir,
    reconstruction_type,
    z_dim,
    pose_estimation,
    tilt,
    cfgs,
    defaults,
):
    tmpdir = tmpdir_factory.mktemp("setup_and_run_commands")
    config_file = os.path.join(tmpdir, "configs.yaml")
    outdir = os.path.join(tmpdir, "out")
    args = parse_setup_request(
        config_file,
        model,
        dataset,
        particles.path,
        ctf.path,
        poses.path,
        indices.path,
        datadir.path,
        reconstruction_type,
        z_dim,
        pose_estimation,
        tilt,
        cfgs,
        defaults,
    )

    parser = argparse.ArgumentParser()
    setup.add_args(parser)
    setup.main(parser.parse_args(args))
    parser = argparse.ArgumentParser()
    train.add_args(parser)
    train.main(parser.parse_args([config_file, "-o", outdir]))

    out_files = os.listdir(outdir)
    assert "weights.2.pkl" in out_files, "Missing output model weights!"
