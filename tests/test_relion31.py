"""Tests of compatibility with the RELION 3.1 format for .star files (optics groups)."""

import pytest
import argparse
import os
import pickle
import numpy as np
from cryodrgn.commands import downsample, parse_ctf_star, parse_pose_star
from cryodrgn.commands_utils import filter_star
from cryodrgn.source import ImageSource, parse_star


# TODO -- convert these starfiles to 3.0 and compare outputs across these tests
@pytest.fixture
def relion_starfile(request):
    return os.path.join(pytest.DATADIR, request.param)


@pytest.mark.parametrize(
    "relion_starfile", ["relion31.star", "relion31.v2.star"], indirect=True
)
@pytest.mark.parametrize("resolution", [32, 64])
def test_downsample(tmpdir, relion_starfile, resolution):
    input_data = ImageSource.from_file(relion_starfile, lazy=False)
    out_mrcs = os.path.join(tmpdir, "temp.mrcs")

    parser = argparse.ArgumentParser()
    downsample.add_args(parser)
    args = parser.parse_args(
        [f"{relion_starfile}", "-D", str(resolution), "-o", out_mrcs]
    )
    downsample.main(args)

    output_data = ImageSource.from_file(out_mrcs, lazy=False)
    in_imgs = input_data.images()
    out_imgs = output_data.images()
    assert out_imgs.shape == (in_imgs.shape[0], resolution, resolution)
    assert np.isclose(in_imgs.sum(), out_imgs.sum())


@pytest.mark.parametrize(
    "relion_starfile, indices",
    [
        ("relion31.star", "just-4"),
        ("relion31.v2.star", "just-4"),
        ("relion31.6opticsgroups.star", "just-4"),
        ("relion31.6opticsgroups.star", "just-5"),
    ],
    indirect=True,
)
def test_filter_star(tmpdir, relion_starfile, indices):
    parser = argparse.ArgumentParser()
    filter_star.add_args(parser)
    starfile = os.path.join(
        tmpdir, f"filtered_{os.path.basename(relion_starfile)}_{indices.label}.star"
    )
    args = [
        f"{relion_starfile}",
        "-o",
        starfile,
        "--ind",
        indices.path,
    ]

    filter_star.main(parser.parse_args(args))
    stardata, data_optics = parse_star(starfile)
    with open(indices.path, "rb") as f:
        ind = pickle.load(f)

    assert data_optics is not None
    assert stardata.shape[0] == len(ind)
    assert data_optics.shape[0] == 1


@pytest.mark.parametrize(
    "relion_starfile, resolution, apix",
    [
        ("relion31.star", None, None),
        ("relion31.6opticsgroups.star", None, None),
        ("relion31.6opticsgroups.star", "256", "1.0"),
        ("relion31.star", "256", "1.0"),
        ("relion31.v2.star", "256", "1.0"),
    ],
    indirect=["relion_starfile"],
)
def test_parse_pose_star(tmpdir, relion_starfile, resolution, apix):
    parser = argparse.ArgumentParser()
    parse_pose_star.add_args(parser)
    pose_file = os.path.join(
        tmpdir, f"pose_{os.path.basename(relion_starfile)}_{resolution}_{apix}.pkl"
    )
    args = [f"{relion_starfile}", "-o", pose_file]
    if resolution is not None:
        args += ["-D", resolution]
    if apix is not None:
        args += ["--Apix", apix]

    parse_pose_star.main(parser.parse_args(args))
    stardata, data_optics = parse_star(relion_starfile)
    with open(pose_file, "rb") as f:
        poses = pickle.load(f)

    assert len(poses) == 2
    assert poses[0].shape[0] == stardata.shape[0]
    assert poses[1].shape[0] == stardata.shape[0]


@pytest.mark.parametrize(
    "relion_starfile", ["relion31.star", "relion31.v2.star"], indirect=True
)
def test_parse_ctf_star(tmpdir, relion_starfile):
    parser = argparse.ArgumentParser()
    parse_ctf_star.add_args(parser)
    args = parser.parse_args(
        [
            f"{relion_starfile}",
            "-D",
            "256",
            "--Apix",
            "1",
            "--kv",
            "300",
            "-w",
            ".1",
            "--ps",
            "0",
            "--cs",
            "2.7",
            "-o",
            os.path.join(tmpdir, "ctf.pkl"),
        ]
    )
    parse_ctf_star.main(args)
