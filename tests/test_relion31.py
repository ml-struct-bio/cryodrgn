"""Tests of compatibility with the RELION 3.1 format for .star files (optics groups)."""

import pytest
import argparse
import os
import pickle
from cryodrgn.starfile import Starfile
from cryodrgn.commands import downsample, parse_ctf_star, parse_pose_star
from cryodrgn.commands_utils import filter_star


# TODO -- convert these starfiles to 3.0 and compare outputs across these tests
@pytest.fixture
def relion_starfile(request):
    return os.path.join(pytest.DATADIR, request.param)


@pytest.mark.xfail(reason="coming soon")
@pytest.mark.parametrize(
    "relion_starfile", ["relion31.star", "relion31.v2.star"], indirect=True
)
def test_downsample(tmpdir, relion_starfile):
    args = downsample.add_args(argparse.ArgumentParser()).parse_args(
        [
            f"{relion_starfile}",
            "-D",
            "32",
            "-o",
            os.path.join(tmpdir, "temp.mrcs"),
        ]
    )
    downsample.main(args)


@pytest.mark.parametrize(
    "relion_starfile, indices",
    [
        ("relion31.star", "just-4"),
        pytest.param(
            "relion31.v2.star",
            "just-4",
            marks=pytest.mark.xfail(
                raises=AssertionError,
                reason="don't yet support relion31 files with optics table last!",
            ),
        ),
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
    stardata = Starfile.load(starfile)
    with open(indices.path, "rb") as f:
        ind = pickle.load(f)

    assert stardata.relion31
    assert stardata.df.shape[0] == len(ind)
    assert stardata.data_optics.df.shape[0] == 1


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
    stardata = Starfile.load(relion_starfile)
    with open(pose_file, "rb") as f:
        poses = pickle.load(f)

    assert len(poses) == 2
    assert poses[0].shape[0] == stardata.df.shape[0]
    assert poses[1].shape[0] == stardata.df.shape[0]


@pytest.mark.xfail(reason="coming soon")
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
