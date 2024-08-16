"""Tests of compatibility with the RELION 3.1 format for .star files (optics groups)."""

import pytest
import argparse
import os
import pickle
from cryodrgn.commands import parse_ctf_star, parse_pose_star
from cryodrgn.commands_utils import filter_star
from cryodrgn.source import parse_star


# TODO -- convert these starfiles to 3.0 and compare outputs across these tests
@pytest.fixture
def relion_starfile(request):
    return os.path.join(pytest.DATADIR, request.param)


pytestmark = pytest.mark.parametrize(
    "relion_starfile",
    ["relion31.star", "relion31.v2.star", "relion31.6opticsgroups.star"],
    indirect=True,
)


@pytest.mark.parametrize("indices", ["just-4"], indirect=True)
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


@pytest.mark.parametrize("resolution, apix", [(None, None), ("256", "1.0")])
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

    with open(os.path.join(tmpdir, "ctf.pkl"), "rb") as f:
        ctf_params = pickle.load(f)

    print(ctf_params.shape)
