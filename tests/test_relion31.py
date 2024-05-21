import argparse
import os
import os.path
import pytest
from cryodrgn.commands import downsample, parse_ctf_star, parse_pose_star


@pytest.fixture
def relion_starfile(request):
    return os.path.join(pytest.data_dir, request.param)


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
    args = [f"{relion_starfile}", "-o", os.path.join(tmpdir, "pose.pkl")]
    if resolution is not None:
        args += ["-D", resolution]
    if apix is not None:
        args += ["--Apix", apix]

    parse_pose_star.main(parser.parse_args(args))


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
            "1",
            "--ps",
            "0",
            "--cs",
            "2.7",
            "-o",
            os.path.join(tmpdir, "ctf.pkl"),
        ]
    )
    parse_ctf_star.main(args)
