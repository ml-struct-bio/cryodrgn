import argparse
import os
import os.path
import pytest
from cryodrgn.commands import downsample, parse_ctf_star, parse_pose_star

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


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
            "--relion31",
            "-D",
            "32",
            "-o",
            os.path.join(tmpdir, "temp.mrcs"),
        ]
    )
    downsample.main(args)


@pytest.mark.xfail(reason="coming soon")
@pytest.mark.parametrize(
    "relion_starfile", ["relion31.star", "relion31.v2.star"], indirect=True
)
def test_parse_pose_star(tmpdir, relion_starfile):
    args = parse_pose_star.add_args(argparse.ArgumentParser()).parse_args(
        [
            f"{relion_starfile}",
            "--relion31",
            "-D",
            "256",
            "--Apix",
            "1",
            "-o",
            os.path.join(tmpdir, "pose1.pkl"),
        ]
    )
    parse_pose_star.main(args)


@pytest.mark.parametrize(
    "relion_starfile", ["relion31.star", "relion31.v2.star"], indirect=True
)
def test_parse_pose_star_ignore_optics(tmpdir, relion_starfile):
    args = [
        f"{relion_starfile}",
        "-D",
        "256",
        "--Apix",
        "1",
        "-o",
        os.path.join(tmpdir, "pose.pkl"),
        "--ignore-optics",
    ]
    parser = argparse.ArgumentParser()
    parse_pose_star.add_args(parser)
    parse_pose_star.main(parser.parse_args(args))


@pytest.mark.xfail(reason="coming soon")
@pytest.mark.parametrize(
    "relion_starfile", ["relion31.star", "relion31.v2.star"], indirect=True
)
def test_parse_ctf_star(tmpdir, relion_starfile):
    args = parse_ctf_star.add_args(argparse.ArgumentParser()).parse_args(
        [
            f"{relion_starfile}",
            "--relion31",
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
