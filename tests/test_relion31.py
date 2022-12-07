import argparse
import os
import os.path
import pytest
from cryodrgn.commands import downsample, parse_ctf_star, parse_pose_star

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.mark.xfail(reason="coming soon")
@pytest.mark.parametrize(
    "relion_starfile",
    [f"{DATA_FOLDER}/relion31.star", f"{DATA_FOLDER}/relion31.v2.star"],
)
def test_downsample(relion_starfile):
    os.makedirs("output", exist_ok=True)
    args = downsample.add_args(argparse.ArgumentParser()).parse_args(
        [
            f"{relion_starfile}",
            "--relion31",
            "-D",
            "32",
            "-o",
            "output/temp.mrcs",
        ]
    )
    downsample.main(args)


@pytest.mark.xfail(reason="coming soon")
@pytest.mark.parametrize(
    "relion_starfile",
    [f"{DATA_FOLDER}/relion31.star", f"{DATA_FOLDER}/relion31.v2.star"],
)
def test_parse_pose_star(relion_starfile):
    os.makedirs("output", exist_ok=True)
    args = parse_pose_star.add_args(argparse.ArgumentParser()).parse_args(
        [
            f"{relion_starfile}",
            "--relion31",
            "-D",
            "256",
            "--Apix",
            "1",
            "-o",
            "output/pose.pkl",
        ]
    )
    parse_pose_star.main(args)


@pytest.mark.xfail(reason="coming soon")
@pytest.mark.parametrize(
    "relion_starfile",
    [f"{DATA_FOLDER}/relion31.star", f"{DATA_FOLDER}/relion31.v2.star"],
)
def test_parse_ctf_star(relion_starfile):
    os.makedirs("output", exist_ok=True)
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
            "output/ctf.pkl",
        ]
    )
    parse_ctf_star.main(args)
