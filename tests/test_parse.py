import argparse
import os
import os.path

import pytest

from cryodrgn.commands import (
    parse_ctf_csparc,
    parse_ctf_star,
    parse_pose_csparc,
    parse_pose_star,
)
from cryodrgn.utils import assert_pkl_close

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture
def particles_starfile():
    return f"{DATA_FOLDER}/FinalRefinement-OriginalParticles-PfCRT.star"


@pytest.fixture
def particles_csfile():
    return f"{DATA_FOLDER}/cryosparc_P12_J24_001_particles.cs"


def test_parse_ctf_star(particles_starfile):
    os.makedirs("output", exist_ok=True)
    args = parse_ctf_star.add_args(argparse.ArgumentParser()).parse_args(
        [
            particles_starfile,
            "-w",
            "0.1",
            "-D",
            "300",
            "--Apix",
            "1.035",
            "-o",
            "output/ctf1.pkl",
            "--png",
            "output/ctf1.png",
        ]
    )
    parse_ctf_star.main(args)
    assert_pkl_close("output/ctf1.pkl", f"{DATA_FOLDER}/ctf1.pkl")


def test_parse_ctf_cs(particles_csfile):
    os.makedirs("output", exist_ok=True)
    args = parse_ctf_csparc.add_args(argparse.ArgumentParser()).parse_args(
        [particles_csfile, "-o", "output/ctf2.pkl", "--png", "output/ctf2.png"]
    )
    parse_ctf_csparc.main(args)
    assert_pkl_close("output/ctf2.pkl", f"{DATA_FOLDER}/ctf2.pkl")


def test_parse_pose_star(particles_starfile):
    os.makedirs("output", exist_ok=True)
    args = parse_pose_star.add_args(argparse.ArgumentParser()).parse_args(
        [
            particles_starfile,
            "-D",
            "300",
            "-o",
            "output/pose1.pkl",
        ]
    )
    parse_pose_star.main(args)
    assert_pkl_close("output/pose1.pkl", f"{DATA_FOLDER}/pose.star.pkl")


def test_parse_pose_cs(particles_csfile):
    os.makedirs("output", exist_ok=True)
    args = parse_pose_csparc.add_args(argparse.ArgumentParser()).parse_args(
        [
            particles_csfile,
            "-D",
            "180",
            "-o",
            "output/pose2.pkl",
        ]
    )
    parse_pose_csparc.main(args)
    assert_pkl_close("output/pose2.pkl", f"{DATA_FOLDER}/pose.cs.pkl")
