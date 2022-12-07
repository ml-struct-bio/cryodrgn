import argparse
import os
import os.path
import pytest
from cryodrgn.commands_utils import write_star

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture
def particles_starfile():
    return f"{DATA_FOLDER}/FinalRefinement-OriginalParticles-PfCRT.star"


def test_writestarfile1():
    os.makedirs("output", exist_ok=True)
    args = write_star.add_args(argparse.ArgumentParser()).parse_args(
        [
            f"{DATA_FOLDER}/hand.5.mrcs",
            "--ctf",
            f"{DATA_FOLDER}/ctf1.pkl",
            "-o",
            "output/star1.star",
        ]
    )
    write_star.main(args)


def test_writestarfile2(particles_starfile):
    os.makedirs("output", exist_ok=True)
    args = write_star.add_args(argparse.ArgumentParser()).parse_args(
        [
            f"{DATA_FOLDER}/hand.5.mrcs",
            "--ctf",
            f"{DATA_FOLDER}/ctf1.pkl",
            "-o",
            "output/star2.star",
        ]
    )
    write_star.main(args)
