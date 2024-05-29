import os
import argparse
from cryodrgn.commands_utils import phase_flip

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


def test_phase_flip(tmpdir):
    args = phase_flip.add_args(argparse.ArgumentParser()).parse_args(
        [
            f"{DATA_FOLDER}/relion31.mrcs",
            f"{DATA_FOLDER}/ctf1.pkl",
            "-o",
            os.path.join(tmpdir, "phase_flipped.mrcs"),
        ]
    )
    phase_flip.main(args)
