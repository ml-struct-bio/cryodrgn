import pytest
import os
import argparse
from cryodrgn.commands_utils import phase_flip


def test_phase_flip(tmpdir):
    args = phase_flip.add_args(argparse.ArgumentParser()).parse_args(
        [
            os.path.join(pytest.DATADIR, "relion31.mrcs"),
            os.path.join(pytest.DATADIR, "ctf1.pkl"),
            "-o",
            os.path.join(tmpdir, "phase_flipped.mrcs"),
        ]
    )
    phase_flip.main(args)
