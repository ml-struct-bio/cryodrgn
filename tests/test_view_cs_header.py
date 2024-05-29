import os.path
import argparse
import pytest
from cryodrgn.commands_utils import view_cs_header


def test_view_cs_header():
    args = view_cs_header.add_args(argparse.ArgumentParser()).parse_args(
        [os.path.join(pytest.DATADIR, "cryosparc_P12_J24_001_particles.cs")]
    )
    view_cs_header.main(args)
