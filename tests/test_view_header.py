import pytest
import argparse
from cryodrgn.commands_utils import view_header


@pytest.mark.parametrize("particles", ["hand", "toy.mrcs"], indirect=True)
def test_view_header(particles):
    args = view_header.add_args(argparse.ArgumentParser()).parse_args([particles.path])
    view_header.main(args)
