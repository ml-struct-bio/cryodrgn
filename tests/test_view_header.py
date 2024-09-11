import pytest
import argparse
from cryodrgn.commands_utils import view_header


@pytest.mark.parametrize("particles", ["hand", "toy.mrcs"], indirect=True)
def test_view_header(particles):
    parser = argparse.ArgumentParser()
    view_header.add_args(parser)
    view_header.main(parser.parse_args([particles.path]))
