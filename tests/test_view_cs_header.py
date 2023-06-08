import os.path
import argparse
import pytest
from cryodrgn.commands_utils import view_cs_header

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture
def particles_csfile():
    return f"{DATA_FOLDER}/cryosparc_P12_J24_001_particles.cs"


def test_view_cs_header():
    args = view_cs_header.add_args(argparse.ArgumentParser()).parse_args(
        [f"{DATA_FOLDER}/cryosparc_P12_J24_001_particles.cs"]
    )
    view_cs_header.main(args)
