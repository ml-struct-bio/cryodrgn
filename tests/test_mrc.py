import os.path

import numpy as np
import pytest

from cryodrgn import dataset, mrc

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture
def mrcs_data():
    return mrc.parse_mrc(f"{DATA_FOLDER}/toy_projections.mrcs", lazy=False)[0]


def test_lazy_loading(mrcs_data):
    data, _ = mrc.parse_mrc(f"{DATA_FOLDER}/toy_projections.mrcs", lazy=True)
    data = np.asarray([x.get() for x in data])
    assert np.allclose(data, mrcs_data)


def test_star(mrcs_data):
    data = dataset.load_particles(f"{DATA_FOLDER}/toy_projections.star")
    assert isinstance(data, np.ndarray)
    assert np.allclose(data, mrcs_data)


def test_txt(mrcs_data):
    data = dataset.load_particles(f"{DATA_FOLDER}/toy_projections.txt")
    assert isinstance(data, np.ndarray)
    assert np.allclose(data, mrcs_data)
