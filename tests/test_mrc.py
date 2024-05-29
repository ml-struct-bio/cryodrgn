import pytest
import os
import numpy as np
import torch
from cryodrgn.source import ImageSource


@pytest.fixture
def mrcs_data():
    return ImageSource.from_file(
        os.path.join(pytest.DATADIR, "toy_projections.mrcs"), lazy=False
    ).images()


def test_lazy_loading(mrcs_data):
    lazy_data = ImageSource.from_file(
        os.path.join(pytest.DATADIR, "toy_projections.mrcs"), lazy=True
    ).images()

    assert isinstance(mrcs_data, torch.Tensor)
    assert isinstance(lazy_data, torch.Tensor)
    assert torch.allclose(mrcs_data, lazy_data)

    lazy_np = np.array(lazy_data)
    busy_np = np.array(mrcs_data[:])
    assert (lazy_np == busy_np).all()


def test_star(mrcs_data):
    star_data = ImageSource.from_file(
        os.path.join(pytest.DATADIR, "toy_projections.star")
    ).images()

    assert isinstance(star_data, torch.Tensor)
    assert torch.allclose(star_data, mrcs_data)

    star_np = np.array(star_data)
    busy_np = np.array(mrcs_data[:])
    assert (star_np == busy_np).all()


def test_txt(mrcs_data):
    txt_data = ImageSource.from_file(
        os.path.join(pytest.DATADIR, "toy_projections.txt")
    ).images()

    assert isinstance(txt_data, torch.Tensor)
    assert torch.allclose(txt_data, mrcs_data)

    txt_np = np.array(txt_data)
    busy_np = np.array(mrcs_data[:])
    assert (txt_np == busy_np).all()
