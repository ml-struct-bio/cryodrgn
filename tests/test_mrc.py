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


@pytest.mark.integration
@pytest.mark.slow
def test_lazy_loading(mrcs_data):
    lazy_data = ImageSource.from_file(
        os.path.join(pytest.DATADIR, "toy_projections.mrcs"), lazy=True
    ).images()

    assert isinstance(mrcs_data, torch.Tensor)
    assert isinstance(lazy_data, torch.Tensor)
    assert torch.allclose(mrcs_data, lazy_data)

    # Use explicit numpy conversion methods to avoid numpy 2.0 warnings
    lazy_np = lazy_data.numpy()
    busy_np = mrcs_data.numpy()
    assert (lazy_np == busy_np).all()


@pytest.mark.integration
@pytest.mark.slow
def test_star(mrcs_data):
    star_data = ImageSource.from_file(
        os.path.join(pytest.DATADIR, "toy_projections.star")
    ).images()

    assert isinstance(star_data, torch.Tensor)
    assert torch.allclose(star_data, mrcs_data)

    # Use explicit numpy conversion methods to avoid numpy 2.0 warnings
    star_np = star_data.numpy()
    busy_np = mrcs_data.numpy()
    assert (star_np == busy_np).all()


@pytest.mark.integration
@pytest.mark.slow
def test_txt(mrcs_data):
    txt_data = ImageSource.from_file(
        os.path.join(pytest.DATADIR, "toy_projections.txt")
    ).images()

    assert isinstance(txt_data, torch.Tensor)
    assert torch.allclose(txt_data, mrcs_data)

    # Use explicit numpy conversion methods to avoid numpy 2.0 warnings
    txt_np = txt_data.numpy()
    busy_np = mrcs_data.numpy()
    assert (txt_np == busy_np).all()
