import os.path
import torch
import pytest
from cryodrgn.source import ImageSource

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture
def mrcs_data():
    return ImageSource.from_file(
        f"{DATA_FOLDER}/toy_projections.mrcs", lazy=False
    ).images()


def test_lazy_loading(mrcs_data):
    data = ImageSource.from_file(
        f"{DATA_FOLDER}/toy_projections.mrcs", lazy=True
    ).images()
    assert isinstance(data, torch.Tensor)
    assert torch.allclose(data, mrcs_data)


def test_star(mrcs_data):
    data = ImageSource.from_file(f"{DATA_FOLDER}/toy_projections.star").images()
    assert isinstance(data, torch.Tensor)
    assert torch.allclose(data, mrcs_data)


def test_txt(mrcs_data):
    data = ImageSource.from_file(f"{DATA_FOLDER}/toy_projections.txt").images()
    assert isinstance(data, torch.Tensor)
    assert torch.allclose(data, mrcs_data)
