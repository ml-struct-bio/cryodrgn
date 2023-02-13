import os.path
import torch
import pytest
from cryodrgn.source import ImageSource

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture
def mrcs_data():
    return ImageSource.from_mrcs(
        f"{DATA_FOLDER}/toy_projections.mrcs", lazy=False
    ).images()


def test_loading_mrcs(mrcs_data):
    src = ImageSource.from_file(f"{DATA_FOLDER}/toy_projections.mrcs")
    # Indexing inside the 'images' attributes causes the hitherto lazy data to turn into a Tensor
    assert torch.allclose(src[:], mrcs_data)
    # We can, of course, do selective indexing to avoid realizing ALL the underlying data to memory
    assert torch.allclose(src[10:15], mrcs_data[10:15])


def test_loading_starfile(mrcs_data):
    src = ImageSource.from_file(
        f"{DATA_FOLDER}/toy_projections_13.star", datadir=DATA_FOLDER
    )
    arr = src[7:12]
    assert arr.shape == (5, 30, 30)


def test_loading_txtfile(mrcs_data):
    src = ImageSource.from_file(f"{DATA_FOLDER}/toy_projections_2.txt")
    # Each line of the txt file points to an .mrcs file with 1000 particles. Try to get a slice across the two.
    arr = src[990:1005]
    assert arr.shape == (15, 30, 30)
    assert torch.allclose(arr[:10], mrcs_data[-10:])
    assert torch.allclose(arr[10:], mrcs_data[:5])


def test_loading_csfile(mrcs_data):
    src = ImageSource.from_file(f"{DATA_FOLDER}/empiar_10076_7.cs")
    arr = src[:]
    assert arr.shape == (7, 320, 320)
    starfile_data = ImageSource.from_star(
        f"{DATA_FOLDER}/empiar_10076_7.star", datadir=DATA_FOLDER, lazy=False
    )[:]

    assert torch.allclose(arr, starfile_data)


def test_source_iteration():
    # An ImageSource can be iterated over, with an optional chunksize (default 1000
    src = ImageSource.from_file(f"{DATA_FOLDER}/toy_projections.mrcs")
    for indices, chunk in src.chunks(chunksize=300):
        assert isinstance(chunk, torch.Tensor)
        assert chunk.ndim == 3
        # chunk will be (300, L, L) in shape except the last iteration where it may be (<300, L, L)
        assert chunk.shape[0] <= 300
