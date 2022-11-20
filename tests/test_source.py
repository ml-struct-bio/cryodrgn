import os.path

import numpy as np
import pytest
from cryodrgn.source import ImageSource
from cryodrgn import mrc, starfile

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture
def mrcs_data():
    return mrc.parse_mrc(f"{DATA_FOLDER}/toy_projections.mrcs", lazy=False)[0]


def test_loading_mrcs(mrcs_data):
    data = ImageSource.from_file(f"{DATA_FOLDER}/toy_projections.mrcs")
    # Indexing inside the 'images' attributes causes the hitherto lazy data to turn into an ndarray
    assert np.allclose(data[:], mrcs_data)
    # We can, of course, do selective indexing to avoid realizing ALL the underlying data to memory
    assert np.allclose(data[10:15], mrcs_data[10:15])


def test_loading_starfile(mrcs_data):
    data = ImageSource.from_file(
        f"{DATA_FOLDER}/toy_projections_13.star", datadir=DATA_FOLDER
    )
    arr = data[7:12]
    assert arr.shape == (5, 30, 30)
    legacy_value = starfile.Starfile.load(
        f"{DATA_FOLDER}/toy_projections_13.star"
    ).get_particles(datadir=DATA_FOLDER, lazy=False)[7:12]
    assert isinstance(legacy_value, np.ndarray)

    assert np.allclose(
        arr,
        legacy_value,
    )


def test_loading_txtfile(mrcs_data):
    data = ImageSource.from_file(f"{DATA_FOLDER}/toy_projections_2.txt")
    # Each line of the txt file points to an .mrcs file with 1000 particles. Try to get a slice across the two.
    arr = data[990:1005]
    assert arr.shape == (15, 30, 30)
    assert np.allclose(arr[:10], mrcs_data[-10:])
    assert np.allclose(arr[10:], mrcs_data[:5])


def test_loading_csfile(mrcs_data):
    data = ImageSource.from_file(f"{DATA_FOLDER}/empiar_10076_7.cs")
    arr = data[:]
    assert arr.shape == (7, 320, 320)
    legacy_value = starfile.Starfile.load(
        f"{DATA_FOLDER}/empiar_10076_7.star"
    ).get_particles(datadir=DATA_FOLDER, lazy=False)
    assert isinstance(legacy_value, np.ndarray)

    assert np.allclose(arr, legacy_value)
