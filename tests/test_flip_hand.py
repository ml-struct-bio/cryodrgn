import pytest
import os
import shutil
import numpy as np
from cryodrgn.source import ImageSource
from cryodrgn.mrcfile import parse_mrc
from cryodrgn.utils import run_command


@pytest.mark.parametrize("volume", ["toy", "hand"], indirect=True)
def test_output(tmpdir, volume):
    """Try different ways of specifying the output file."""

    vol_file = os.path.join(tmpdir, "hand-vol.mrc")
    shutil.copyfile(volume.path, vol_file)
    flipped_file = os.path.join(
        tmpdir, os.path.basename(vol_file).replace(".mrc", "_flipped.mrc")
    )

    out, err = run_command(f"cryodrgn_utils flip_hand {vol_file}")
    assert err == ""

    mrcs_data, _ = parse_mrc(vol_file)
    flipped_data, _ = parse_mrc(flipped_file)
    assert np.allclose(flipped_data, mrcs_data[::-1])

    flipped_file = os.path.join(
        tmpdir, "vols", os.path.basename(vol_file).replace(".mrc", "_flipped.mrc")
    )

    out, err = run_command(f"cryodrgn_utils flip_hand {vol_file} -o {flipped_file}")
    assert err == ""

    mrcs_data, _ = parse_mrc(vol_file)
    flipped_data, _ = parse_mrc(flipped_file)
    assert np.allclose(flipped_data, mrcs_data[::-1])


@pytest.mark.parametrize("volume", ["hand"], indirect=True)
def test_mrc_file(tmpdir, volume):
    flipped_file = os.path.join(
        tmpdir, os.path.basename(volume.path).replace(".mrc", "_flipped.mrc")
    )

    out, err = run_command(f"cryodrgn_utils flip_hand {volume.path} -o {flipped_file}")
    assert err == ""

    mrcs_data, _ = parse_mrc(volume.path)
    flipped_data, _ = parse_mrc(flipped_file)
    assert np.allclose(flipped_data, mrcs_data[::-1])


@pytest.mark.parametrize("volume", ["toy"], indirect=True)
def test_image_source(tmpdir, volume):
    flipped_file = os.path.join(
        tmpdir, os.path.basename(volume.path).replace(".mrc", "_flipped.mrc")
    )

    out, err = run_command(f"cryodrgn_utils flip_hand {volume.path} -o {flipped_file}")
    assert err == ""

    mrcs_data = ImageSource.from_file(volume.path).images()
    flipped_data = ImageSource.from_file(flipped_file).images()
    # torch doesn't let us use a -ve stride, hence the conversion below
    assert np.allclose(np.array(flipped_data.cpu()), np.array(mrcs_data.cpu())[::-1])
