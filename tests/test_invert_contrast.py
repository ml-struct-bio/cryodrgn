import pytest
import os
import shutil
import numpy as np
import torch
from cryodrgn.source import ImageSource
from cryodrgn.mrcfile import parse_mrc
from cryodrgn.utils import run_command


@pytest.mark.parametrize("volume", ["hand", "toy"], indirect=True)
def test_output(tmpdir, volume):
    """Try different ways of specifying the output file."""
    vol_file = os.path.join(tmpdir, "vol.mrc")
    inv_file = os.path.join(
        tmpdir, os.path.basename(vol_file).replace(".mrc", "_inverted.mrc")
    )
    shutil.copyfile(volume.path, vol_file)

    out, err = run_command(f"cryodrgn_utils invert_contrast {vol_file}")
    assert err == ""

    mrcs_data, _ = parse_mrc(vol_file)
    flipped_data, _ = parse_mrc(inv_file)
    # torch doesn't let us use a -ve stride, hence the conv
    assert np.allclose(flipped_data, -mrcs_data)

    inv_file = os.path.join(
        tmpdir, "vols", os.path.basename(vol_file).replace(".mrc", "_inverted.mrc")
    )

    out, err = run_command(f"cryodrgn_utils invert_contrast {vol_file} -o {inv_file}")
    assert err == ""

    mrcs_data, _ = parse_mrc(volume.path)
    flipped_data, _ = parse_mrc(inv_file)
    assert np.allclose(flipped_data, -mrcs_data)


@pytest.mark.parametrize("volume", ["hand"], indirect=True)
def test_mrc_file(tmpdir, volume):
    inv_file = os.path.join(
        tmpdir, os.path.basename(volume.path).replace(".mrc", "_inverted.mrc")
    )

    out, err = run_command(
        f"cryodrgn_utils invert_contrast {volume.path} -o {inv_file}"
    )
    assert err == ""

    mrcs_data, _ = parse_mrc(volume.path)
    flipped_data, _ = parse_mrc(inv_file)
    assert np.allclose(flipped_data, -mrcs_data)


@pytest.mark.parametrize("volume", ["hand"], indirect=True)
def test_image_source(tmpdir, volume):
    inv_file = os.path.join(
        tmpdir, os.path.basename(volume.path).replace(".mrc", "_inverted.mrc")
    )

    out, err = run_command(
        f"cryodrgn_utils invert_contrast {volume.path} -o {inv_file}"
    )
    assert err == ""

    mrcs_data = ImageSource.from_file(volume.path).images()
    flipped_data = ImageSource.from_file(inv_file).images()
    assert torch.allclose(flipped_data, -mrcs_data)
