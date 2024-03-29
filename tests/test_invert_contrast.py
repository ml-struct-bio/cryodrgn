import os
import shutil
import numpy as np
import torch
from cryodrgn.source import ImageSource
from cryodrgn.mrc import MRCFile
from cryodrgn.utils import run_command

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


def test_output():
    """Try different ways of specifying the output file."""
    os.makedirs("output", exist_ok=True)

    vol_file = os.path.join("output", "hand-vol.mrc")
    shutil.copyfile(os.path.join(DATA_FOLDER, "hand-vol.mrc"), vol_file)
    inv_file = os.path.join(
        "output", os.path.basename(vol_file).replace(".mrc", "_inverted.mrc")
    )

    out, err = run_command(f"cryodrgn_utils invert_contrast {vol_file}")
    assert err == ""

    mrcs_data, _ = MRCFile.parse(vol_file)
    flipped_data, _ = MRCFile.parse(inv_file)
    # torch doesn't let us use a -ve stride, hence the conv
    assert np.allclose(flipped_data, -mrcs_data)

    inv_file = os.path.join(
        "output", "vols", os.path.basename(vol_file).replace(".mrc", "_inverted.mrc")
    )

    out, err = run_command(f"cryodrgn_utils invert_contrast {vol_file} -o {inv_file}")
    assert err == ""

    mrcs_data, _ = MRCFile.parse(vol_file)
    flipped_data, _ = MRCFile.parse(inv_file)
    assert np.allclose(flipped_data, -mrcs_data)


def test_mrc_file():
    vol_file = os.path.join(DATA_FOLDER, "hand-vol.mrc")
    inv_file = os.path.join(
        "output", os.path.basename(vol_file).replace(".mrc", "_inverted.mrc")
    )

    out, err = run_command(f"cryodrgn_utils invert_contrast {vol_file} -o {inv_file}")
    assert err == ""

    mrcs_data, _ = MRCFile.parse(vol_file)
    flipped_data, _ = MRCFile.parse(inv_file)
    assert np.allclose(flipped_data, -mrcs_data)


def test_image_source():
    vol_file = os.path.join(DATA_FOLDER, "toy_projections.mrc")
    inv_file = os.path.join(
        "output", os.path.basename(vol_file).replace(".mrc", "_inverted.mrc")
    )

    out, err = run_command(f"cryodrgn_utils invert_contrast {vol_file} -o {inv_file}")
    assert err == ""

    mrcs_data = ImageSource.from_file(vol_file).images()
    flipped_data = ImageSource.from_file(inv_file).images()
    assert torch.allclose(flipped_data, -mrcs_data)
