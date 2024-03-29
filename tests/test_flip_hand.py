import os
import shutil
import numpy as np
from cryodrgn.source import ImageSource
from cryodrgn.mrc import MRCFile
from cryodrgn.utils import run_command

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


def test_output():
    """Try different ways of specifying the output file."""
    os.makedirs("output", exist_ok=True)

    vol_file = os.path.join("output", "hand-vol.mrc")
    shutil.copyfile(os.path.join(DATA_FOLDER, "hand-vol.mrc"), vol_file)
    flipped_file = os.path.join(
        "output", os.path.basename(vol_file).replace(".mrc", "_flipped.mrc")
    )

    out, err = run_command(f"cryodrgn_utils flip_hand {vol_file}")
    assert err == ""

    mrcs_data, _ = MRCFile.parse(vol_file)
    flipped_data, _ = MRCFile.parse(flipped_file)
    assert np.allclose(flipped_data, mrcs_data[::-1])

    flipped_file = os.path.join(
        "output", "vols", os.path.basename(vol_file).replace(".mrc", "_flipped.mrc")
    )

    out, err = run_command(f"cryodrgn_utils flip_hand {vol_file} -o {flipped_file}")
    assert err == ""

    mrcs_data, _ = MRCFile.parse(vol_file)
    flipped_data, _ = MRCFile.parse(flipped_file)
    assert np.allclose(flipped_data, mrcs_data[::-1])


def test_mrc_file():
    vol_file = os.path.join(DATA_FOLDER, "hand-vol.mrc")
    flipped_file = os.path.join(
        "output", os.path.basename(vol_file).replace(".mrc", "_flipped.mrc")
    )

    out, err = run_command(f"cryodrgn_utils flip_hand {vol_file} -o {flipped_file}")
    assert err == ""

    mrcs_data, _ = MRCFile.parse(vol_file)
    flipped_data, _ = MRCFile.parse(flipped_file)
    assert np.allclose(flipped_data, mrcs_data[::-1])


def test_image_source():
    vol_file = os.path.join(DATA_FOLDER, "toy_projections.mrc")
    flipped_file = os.path.join(
        "output", os.path.basename(vol_file).replace(".mrc", "_flipped.mrc")
    )

    out, err = run_command(f"cryodrgn_utils flip_hand {vol_file} -o {flipped_file}")
    assert err == ""

    mrcs_data = ImageSource.from_file(vol_file).images()
    flipped_data = ImageSource.from_file(flipped_file).images()
    # torch doesn't let us use a -ve stride, hence the conversion below
    assert np.allclose(np.array(flipped_data.cpu()), np.array(mrcs_data.cpu())[::-1])
