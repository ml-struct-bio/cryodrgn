"""Unit tests of the cryodrgn_utils gen_mask command."""

import pytest
import os
from hashlib import md5
from cryodrgn.utils import run_command


def hash_file(filename: str) -> str:
    assert os.path.exists(filename)

    with open(filename, "rb") as f:
        file_hash = md5()
        while chunk := f.read(8192):
            file_hash.update(chunk)

    return file_hash.hexdigest()


@pytest.mark.parametrize(
    "train_dir",
    [{"train_cmd": "train_nn", "epochs": 3, "seed": 987}],
    indirect=True,
)
@pytest.mark.parametrize("vol_idx", [1, 2])
def test_fidelity(trained_dir, vol_idx) -> None:
    """Test that we can compare two volumes produced during reconstruction training."""

    vol_file = os.path.join(trained_dir.outdir, f"reconstruct.{vol_idx}.mrc")
    mask_file = os.path.join(trained_dir.outdir, f"mask.{vol_idx}.mrc")
    out0, err = run_command(f"cryodrgn_utils gen_mask {vol_file} {mask_file}")
    assert err == ""

    thresh_vals = {1: 0.3228, 2: 0.3134}
    assert (
        round(float(out0.split("\n")[0].split("Threshold: ")[1]), 4)
        == thresh_vals[vol_idx]
    )

    mask_hashes = {
        1: "df67e4492a152fda1410eb52dfd3b967",
        2: "31a7cb816e694b9d8b307d8bbcfbe199",
    }
    assert hash_file(mask_file) == mask_hashes[vol_idx]


@pytest.mark.parametrize(
    "train_dir",
    [{"train_cmd": "train_nn", "epochs": 3, "seed": 5555}],
    indirect=True,
)
def test_png_output_file(trained_dir) -> None:
    vol_file = os.path.join(trained_dir.outdir, "reconstruct.2.mrc")
    mask_file = os.path.join(trained_dir.outdir, "mask.mrc")
    plot_file = os.path.join(trained_dir.outdir, "slices.png")

    out0, err = run_command(
        f"cryodrgn_utils gen_mask {vol_file} {mask_file} -p {plot_file}"
    )
    assert err == ""

    assert round(float(out0.split("\n")[0].split("Threshold: ")[1]), 4) == 0.477
    assert hash_file(mask_file) == "06b68d753304c5ae245db037e3c3e681"
    assert hash_file(plot_file) == "9e56fe7f013fa57619257a4c0b13f757"
