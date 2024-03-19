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
        round(float(out0.split("\n")[0].split("Threshold=")[1]), 4)
        == thresh_vals[vol_idx]
    )

    mask_hashes = {
        1: "8c489a79e9bf4adddc9b8f47508f10e7",
        2: "310ed7776301d3785aea1fe88bd7daa0",
    }
    assert hash_file(mask_file) == mask_hashes[vol_idx]


@pytest.mark.parametrize(
    "train_dir",
    [{"train_cmd": "train_nn", "epochs": 3, "seed": 5555}],
    indirect=True,
)
@pytest.mark.parametrize("dist_val", [10, 20])
def test_png_output_file(trained_dir, dist_val) -> None:
    vol_file = os.path.join(trained_dir.outdir, "reconstruct.2.mrc")
    mask_file = os.path.join(trained_dir.outdir, "mask.mrc")
    plot_file = os.path.join(trained_dir.outdir, "slices.png")

    out0, err = run_command(
        f"cryodrgn_utils gen_mask {vol_file} {mask_file} "
        f"-p {plot_file} --dist {dist_val}"
    )
    assert err == ""
    assert round(float(out0.split("\n")[0].split("Threshold=")[1]), 4) == 0.477

    mask_hashes = {
        10: "1d970ec46645a4b9953d4f1bc0c2dfe9",
        20: "d880019cae20e440b257d77aa331aaa1",
    }
    assert hash_file(mask_file) == mask_hashes[dist_val]

    plot_hashes = {
        10: "71e898a77ce2913cc4d755e00c8bfd68",
        20: "6f6a6ce284134fa43478a220d271f5f2",
    }
    assert hash_file(plot_file) == plot_hashes[dist_val]
