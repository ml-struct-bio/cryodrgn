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


@pytest.fixture
def volume(request) -> tuple[str, str]:
    return (
        os.path.join(
            os.path.dirname(__file__),
            "..",
            "testing",
            "data",
            ".".join([request.param, "mrc"]),
        ),
        request.param,
    )


@pytest.mark.parametrize(
    "volume", ["toymodel_small_nocenter", "spike-vol", "50S-vol"], indirect=True
)
def test_mask_fidelity(volume) -> None:
    """Test that we can compare two volumes produced during reconstruction training."""
    os.makedirs("output", exist_ok=True)
    vol_file, vol_lbl = volume
    mask_file = os.path.join("output", f"{vol_lbl}_mask.mrc")
    out0, err = run_command(
        "cryodrgn_utils gen_mask " f"{vol_file} {mask_file} --dist 2 --dilate 3"
    )
    assert err == ""

    thresh_vals = {
        "toymodel_small_nocenter": 0.5,
        "spike-vol": 12.472,
        "50S-vol": 1.2661,
    }
    assert (
        round(float(out0.split("\n")[0].split("Threshold=")[1]), 4)
        == thresh_vals[vol_lbl]
    )

    mask_hashes = {
        "toymodel_small_nocenter": "85f8073b2a933f7d3fb0f890d8630af8",
        "spike-vol": "e0dc53ecc7566085daee70402672426f",
        "50S-vol": "a1ab19e933b6bee5cf1119a01f2cb3de",
    }
    assert hash_file(mask_file) == mask_hashes[vol_lbl]


@pytest.mark.parametrize("volume", ["toymodel_small_nocenter"], indirect=True)
@pytest.mark.parametrize("dist_val", [3, 5])
def test_png_output_file(volume, dist_val) -> None:
    os.makedirs("output", exist_ok=True)
    vol_file, vol_lbl = volume
    mask_file = os.path.join("output", f"{vol_lbl}_{dist_val}_mask.mrc")
    plot_file = os.path.join("output", f"{vol_lbl}_{dist_val}slices.png")

    out0, err = run_command(
        f"cryodrgn_utils gen_mask {vol_file} {mask_file} "
        f"-p {plot_file} --dist {dist_val}"
    )
    assert err == ""

    thresh_vals = {
        "toymodel_small_nocenter": 0.5,
        "spike-vol": 12.472,
        "50S-vol": 1.2661,
    }
    assert (
        round(float(out0.split("\n")[0].split("Threshold=")[1]), 4)
        == thresh_vals[vol_lbl]
    )

    mask_hashes = {
        3: {
            "toymodel_small_nocenter": "eafaaafd35bdbbdc880802367f892921",
        },
        5: {
            "toymodel_small_nocenter": "3ddb1ca57d656d9b8bbc2cf2b045c3b8",
        },
    }
    assert hash_file(mask_file) == mask_hashes[dist_val][vol_lbl]

    plot_hashes = {
        3: {
            "toymodel_small_nocenter": "88d713543bd093d3cba30a6b27ef0a34",
        },
        5: {
            "toymodel_small_nocenter": "f112135a63307a8baadba51c6f97150e",
        },
    }
    assert hash_file(plot_file) == plot_hashes[dist_val][vol_lbl]
