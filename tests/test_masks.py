"""Unit tests of the cryodrgn_utils gen_mask command."""

import pytest
import os
from hashlib import md5
from cryodrgn.utils import run_command
from cryodrgn.mrcfile import parse_mrc


def hash_file(filename: str) -> str:
    assert os.path.exists(filename)

    with open(filename, "rb") as f:
        file_hash = md5()
        while chunk := f.read(8192):
            file_hash.update(chunk)

    return file_hash.hexdigest()


@pytest.mark.parametrize("volume", ["toy-small", "spike"], indirect=True)
@pytest.mark.parametrize("dist", [2, 5])
@pytest.mark.parametrize("dilate", [3, 7])
@pytest.mark.parametrize("apix", [None, 1, 2.79])
def test_mask_fidelity(tmpdir, volume, dist, dilate, apix) -> None:
    """Test that we can compare two volumes produced during reconstruction training."""
    mask_file = os.path.join(tmpdir, f"{volume.label}_mask.mrc")
    cmd_str = f"cryodrgn_utils gen_mask {volume.path} {mask_file} "
    cmd_str += f"--dist {dist} --dilate {dilate} "
    if apix is not None:
        cmd_str += f"--Apix {apix}"

    out0, err = run_command(cmd_str)
    assert err == ""

    thresh_vals = {"toy-small": 0.5, "spike": 12.472}
    assert (
        round(float(out0.split("\n")[0].split("Threshold=")[1]), 4)
        == thresh_vals[volume.label]
    )

    assert os.path.exists(mask_file)
    vol_arr, vol_header = parse_mrc(volume.path)
    mask_arr, mask_header = parse_mrc(mask_file)

    assert mask_arr.shape == vol_arr.shape
    if apix is None:
        assert mask_header.apix == vol_header.apix
    else:
        assert mask_header.apix == apix

    mask_sums = {
        "toy-small": {
            2: {
                3: {None: 1760.1, 1: 1760.1, 2.79: 242.0},
                7: {None: 5571.7, 1: 5571.7, 2.79: 592.0},
            },
            5: {
                3: {None: 4256.2, 1: 4256.2, 2.79: 543.9},
                7: {None: 9545.2, 1: 9545.2, 2.79: 1021.3},
            },
        },
        "spike": {
            2: {
                3: {None: 0.0, 1: 3042.6, 2.79: 450.0},
                7: {None: 0.0, 1: 8542.9, 2.79: 1069.0},
            },
            5: {
                3: {None: 0.0, 1: 6853.3, 2.79: 1011.0},
                7: {None: 0.0, 1: 11745.8, 2.79: 1821.7},
            },
        },
    }

    # the mask is always at max value for the volume itself
    assert not ((mask_arr != 1) & (vol_arr > thresh_vals[volume.label])).any()
    out_sum = mask_arr[vol_arr <= thresh_vals[volume.label]].sum()
    assert round(float(out_sum), 1) == mask_sums[volume.label][dist][dilate][apix]


@pytest.mark.parametrize("volume", ["toy-small"], indirect=True)
@pytest.mark.parametrize("dist_val", [3, 5])
def test_png_output_file(tmpdir, volume, dist_val) -> None:
    mask_file = os.path.join(tmpdir, f"{volume.label}_{dist_val}_mask.mrc")
    plot_file = os.path.join(tmpdir, f"{volume.label}_{dist_val}_slices.png")

    out0, err = run_command(
        f"cryodrgn_utils gen_mask {volume.path} {mask_file} "
        f"-p {plot_file} --dist {dist_val}"
    )
    assert err == ""

    thresh_vals = {"toy-small": 0.5, "spike": 12.472, "50S-vol": 1.2661}
    assert (
        round(float(out0.split("\n")[0].split("Threshold=")[1]), 4)
        == thresh_vals[volume.label]
    )

    mask_hashes = {
        3: {
            "toy-small": "84b23b71ef048218874f9b468cee6abf",
        },
        5: {
            "toy-small": "84b9810568cc8c2d00b320d2dc24564e",
        },
    }
    assert hash_file(mask_file) == mask_hashes[dist_val][volume.label]

    plot_hashes = {
        3: {
            "toy-small": "88d713543bd093d3cba30a6b27ef0a34",
        },
        5: {
            "toy-small": "f112135a63307a8baadba51c6f97150e",
        },
    }
    assert hash_file(plot_file) == plot_hashes[dist_val][volume.label]
