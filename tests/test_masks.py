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
    mask_hashes = {
        "toy-small": {
            2: {
                3: {
                    None: "85f8073b2a933f7d3fb0f890d8630af8",
                    1: "85f8073b2a933f7d3fb0f890d8630af8",
                    2.79: "bb4207c6d0610923c368d193fffa961b",
                },
                7: {
                    None: "b073bfd7a9ae12a6dbf61fc822845524",
                    1: "b073bfd7a9ae12a6dbf61fc822845524",
                    2.79: "e57542e0f6ec7a9df2152df99fdeeef3",
                },
            },
            5: {
                3: {
                    None: "71b5eb59642a8ffd012d47d805c3de5c",
                    1: "71b5eb59642a8ffd012d47d805c3de5c",
                    2.79: "b5823d127f0c7171f550f8f6b65ffc55",
                },
                7: {
                    None: "1037b5733e3c499d932592ac8574496d",
                    1: "1037b5733e3c499d932592ac8574496d",
                    2.79: "fd48d91f9de4ea2d4f3e4ad00bbc3f3a",
                },
            },
        },
        "spike": {
            2: {
                3: {
                    None: "e0dc53ecc7566085daee70402672426f",
                    1: "1a4e2f540a3878e1a1711a5ace8e7f32",
                    2.79: "3a09e498807a44316a87609a1dd06b2e",
                },
                7: {
                    None: "e0dc53ecc7566085daee70402672426f",
                    1: "e907111fa79e20aa2eeef7ca41a57d1b",
                    2.79: "72d1c7606076633a4eb3fd27391d21a4",
                },
            },
            5: {
                3: {
                    None: "e0dc53ecc7566085daee70402672426f",
                    1: "acf89e2e2d1583ac037bab0f50ce772c",
                    2.79: "fcc9b4fc697a7584c9fa6ff2b952439d",
                },
                7: {
                    None: "e0dc53ecc7566085daee70402672426f",
                    1: "7fa0570251fefeb807c56223af2ab58c",
                    2.79: "a659ea05459ed429ce8c7fbd7274cf2a",
                },
            },
        },
    }
    assert hash_file(mask_file) == mask_hashes[volume.label][dist][dilate][apix]


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
            "toy-small": "eafaaafd35bdbbdc880802367f892921",
        },
        5: {
            "toy-small": "3ddb1ca57d656d9b8bbc2cf2b045c3b8",
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
