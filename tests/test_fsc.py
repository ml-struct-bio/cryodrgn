"""Unit tests of the cryodrgn fsc command."""

import pytest
import os
import numpy as np
from cryodrgn.utils import run_command


@pytest.mark.parametrize(
    "train_dir",
    [{"train_cmd": "train_nn", "epochs": 5, "seed": 1030}],
    indirect=True,
)
def test_fidelity(trained_dir) -> None:
    """Test that we can compare two volumes produced during reconstruction training."""

    vol_file1 = os.path.join(trained_dir.outdir, "reconstruct.0.mrc")
    vol_file2 = os.path.join(trained_dir.outdir, "reconstruct.4.mrc")
    fsc_file = os.path.join(trained_dir.outdir, "fsc.txt")

    out0, err = run_command(f"cryodrgn_utils fsc {vol_file1} {vol_file2}")

    assert err == "", err
    assert float(out0.split()[4]) == 0.0
    assert float(out0.split()[5]) == 1.0
    assert float(out0.split()[8]) == 0.031250
    assert float(out0.split()[9]) == 0.893017

    assert out0.split()[-8] == "0.5:"
    assert round(float(out0.split()[-7]), 6) == 2.064516
    assert out0.split()[-2] == "0.143:"
    assert round(float(out0.split()[-1]), 6) == 2.064516

    assert (
        round(
            sum(float(x) for x in out0.split() if x[:2] == "0." and x[-1].isnumeric()),
            3,
        )
        == 34.584
    )
    assert len(out0.split("\n")) == 36

    out, err = run_command(f"cryodrgn_utils fsc {vol_file1} {vol_file2} -o {fsc_file}")
    fsc_vals = np.fromfile(fsc_file, dtype=float, sep=" ")

    assert err == "", err
    assert out.split()[4] == "0.5:"
    assert round(float(out.split()[5]), 6) == 2.064516
    assert out.split()[10] == "0.143:"
    assert round(float(out.split()[11]), 6) == 2.064516
    assert len(out.split("\n")) == 3
    assert round(10 ** (fsc_vals @ np.tile([1, -1], len(fsc_vals) // 2)), 5) == 0.99311


@pytest.mark.parametrize(
    "train_dir",
    [
        {"dataset": "hand", "train_cmd": "train_nn", "epochs": 3, "seed": 2456},
        {"dataset": "toy", "train_cmd": "train_nn", "epochs": 3, "seed": 2456},
    ],
    indirect=True,
)
@pytest.mark.parametrize("epochs", [(0, 1), (0, 2), (1, 2)])
def test_output_file(trained_dir, epochs: tuple[int, int]) -> None:
    vol_file1 = os.path.join(trained_dir.outdir, f"reconstruct.{epochs[0]}.mrc")
    vol_file2 = os.path.join(trained_dir.outdir, f"reconstruct.{epochs[1]}.mrc")
    fsc_file = os.path.join(trained_dir.outdir, "fsc.txt")

    out0, err = run_command(f"cryodrgn_utils fsc {vol_file1} {vol_file2}")
    assert err == "", err
    assert out0.split()[-8] == "0.5:"
    assert out0.split()[-2] == "0.143:"

    out, err = run_command(f"cryodrgn_utils fsc {vol_file1} {vol_file2} -o {fsc_file}")
    assert err == "", err
    fsc_vals = np.fromfile(fsc_file, dtype=float, sep=" ")
    for i, fsc_val in enumerate(fsc_vals):
        assert round(fsc_val, 6) == round(float(out0.split()[4 + 2 * i + 1]), 6)
