"""Unit tests of the cryodrgn fsc command."""
import pandas as pd
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
    out, err = run_command(f"cryodrgn_utils fsc {vol_file1} {vol_file2}")
    assert err == ""

    assert np.array_equal(
        np.array(out.split("\n")[2].split(), dtype="float"),
        np.array([0.0, 1.0]),
    )
    assert np.array_equal(
        np.array(out.split("\n")[4].split(), dtype="float"),
        np.array([0.0312, 0.8930]),
    )

    assert out.split("\n")[-3].split()[4] == "0.5:"
    assert round(float(out.split("\n")[-3].split()[5]), 6) == 2.064516
    assert out.split("\n")[-2].split()[4] == "0.143:"
    assert round(float(out.split("\n")[-3].split()[5]), 6) == 2.064516

    assert (
        round(
            sum(float(x) for x in out.split() if x[:2] == "0." and x[-1].isnumeric()),
            3,
        )
        == 34.584
    )
    assert len(out.split("\n")) == 38


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
    assert err == ""
    assert out0.split("\n")[-3].split()[4] == "0.5:"
    assert out0.split("\n")[-2].split()[4] == "0.143:"

    out, err = run_command(f"cryodrgn_utils fsc {vol_file1} {vol_file2} -o {fsc_file}")
    assert err == ""

    fsc_vals = pd.read_csv(fsc_file, dtype=float, sep=" ")
    assert fsc_vals.shape[0] * 2 + 20 == len(out0.split())
    for i, (_, (res_val, fsc_val)) in enumerate(fsc_vals.iterrows()):
        assert round(res_val, 4) == round(float(out0.split()[6 + 2 * i]), 6)
        assert round(fsc_val, 4) == round(float(out0.split()[7 + 2 * i]), 6)


@pytest.mark.parametrize(
    "train_dir",
    [
        {"dataset": "toy", "train_cmd": "train_nn", "epochs": 5, "seed": 2456},
    ],
    indirect=True,
)
@pytest.mark.parametrize("epochs", [(1, 2), (3, 4)])
def test_apply_mask(trained_dir, epochs: tuple[int, int]) -> None:
    vol_file1 = os.path.join(trained_dir.outdir, f"reconstruct.{epochs[0]}.mrc")
    vol_file2 = os.path.join(trained_dir.outdir, f"reconstruct.{epochs[1]}.mrc")
    mask_file = os.path.join(trained_dir.outdir, "mask.mrc")

    out, err = run_command(f"cryodrgn_utils gen_mask {vol_file1} {mask_file} --dist 3")
    assert err == ""
    assert 0.13 < float(out.split("\n")[0].split("Threshold=")[1]) < 0.16

    out, err = run_command(
        f"cryodrgn_utils fsc {vol_file1} {vol_file2} --mask {mask_file}"
    )
    assert err == ""

    assert out.split("\n")[-3].split()[4] == "0.5:"
    assert out.split("\n")[-2].split()[4] == "0.143:"
    assert round(float(out.split("\n")[7].split()[0]), 3) == 0.167
    assert 0.97 < float(out.split("\n")[7].split()[1]) < 0.99


@pytest.mark.parametrize(
    "train_dir",
    [
        {"dataset": "toy", "train_cmd": "train_nn", "epochs": 5, "seed": 2456},
    ],
    indirect=True,
)
@pytest.mark.parametrize("epochs", [(3, 4)])
def test_plotting(trained_dir, epochs: tuple[int, int]) -> None:
    vol_file1 = os.path.join(trained_dir.outdir, f"reconstruct.{epochs[0]}.mrc")
    vol_file2 = os.path.join(trained_dir.outdir, f"reconstruct.{epochs[1]}.mrc")

    plot_file = "fsc-plot.png"
    if os.path.exists(plot_file):
        os.unlink(plot_file)
    out, err = run_command(f"cryodrgn_utils fsc {vol_file1} {vol_file2} -p")
    assert err == ""
    assert os.path.exists(plot_file)
    os.unlink(plot_file)

    plot_file = os.path.join(trained_dir.outdir, "fsc-compare.png")
    if os.path.exists(plot_file):
        os.unlink(plot_file)
    out, err = run_command(f"cryodrgn_utils fsc {vol_file1} {vol_file2} -p {plot_file}")
    assert err == ""
    assert os.path.exists(plot_file)
    os.unlink(plot_file)

    fsc_file = os.path.join(trained_dir.outdir, "fsc.txt")
    out, err = run_command(
        f"cryodrgn_utils fsc {vol_file1} {vol_file2} -o {fsc_file} -p"
    )
    assert err == ""
    assert os.path.exists(fsc_file.replace(".txt", ".png"))
