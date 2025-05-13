"""Unit tests of the cryodrgn fsc command."""

import pytest
import os
import pandas as pd
from cryodrgn.utils import run_command


@pytest.mark.parametrize(
    "train_dir",
    [
        {"dataset": "hand", "train_cmd": "train_nn", "epochs": 3, "seed": 2456},
        {"dataset": "toy", "train_cmd": "train_nn", "epochs": 3, "seed": 2456},
    ],
    indirect=True,
)
@pytest.mark.parametrize("epochs", [(1, 2), (1, 3), (2, 3)])
def test_output_file(trained_dir, epochs: tuple[int, int]) -> None:
    vol_file1 = trained_dir.volume_files[epochs[0]]
    vol_file2 = trained_dir.volume_files[epochs[1]]
    fsc_file = os.path.join(trained_dir.outdir, "fsc.txt")

    out0, err = run_command(f"cryodrgn_utils fsc {vol_file1} {vol_file2}")
    assert err == ""
    assert out0.split("\n")[1].split()[-3] == "FSC=0.5:"
    assert out0.split("\n")[2].split()[-3] == "FSC=0.143:"

    out, err = run_command(f"cryodrgn_utils fsc {vol_file1} {vol_file2} -o {fsc_file}")
    assert err == ""

    fsc_vals = pd.read_csv(fsc_file, dtype=float, sep=" ")
    assert fsc_vals.shape[0] * 2 + 34 == len(out0.split())
    for i, (_, (res_val, fsc_val)) in enumerate(fsc_vals.iterrows()):
        assert round(res_val, 4) == round(float(out0.split()[34 + 2 * i]), 6)
        assert round(fsc_val, 4) == round(float(out0.split()[35 + 2 * i]), 6)


@pytest.mark.parametrize(
    "train_dir",
    [
        {"dataset": "toy", "train_cmd": "train_nn", "epochs": 5, "seed": 2456},
    ],
    indirect=True,
)
@pytest.mark.parametrize("epochs", [(1, 2), (3, 4)])
def test_apply_mask(trained_dir, epochs: tuple[int, int]) -> None:
    vol_file1 = trained_dir.volume_files[epochs[0]]
    vol_file2 = trained_dir.volume_files[epochs[1]]
    mask_file = os.path.join(trained_dir.outdir, "mask.mrc")

    out, err = run_command(f"cryodrgn_utils gen_mask {vol_file1} {mask_file} --dist 3")
    assert err == ""
    assert float(out.split("\n")[0].split("Threshold=")[1]) < 0.17

    out, err = run_command(
        f"cryodrgn_utils fsc {vol_file1} {vol_file2} --mask {mask_file}"
    )
    assert err == ""

    assert out.split("\n")[1].split()[-3] == "FSC=0.5:"
    assert out.split("\n")[2].split()[-3] == "FSC=0.143:"
    assert 0.97 < float(out.split("\n")[10].split()[1]) < 0.99


@pytest.mark.parametrize(
    "train_dir",
    [
        {"dataset": "toy", "train_cmd": "train_nn", "epochs": 4, "seed": 2456},
    ],
    indirect=True,
)
def test_apply_phase_randomization(trained_dir) -> None:
    vol_file1 = trained_dir.volume_files[2]
    vol_file2 = trained_dir.volume_files[3]

    out, err = run_command(f"cryodrgn_utils fsc {vol_file1} {vol_file2} --corrected 16")
    assert err == ""
    assert out.split("\n")[1].split()[-3] == "FSC=0.5:"
    assert out.split("\n")[2].split()[-3] == "FSC=0.143:"
    assert 0.95 < float(out.split("\n")[10].split()[1])


@pytest.mark.parametrize(
    "train_dir",
    [
        {"dataset": "toy", "train_cmd": "train_nn", "epochs": 4, "seed": 1011},
    ],
    indirect=True,
)
def test_use_cryosparc_correction(trained_dir) -> None:
    vol_file1 = trained_dir.volume_files[2]
    vol_file2 = trained_dir.volume_files[3]
    ref_vol = trained_dir.volume_files[4]
    mask_file = os.path.join(trained_dir.outdir, "mask.mrc")

    out, err = run_command(
        f"cryodrgn_utils fsc {vol_file1} {vol_file2} --ref-volume {ref_vol} "
    )
    assert err == ""
    assert out.split("\n")[8].split()[-3] == "FSC=0.5:"
    assert out.split("\n")[9].split()[-3] == "FSC=0.143:"
    fsc_base = float(out.split("\n")[15].split()[3])

    out, err = run_command(
        f"cryodrgn_utils gen_mask {vol_file1} {mask_file} --dist 6 --dilate 6"
    )
    assert err == ""

    out, err = run_command(
        f"cryodrgn_utils fsc {vol_file1} {vol_file2} "
        f"--ref-volume {ref_vol} --mask {mask_file}"
    )
    assert err == ""
    assert out.split("\n")[5].split()[-3] == "FSC=0.5:"
    assert out.split("\n")[6].split()[-3] == "FSC=0.143:"
    assert fsc_base == float(out.split("\n")[12].split()[3])


@pytest.mark.parametrize(
    "train_dir",
    [
        {"dataset": "toy", "train_cmd": "train_nn", "epochs": 4, "seed": 2456},
    ],
    indirect=True,
)
@pytest.mark.parametrize("epochs", [(3, 4)])
def test_plotting(trained_dir, epochs: tuple[int, int]) -> None:
    vol_file1 = trained_dir.volume_files[epochs[0]]
    vol_file2 = trained_dir.volume_files[epochs[1]]

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
