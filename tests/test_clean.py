"""Unit tests of the cryodrgn clean command."""

import pytest
import os
from cryodrgn.utils import run_command


@pytest.mark.parametrize("train_dir", [{"train_cmd": "train_nn"}], indirect=True)
@pytest.mark.parametrize("every_n", [1, 2, 3, 5])
def test_clean_here(train_dir, every_n: int) -> None:
    """Test that we can clean the output of the current directory."""

    os.chdir(train_dir.out_lbl)
    out, err = run_command(f"cryodrgn_utils clean -n {every_n} -d")
    os.chdir("..")
    assert out == f"\tWould remove {2 * (9 - 9 // every_n)} files!\n"
    assert err == ""
    assert train_dir.all_files_present

    os.chdir(train_dir.out_lbl)
    out, err = run_command(f"cryodrgn_utils clean -n {every_n}")
    os.chdir("..")

    assert out == f"\tRemoved {2 * (9 - 9 // every_n)} files!\n"
    assert err == ""

    assert not train_dir.epoch_cleaned(None)
    for epoch in range(10):
        assert train_dir.epoch_cleaned(epoch) == (epoch % every_n != 0)

    train_dir.replace_files()


@pytest.mark.parametrize(
    "train_dir",
    [
        {"train_cmd": "train_nn"},
        {"train_cmd": "train_vae"},
    ],
    indirect=True
)
@pytest.mark.parametrize("every_n", [1, 2, 5])
def test_clean_one(train_dir, every_n: int) -> None:
    """Test that we can clean the output of one directory."""

    out, err = run_command(f"cryodrgn_utils clean {train_dir.out_lbl} -n {every_n} -d")
    assert out == f"\tWould remove {2 * (9 - (9 // every_n))} files!\n"
    assert err == ""
    assert train_dir.all_files_present

    out, err = run_command(f"cryodrgn_utils clean {train_dir.out_lbl} -n {every_n}")
    assert out == f"\tRemoved {2 * (9 - (9 // every_n))} files!\n"
    assert err == ""

    assert not train_dir.epoch_cleaned(None)
    for epoch in range(10):
        assert train_dir.epoch_cleaned(epoch) == (epoch % every_n != 0)

    train_dir.replace_files()


@pytest.mark.parametrize(
    "train_dirs",
    [
        (("hand", "train_nn"), ("toy", "train_vae")),
    ],
    indirect=True,
)
@pytest.mark.parametrize("every_n", [2, 4])
def test_clean_two(train_dirs, every_n: int) -> None:
    """Test that we can clean the output of two directories."""

    dir_str = " ".join([train_dir.out_lbl for train_dir in train_dirs])
    out, err = run_command(f"cryodrgn_utils clean {dir_str} -n {every_n} -d")

    assert out == (
        f"\tWould remove {2 * (9 - (9 // every_n))} files!\n"
        f"\tWould remove {2 * (9 - (9 // every_n))} files!\n"
    )
    assert err == ""
    for train_dir in train_dirs:
        assert train_dir.all_files_present

    out, err = run_command(f"cryodrgn_utils clean {dir_str} -n {every_n}")
    assert out == (
        f"\tRemoved {2 * (9 - (9 // every_n))} files!\n"
        f"\tRemoved {2 * (9 - (9 // every_n))} files!\n"
    )
    assert err == ""

    for train_dir in train_dirs:
        assert not train_dir.epoch_cleaned(None)
        for epoch in range(10):
            assert train_dir.epoch_cleaned(epoch) == (epoch % every_n != 0)

        train_dir.replace_files()
