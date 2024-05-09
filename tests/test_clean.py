"""Unit tests of the cryodrgn clean command."""

import pytest
import os
from cryodrgn.utils import run_command


@pytest.mark.parametrize("train_dir", [{"train_cmd": "train_nn"}], indirect=True)
@pytest.mark.parametrize("every_n", [1, 2, 3, 5])
def test_clean_here(trained_dir, every_n: int) -> None:
    """Test that we can clean the output of the current directory."""

    os.chdir(trained_dir.out_lbl)
    out, err = run_command(f"cryodrgn_utils clean -n {every_n} -d")
    os.chdir("..")
    assert out == f"\tWould remove {2 * (9 - 9 // every_n)} files!\n"
    assert err == ""
    assert trained_dir.all_files_present

    os.chdir(trained_dir.out_lbl)
    out, err = run_command(f"cryodrgn_utils clean -n {every_n}")
    os.chdir("..")

    assert out == f"\tRemoved {2 * (9 - 9 // every_n)} files!\n"
    assert err == ""

    assert not trained_dir.epoch_cleaned(None)
    for epoch in range(10):
        assert trained_dir.epoch_cleaned(epoch) == (epoch % every_n != 0)


@pytest.mark.parametrize(
    "train_dir", [{"train_cmd": "train_nn"}, {"train_cmd": "train_vae"}], indirect=True
)
@pytest.mark.parametrize("every_n", [1, 2, 5])
def test_clean_one(trained_dir, every_n: int) -> None:
    """Test that we can clean the output of one directory."""

    out, err = run_command(
        f"cryodrgn_utils clean {os.path.relpath(trained_dir.out_lbl)} -n {every_n} -d"
    )

    rmv_count = 2 * (trained_dir.epochs - 1 - ((trained_dir.epochs - 1) // every_n))
    assert out == f"\tWould remove {rmv_count} files!\n"
    assert err == ""
    assert trained_dir.all_files_present

    out, err = run_command(
        f"cryodrgn_utils clean {os.path.relpath(trained_dir.out_lbl)} -n {every_n}"
    )
    assert out == f"\tRemoved {rmv_count} files!\n"
    assert err == ""

    assert not trained_dir.epoch_cleaned(None)
    for epoch in range(trained_dir.epochs):
        assert trained_dir.epoch_cleaned(epoch) == (epoch % every_n != 0)


@pytest.mark.parametrize(
    "train_dirs",
    [
        (
            {"dataset": "hand", "train_cmd": "train_vae", "epochs": 4},
            {"dataset": "toy", "train_cmd": "train_nn", "epochs": 11},
        )
    ],
    indirect=True,
)
@pytest.mark.parametrize("every_n", [2, 3])
def test_clean_two(trained_dirs, every_n: int) -> None:
    """Test that we can clean the output of two directories."""

    dir_str = " ".join(
        [os.path.relpath(trained_dir.out_lbl) for trained_dir in trained_dirs]
    )
    out, err = run_command(f"cryodrgn_utils clean {dir_str} -n {every_n} -d")

    rmv_counts = [
        2 * (trained_dir.epochs - 1 - ((trained_dir.epochs - 1) // every_n))
        for trained_dir in trained_dirs
    ]

    assert out == (
        f"\tWould remove {rmv_counts[0]} files!\n"
        f"\tWould remove {rmv_counts[1]} files!\n"
    )
    assert err == ""
    for trained_dir in trained_dirs:
        assert trained_dir.all_files_present

    out, err = run_command(f"cryodrgn_utils clean {dir_str} -n {every_n}")
    assert out == (
        f"\tRemoved {rmv_counts[0]} files!\n" f"\tRemoved {rmv_counts[1]} files!\n"
    )
    assert err == ""

    for trained_dir in trained_dirs:
        assert not trained_dir.epoch_cleaned(None)
        for epoch in range(trained_dir.epochs):
            assert trained_dir.epoch_cleaned(epoch) == (epoch % every_n != 0)
