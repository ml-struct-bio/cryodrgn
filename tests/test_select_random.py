import pytest
import os
import argparse
from cryodrgn.commands_utils import select_random
from cryodrgn.utils import load_pkl


@pytest.mark.parametrize(
    "total_count, chosen_count", [(2, 1), (100, 1), (100, 42), (200, 177)]
)
def test_select_random_n(tmpdir, total_count, chosen_count):
    sel_file = os.path.join(tmpdir, "sel_random.pkl")
    rej_file = os.path.join(tmpdir, "sel_random_inverted.pkl")

    parser = argparse.ArgumentParser()
    select_random.add_args(parser)
    select_random.main(
        parser.parse_args(
            [str(total_count), "-o", sel_file, "-n", str(chosen_count), "-s", rej_file]
        )
    )

    selected = load_pkl(sel_file)
    rejected = load_pkl(rej_file)
    assert len(selected) == chosen_count  # ensure no duplicates
    assert len(rejected) == total_count - chosen_count  # ensure no duplicates


@pytest.mark.parametrize("total_count", [2, 3, 100, 500])
@pytest.mark.parametrize("chosen_frac", [0.01, 0.02, 0.5, 0.77])
def test_select_random_frac(tmpdir, total_count, chosen_frac):
    sel_file = os.path.join(tmpdir, "sel_random.pkl")
    rej_file = os.path.join(tmpdir, "sel_random_inverted.pkl")

    parser = argparse.ArgumentParser()
    select_random.add_args(parser)
    select_random.main(
        parser.parse_args(
            [
                str(total_count),
                "-o",
                sel_file,
                "--frac",
                str(chosen_frac),
                "-s",
                rej_file,
            ]
        )
    )

    selected = load_pkl(sel_file)
    rejected = load_pkl(rej_file)
    sel_count = int(total_count * chosen_frac)
    assert len(selected) == sel_count  # ensure no duplicates
    assert len(rejected) == total_count - sel_count  # ensure no duplicates
