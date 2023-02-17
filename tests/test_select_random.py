import os
import os.path
import argparse
from cryodrgn.commands_utils import select_random
from cryodrgn.utils import load_pkl

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


def test_select_random_n():
    os.makedirs("output", exist_ok=True)

    args = select_random.add_args(argparse.ArgumentParser()).parse_args(
        [
            "100",
            "-o",
            "output/selected_random.pkl",
            "-n",
            "42",
            "-s",
            "output/selected_random_inverted.pkl",
        ]
    )
    select_random.main(args)

    selected = load_pkl("output/selected_random.pkl")
    rejected = load_pkl("output/selected_random_inverted.pkl")
    assert len(selected) == 42  # ensure no duplicates
    assert len(rejected) == 58  # ensure no duplicates


def test_select_random_frac():
    os.makedirs("output", exist_ok=True)

    args = select_random.add_args(argparse.ArgumentParser()).parse_args(
        [
            "200",
            "-o",
            "output/selected_random.pkl",
            "--frac",
            "0.3",
            "-s",
            "output/selected_random_inverted.pkl",
        ]
    )
    select_random.main(args)

    selected = load_pkl("output/selected_random.pkl")
    rejected = load_pkl("output/selected_random_inverted.pkl")
    assert len(selected) == 60  # ensure no duplicates
    assert len(rejected) == 140  # ensure no duplicates
