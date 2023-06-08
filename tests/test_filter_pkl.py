import os
import os.path
import argparse
import numpy as np
from cryodrgn.commands_utils import filter_pkl
from cryodrgn.utils import load_pkl, save_pkl

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


def test_select_clusters():
    os.makedirs("output", exist_ok=True)

    # 100 random indices in [0, 9]
    indices = np.random.randint(10, size=100)
    save_pkl(indices, "output/indices.pik")

    # 42 selected 0-indexed positions in [0, 99]
    selected_positions = np.random.randint(100, size=42)
    save_pkl(selected_positions, "output/selected_positions.pik")

    args = filter_pkl.add_args(argparse.ArgumentParser()).parse_args(
        [
            "output/indices.pik",
            "--ind",
            "output/selected_positions.pik",
            "-o",
            "output/filtered_indices.pik",
        ]
    )
    filter_pkl.main(args)

    # Indices values corresponding to the selected positions
    x = load_pkl("output/filtered_indices.pik")
    assert np.allclose(x, indices[selected_positions])
