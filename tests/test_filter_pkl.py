import pytest
import os
import argparse
import numpy as np
from cryodrgn.commands_utils import filter_pkl
from cryodrgn.utils import load_pkl, save_pkl


@pytest.mark.parametrize("total_size", [100, 1000])
@pytest.mark.parametrize("ind_size", [10, 20, 80])
@pytest.mark.parametrize("selected_size", [4.2, 6.1])
@pytest.mark.parametrize("random_seed", [88, 99])
def test_select_clusters(tmpdir, total_size, ind_size, selected_size, random_seed):
    np.random.seed(random_seed)
    indices = np.random.randint(ind_size, size=total_size)
    test_lbl = "-".join(
        [str(total_size), str(ind_size), str(selected_size), str(random_seed)]
    )
    ind_fl = os.path.join(tmpdir, f"indices_{test_lbl}.pik")
    save_pkl(indices, ind_fl)

    # 42 selected 0-indexed positions in [0, 99]
    selected_positions = np.random.randint(
        total_size, size=int(selected_size * ind_size)
    )
    sel_fl = os.path.join(tmpdir, f"selected-positions_{test_lbl}.pik")
    save_pkl(selected_positions, sel_fl)

    out_fl = os.path.join(tmpdir, f"filtered-indices_{test_lbl}.pkl")
    args = filter_pkl.add_args(argparse.ArgumentParser()).parse_args(
        [ind_fl, "--ind", sel_fl, "-o", out_fl]
    )
    filter_pkl.main(args)

    # Indices values corresponding to the selected positions
    x = load_pkl(out_fl)
    assert np.allclose(x, indices[selected_positions])
