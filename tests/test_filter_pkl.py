import pytest
import os
import argparse
import numpy as np
from cryodrgn.commands_utils import filter_pkl, concat_pkls
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
    parser = argparse.ArgumentParser()
    filter_pkl.add_args(parser)
    args = parser.parse_args([ind_fl, "--ind", sel_fl, "-o", out_fl])
    filter_pkl.main(args)

    # Indices values corresponding to the selected positions
    x = load_pkl(out_fl)
    assert np.allclose(x, indices[selected_positions])


@pytest.mark.parametrize("ctf", ["CTF-Test"], indirect=True)
@pytest.mark.parametrize(
    "ind1, ind2", [([13, 21, 71], [1, 2, 3]), ([5, 10, 15], [4, 5, 6])]
)
def test_filter_ctf_pkl(tmpdir, ctf, ind1, ind2):
    ind1_fl = os.path.join(tmpdir, "ind1.pkl")
    save_pkl(ind1, ind1_fl)
    ind2_fl = os.path.join(tmpdir, "ind2.pkl")
    save_pkl(ind2, ind2_fl)
    out_pkl1 = os.path.join(tmpdir, "filtered-ctf1.pkl")
    out_pkl2 = os.path.join(tmpdir, "filtered-ctf2.pkl")

    parser = argparse.ArgumentParser()
    filter_pkl.add_args(parser)
    args = parser.parse_args([ctf.path, "--ind", ind1_fl, "-o", out_pkl1])
    filter_pkl.main(args)
    args = parser.parse_args([ctf.path, "--ind", ind2_fl, "-o", out_pkl2])
    filter_pkl.main(args)

    ctf_params = load_pkl(ctf.path)
    x1, x2 = load_pkl(out_pkl1), load_pkl(out_pkl2)
    assert isinstance(x1, np.ndarray)
    assert x1.shape == (len(ind1), ctf_params.shape[1])
    assert isinstance(x2, np.ndarray)
    assert x2.shape == (len(ind2), ctf_params.shape[1])

    concat_pkl = os.path.join(tmpdir, "concat-ctf.pkl")
    parser = argparse.ArgumentParser()
    concat_pkls.add_args(parser)
    args = parser.parse_args([out_pkl1, out_pkl2, "-o", concat_pkl])
    concat_pkls.main(args)
    x = load_pkl(concat_pkl)
    assert isinstance(x, np.ndarray)
    assert x.shape == (len(ind1) + len(ind2), ctf_params.shape[1])

    assert np.allclose(x, np.concatenate([x1, x2]))


@pytest.mark.parametrize("poses", ["toy-poses"], indirect=True)
@pytest.mark.parametrize(
    "ind1, ind2", [([10, 20, 80], [1, 2, 3]), ([5, 10, 15], [4, 5, 6])]
)
def test_filter_pose_pkl(tmpdir, poses, ind1, ind2):
    ind1_fl = os.path.join(tmpdir, "ind1.pkl")
    save_pkl(ind1, ind1_fl)
    ind2_fl = os.path.join(tmpdir, "ind2.pkl")
    save_pkl(ind2, ind2_fl)
    out_pkl1 = os.path.join(tmpdir, "filtered-pose1.pkl")
    out_pkl2 = os.path.join(tmpdir, "filtered-pose2.pkl")

    parser = argparse.ArgumentParser()
    filter_pkl.add_args(parser)
    args = parser.parse_args([poses.path, "--ind", ind1_fl, "-o", out_pkl1])
    filter_pkl.main(args)
    args = parser.parse_args([poses.path, "--ind", ind2_fl, "-o", out_pkl2])
    filter_pkl.main(args)

    x1, x2 = load_pkl(out_pkl1), load_pkl(out_pkl2)
    assert isinstance(x1, tuple)
    assert len(x1) == 2
    assert x1[0].shape == (len(ind1), 3, 3)
    assert x1[1].shape == (len(ind1), 2)
    assert isinstance(x2, tuple)
    assert len(x2) == 2
    assert x2[0].shape == (len(ind2), 3, 3)
    assert x2[1].shape == (len(ind2), 2)

    concat_pkl = os.path.join(tmpdir, "concat-pose.pkl")
    parser = argparse.ArgumentParser()
    concat_pkls.add_args(parser)
    args = parser.parse_args([out_pkl1, out_pkl2, "-o", concat_pkl])
    concat_pkls.main(args)
    x = load_pkl(concat_pkl)
    assert isinstance(x, tuple)
    assert len(x) == 2
    assert x[0].shape == (len(ind1) + len(ind2), 3, 3)
    assert x[1].shape == (len(ind1) + len(ind2), 2)

    assert np.allclose(x[0], np.concatenate([x1[0], x2[0]]))
    assert np.allclose(x[1], np.concatenate([x1[1], x2[1]]))
