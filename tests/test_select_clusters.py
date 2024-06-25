import pytest
import os
import argparse
import numpy as np
from cryodrgn.commands_utils import select_clusters
from cryodrgn.utils import load_pkl, save_pkl


@pytest.mark.parametrize(
    "cluster_count, chosen_count", [(2, 1), (3, 1), (5, 3), (10, 3)]
)
def test_select_clusters(tmpdir, cluster_count, chosen_count):

    # 100 labels for clusters 0-9
    labels = np.random.randint(cluster_count, size=100)
    test_lbl = "-".join([str(cluster_count), str(chosen_count)])
    lbl_file = os.path.join(tmpdir, f"labels_{test_lbl}.pik")
    save_pkl(labels, lbl_file)

    selected_clusters = np.random.randint(cluster_count, size=chosen_count).tolist()
    selected_clusters_list = [str(c) for c in selected_clusters]
    selected_file = os.path.join(tmpdir, f"selected_clusters_{test_lbl}.pik")
    args = select_clusters.add_args(argparse.ArgumentParser()).parse_args(
        [lbl_file, "--sel"] + selected_clusters_list + ["-o", selected_file]
    )
    select_clusters.main(args)

    # where (i.e. 0-indexed positions) in the labels array
    # are the selected cluster indices found?
    x = load_pkl(selected_file)
    np.allclose(x, np.where(np.isin(labels, selected_clusters))[0])


@pytest.mark.parametrize(
    "cluster_count, chosen_count", [(2, 1), (3, 1), (5, 3), (10, 3)]
)
def test_select_clusters_parent_ind(tmpdir, cluster_count, chosen_count):
    # parent indices
    parent = np.random.randint(cluster_count, size=500)
    test_lbl = "-".join([str(cluster_count), str(chosen_count)])
    parent_file = os.path.join(tmpdir, f"parent_{test_lbl}.pik")
    save_pkl(parent, parent_file)

    # 100 labels for clusters 0-9
    labels = np.random.randint(cluster_count, size=500)
    label_file = os.path.join(tmpdir, f"labels_{test_lbl}.pik")
    save_pkl(labels, label_file)

    selected_clusters = np.random.randint(cluster_count, size=chosen_count).tolist()
    selected_clusters_list = [str(c) for c in selected_clusters]
    selected_file = os.path.join(tmpdir, f"selected_labels_{test_lbl}.pik")
    args = select_clusters.add_args(argparse.ArgumentParser()).parse_args(
        [label_file, "--sel"]
        + selected_clusters_list
        + [
            "-o",
            selected_file,
            "--parent-ind",
            parent_file,
            "--N-orig",
            str(cluster_count),
        ]
    )
    select_clusters.main(args)

    # What are the indices in the parent array corresponding to the
    # 0-indexed positions in the labels array where the selected cluster indices are found?
    x = load_pkl(selected_file)
    np.allclose(x, parent[np.where(np.isin(labels, selected_clusters))[0]])
