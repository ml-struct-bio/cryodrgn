import os
import os.path
import argparse
import numpy as np
from cryodrgn.commands_utils import select_clusters
from cryodrgn.utils import load_pkl, save_pkl

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


def test_select_clusters():
    os.makedirs("output", exist_ok=True)

    # 100 labels for clusters 0-9
    labels = np.random.randint(10, size=100)
    save_pkl(labels, "output/labels.pik")

    selected_clusters = np.random.randint(10, size=3).tolist()
    selected_clusters_list = [str(c) for c in selected_clusters]
    args = select_clusters.add_args(argparse.ArgumentParser()).parse_args(
        ["output/labels.pik", "--sel"]
        + selected_clusters_list
        + ["-o", "output/selected_labels.pik"]
    )
    select_clusters.main(args)

    # Where (i.e. 0-indexed positions) in the labels array are the selected cluster indices found?
    x = load_pkl("output/selected_labels.pik")
    np.allclose(x, np.where(np.isin(labels, selected_clusters))[0])


def test_select_clusters_parent_ind():
    os.makedirs("output", exist_ok=True)

    # parent indices
    parent = np.random.randint(10, size=100)
    save_pkl(parent, "output/parent.pik")

    # 100 labels for clusters 0-9
    labels = np.random.randint(10, size=100)
    save_pkl(labels, "output/labels.pik")

    selected_clusters = np.random.randint(10, size=3).tolist()
    selected_clusters_list = [str(c) for c in selected_clusters]
    args = select_clusters.add_args(argparse.ArgumentParser()).parse_args(
        [
            "output/labels.pik",
            "--sel",
        ]
        + selected_clusters_list
        + [
            "-o",
            "output/selected_labels_parent.pik",
            "--parent-ind",
            "output/parent.pik",
            "--N-orig",
            "10",
        ]
    )
    select_clusters.main(args)

    # What are the indices in the parent array corresponding to the
    # 0-indexed positions in the labels array where the selected cluster indices are found?
    x = load_pkl("output/selected_labels_parent.pik")
    np.allclose(x, parent[np.where(np.isin(labels, selected_clusters))[0]])
