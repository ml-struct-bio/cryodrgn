import argparse
import os.path
import shutil

import pytest

from cryodrgn.commands import (
    analyze,
    analyze_landscape,
    analyze_landscape_full,
    eval_vol,
    graph_traversal,
    train_vae,
)

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture
def mrcs_file():
    return f"{DATA_FOLDER}/hand.mrcs"


@pytest.fixture
def poses_file():
    return f"{DATA_FOLDER}/hand_rot.pkl"


def test_run(mrcs_file, poses_file):
    args = train_vae.add_args(argparse.ArgumentParser()).parse_args(
        [
            mrcs_file,
            "-o",
            "output",
            "--lr",
            ".0001",
            "--num-epochs",
            "3",
            "--seed",
            "0",
            "--poses",
            poses_file,
            "--zdim",
            "10",
            "--pe-type",
            "gaussian",
            "--multigpu",
        ]
    )
    train_vae.main(args)

    # Load a check-pointed model and run for another epoch, this time without --multigpu
    train_vae.main(
        train_vae.add_args(argparse.ArgumentParser()).parse_args(
            [
                mrcs_file,
                "-o",
                "output",
                "--lr",
                ".0001",
                "--num-epochs",
                "4",
                "--seed",
                "0",
                "--poses",
                poses_file,
                "--zdim",
                "10",
                "--pe-type",
                "gaussian",
                "--load",
                "output/weights.2.pkl",
            ]
        )
    )

    args = analyze.add_args(argparse.ArgumentParser()).parse_args(
        [
            "output",
            "2",  # Epoch number to analyze - 0-indexed
            "--pc",
            "3",  # Number of principal component traversals to generate
            "--ksample",
            "14",  # Number of kmeans samples to generate
            "--vol-start-index",
            "1",
        ]
    )
    analyze.main(args)

    shutil.rmtree("output/landscape.3", ignore_errors=True)
    args = analyze_landscape.add_args(argparse.ArgumentParser()).parse_args(
        [
            "output",
            "2",  # Epoch number to analyze - 0-indexed
            "--sketch-size",
            "10",  # Number of volumes to generate for analysis
            "--downsample",
            "64",
            "--pc-dim",
            "5",
            "--vol-start-index",
            "1",
        ]
    )
    shutil.rmtree("output/landscape.2", ignore_errors=True)
    analyze_landscape.main(args)

    args = analyze_landscape_full.add_args(argparse.ArgumentParser()).parse_args(
        [
            "output",
            "2",  # Epoch number to analyze - 0-indexed
            "-N",
            "10",  # Number of training volumes to generate
            "--downsample",
            "64",
        ]
    )
    shutil.rmtree("output/landscape.2/landscape_full", ignore_errors=True)
    analyze_landscape_full.main(args)

    args = graph_traversal.add_args(argparse.ArgumentParser()).parse_args(
        [
            "output/z.3.pkl",
            "--anchors",
            "22",
            "49",
            "53",
            "6",
            "27",
            "95",
            "64",
            "81",
            "44",
            "58",
            "75",
            "67",
            "9",
            "89",
            "-o",
            "output/graph_traversal_path.txt",
            "--out-z",
            "output/graph_traversal_zpath.txt",
        ]
    )
    graph_traversal.main(args)

    args = eval_vol.add_args(argparse.ArgumentParser()).parse_args(
        [
            "output/weights.3.pkl",
            "--config",
            "output/config.yaml",
            "--zfile",
            "output/graph_traversal_zpath.txt",
            "-o",
            "output/eval_vols",
        ]
    )
    eval_vol.main(args)
