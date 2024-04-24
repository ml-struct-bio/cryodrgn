"""Running an experiment of training followed by downstream analyses."""

import pytest
import argparse
import os.path
import shutil
from cryodrgn.commands import (
    analyze,
    analyze_landscape,
    analyze_landscape_full,
    eval_images,
    eval_vol,
    graph_traversal,
    train_vae,
)

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


class TestTrainVAE:

    mrcs_file = f"{DATA_FOLDER}/hand.mrcs"
    poses_file = f"{DATA_FOLDER}/hand_rot.pkl"

    def test_train_model(self):
        """Train the initial heterogeneous model."""

        args = train_vae.add_args(argparse.ArgumentParser()).parse_args(
            [
                self.mrcs_file,
                "-o",
                "output",
                "--lr",
                ".0001",
                "--num-epochs",
                "3",
                "--seed",
                "0",
                "--poses",
                self.poses_file,
                "--zdim",
                "10",
                "--pe-type",
                "gaussian",
                "--multigpu",
            ]
        )
        train_vae.main(args)

    def test_train_from_checkpoint(self):
        """Load a cached model and run for another epoch, now without --multigpu."""

        train_vae.main(
            train_vae.add_args(argparse.ArgumentParser()).parse_args(
                [
                    self.mrcs_file,
                    "-o",
                    "output",
                    "--lr",
                    ".0001",
                    "--num-epochs",
                    "4",
                    "--seed",
                    "0",
                    "--poses",
                    self.poses_file,
                    "--zdim",
                    "10",
                    "--pe-type",
                    "gaussian",
                    "--load",
                    "output/weights.2.pkl",
                ]
            )
        )

    def test_analyze(self):
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

    def test_landscape(self):
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

    def test_landscape_full(self):
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

    def test_graph_traversal(self):
        parser = argparse.ArgumentParser()
        graph_traversal.add_args(parser)
        args = parser.parse_args(
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
                "--outind",
                "output/graph_traversal_path.txt",
                "--outtxt",
                "output/graph_traversal_zpath.txt",
            ]
        )
        graph_traversal.main(args)

    def test_eval_volume(self):
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

    def test_eval_images(self):
        args = eval_images.add_args(argparse.ArgumentParser()).parse_args(
            [
                self.mrcs_file,
                "output/weights.3.pkl",
                "--config",
                "output/config.yaml",
                "-o",
                "output/out_eval_images_losses.pkl",
                "--out-z",
                "output/out_eval_images_z.pkl",
                "--poses",
                self.poses_file,
                "--log-interval",
                "1",
                "--verbose",
            ]
        )
        eval_images.main(args)


@pytest.mark.parametrize(
    "star_particles", [os.path.join(DATA_FOLDER, "sta_testing_bin8.star")]
)
class TestSta:

    poses_file = os.path.join(DATA_FOLDER, "sta_pose.pkl")
    ctf_file = os.path.join(DATA_FOLDER, "sta_ctf.pkl")
    outdir = os.path.join(DATA_FOLDER, "output", "sta")

    def test_train_vae(self, star_particles):
        args = train_vae.add_args(argparse.ArgumentParser()).parse_args(
            [
                star_particles,
                "--datadir",
                DATA_FOLDER,
                "--encode-mode",
                "tilt",
                "--poses",
                self.poses_file,
                "--ctf",
                self.ctf_file,
                "--zdim",
                "8",
                "-o",
                self.outdir,
                "--tdim",
                "256",
                "--enc-dim",
                "256",
                "--dec-dim",
                "256",
            ]
        )
        train_vae.main(args)
