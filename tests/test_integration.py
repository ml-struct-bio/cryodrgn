"""Running an experiment of training followed by downstream analyses."""

import pytest
import argparse
import os.path
import shutil
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

from cryodrgn.commands import (
    analyze,
    analyze_landscape,
    analyze_landscape_full,
    eval_images,
    eval_vol,
    graph_traversal,
    train_nn,
    train_vae,
    abinit_homo,
    abinit_het,
)

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


class TestFixedHetero:

    mrcs_file = f"{DATA_FOLDER}/hand.mrcs"
    poses_file = f"{DATA_FOLDER}/hand_rot_trans.pkl"

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
        """Produce standard analyses for a particular epoch."""
        args = analyze.add_args(argparse.ArgumentParser()).parse_args(
            [
                "output",
                "2",  # Epoch number to analyze - 0-indexed
                "--pc",
                "3",  # Number of principal component traversals to generate
                "--ksample",
                "20",  # Number of kmeans samples to generate
                "--vol-start-index",
                "1",
            ]
        )
        analyze.main(args)
        assert os.path.exists(os.path.join("output", "analyze.2"))

    def test_notebooks(self):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        os.chdir(os.path.join("output", "analyze.2"))

        # TODO: other notebooks don't work if --ctf not specified
        for ipynb in ["cryoDRGN_figures"]:
            assert os.path.exists(f"{ipynb}.ipynb")
            with open(f"{ipynb}.ipynb") as ff:
                nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

            ep = ExecutePreprocessor(timeout=600, kernel_name="python3")
            ep.preprocess(nb_in)

        os.chdir(os.path.join("..", ".."))

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


class TestAbinitHetero:

    mrcs_file = f"{DATA_FOLDER}/toy_projections.mrcs"
    ctf_file = f"{DATA_FOLDER}/test_ctf.pkl"

    def test_train_model(self):
        """Train the initial heterogeneous model."""

        args = abinit_het.add_args(argparse.ArgumentParser()).parse_args(
            [
                self.mrcs_file,
                "--ctf",
                self.ctf_file,
                "-o",
                "output",
                "--zdim",
                "4",
                "--lr",
                ".0001",
                "--enc-dim",
                "8",
                "--enc-layers",
                "2",
                "--dec-dim",
                "8",
                "--dec-layers",
                "2",
                "--pe-dim",
                "8",
                "--enc-only",
                "--t-extent",
                "4.0",
                "--t-ngrid",
                "2",
                "--pretrain",
                "1",
                "--num-epochs",
                "3",
                "--ps-freq",
                "2",
            ]
        )
        abinit_het.main(args)

    def test_analyze(self):
        """Produce standard analyses for a particular epoch."""
        args = analyze.add_args(argparse.ArgumentParser()).parse_args(
            [
                "output",
                "1",  # Epoch number to analyze - 0-indexed
                "--pc",
                "3",  # Number of principal component traversals to generate
                "--ksample",
                "20",  # Number of kmeans samples to generate
                "--vol-start-index",
                "1",
            ]
        )
        analyze.main(args)
        assert os.path.exists(os.path.join("output", "analyze.1"))

    def test_notebooks(self):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        os.chdir(os.path.join("output", "analyze.1"))

        for ipynb in ["cryoDRGN_figures"]:
            assert os.path.exists(f"{ipynb}.ipynb")
            with open(f"{ipynb}.ipynb") as ff:
                nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

            ep = ExecutePreprocessor(timeout=600, kernel_name="python3")
            ep.preprocess(nb_in)

        os.chdir(os.path.join("..", ".."))


@pytest.mark.parametrize(
    "star_particles",
    [os.path.join(DATA_FOLDER, "sta_testing_bin8.star")],
    ids=("sta-bin8",),
)
class TestSta:
    """Run reconstructions using particles from a .star file as input."""

    poses_file = os.path.join(DATA_FOLDER, "sta_pose.pkl")
    ctf_file = os.path.join(DATA_FOLDER, "sta_ctf.pkl")
    outdir = os.path.join(DATA_FOLDER, "output", "sta")

    def test_train_nn(self, star_particles):
        args = train_nn.add_args(argparse.ArgumentParser()).parse_args(
            [
                star_particles,
                "--datadir",
                DATA_FOLDER,
                "--poses",
                self.poses_file,
                "--ctf",
                self.ctf_file,
                "-o",
                self.outdir,
                "--dim",
                "256",
            ]
        )
        train_nn.main(args)

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

    def test_abinit_homo(self, star_particles):
        args = abinit_homo.add_args(argparse.ArgumentParser()).parse_args(
            [
                star_particles,
                "--datadir",
                DATA_FOLDER,
                "--ctf",
                self.ctf_file,
                "-o",
                self.outdir,
                "--dim",
                "4",
                "--layers",
                "2",
                "--t-extent",
                "4.0",
                "--t-ngrid",
                "2",
                "--pretrain=1",
                "--num-epochs",
                "3",
                "--ps-freq",
                "2",
            ]
        )
        abinit_homo.main(args)

    def test_abinit_het(self, star_particles):
        args = abinit_het.add_args(argparse.ArgumentParser()).parse_args(
            [
                star_particles,
                "--datadir",
                DATA_FOLDER,
                "--ctf",
                self.ctf_file,
                "--zdim",
                "8",
                "-o",
                self.outdir,
                "--enc-dim",
                "4",
                "--enc-layers",
                "2",
                "--dec-dim",
                "4",
                "--dec-layers",
                "2",
                "--pe-dim",
                "4",
                "--enc-only",
                "--t-extent",
                "4.0",
                "--t-ngrid",
                "2",
                "--pretrain=1",
                "--num-epochs",
                "3",
                "--ps-freq",
                "2",
            ]
        )
        abinit_het.main(args)
