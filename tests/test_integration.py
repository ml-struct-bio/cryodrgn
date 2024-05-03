"""Running an experiment of training followed by downstream analyses."""

import pytest
import argparse
import os.path
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

    mrcs_file = os.path.join(DATA_FOLDER, "hand.mrcs")
    poses_file = os.path.join(DATA_FOLDER, "hand_rot_trans.pkl")

    def test_train_model(self, outdir):
        """Train the initial heterogeneous model."""

        args = train_vae.add_args(argparse.ArgumentParser()).parse_args(
            [
                self.mrcs_file,
                "-o",
                outdir,
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

    def test_train_from_checkpoint(self, outdir):
        """Load a cached model and run for another epoch, now without --multigpu."""

        train_vae.main(
            train_vae.add_args(argparse.ArgumentParser()).parse_args(
                [
                    self.mrcs_file,
                    "-o",
                    outdir,
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
                    os.path.join(outdir, "weights.2.pkl"),
                ]
            )
        )

    def test_analyze(self, outdir):
        """Produce standard analyses for a particular epoch."""
        args = analyze.add_args(argparse.ArgumentParser()).parse_args([outdir, "2"])
        analyze.main(args)
        assert os.path.exists(os.path.join(outdir, "analyze.2"))

    @pytest.mark.parametrize("nb_lbl", ["cryoDRGN_figures"])
    def test_notebooks(self, outdir, nb_lbl):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        os.chdir(os.path.join(outdir, "analyze.2"))
        assert os.path.exists(f"{nb_lbl}.ipynb")

        with open(f"{nb_lbl}.ipynb") as ff:
            nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

        ExecutePreprocessor(timeout=600, kernel_name="python3").preprocess(nb_in)
        os.chdir(os.path.join("..", ".."))

    def test_landscape(self, outdir):
        args = analyze_landscape.add_args(argparse.ArgumentParser()).parse_args(
            [
                outdir,
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
        analyze_landscape.main(args)

    def test_landscape_full(self, outdir):
        args = analyze_landscape_full.add_args(argparse.ArgumentParser()).parse_args(
            [
                outdir,
                "2",  # Epoch number to analyze - 0-indexed
                "-N",
                "10",  # Number of training volumes to generate
                "--downsample",
                "64",
            ]
        )
        analyze_landscape_full.main(args)

    def test_graph_traversal(self, outdir):
        parser = argparse.ArgumentParser()
        graph_traversal.add_args(parser)
        args = parser.parse_args(
            [
                os.path.join(outdir, "z.3.pkl"),
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
                os.path.join(outdir, "graph_traversal_path.txt"),
                "--outtxt",
                os.path.join(outdir, "graph_traversal_zpath.txt"),
            ]
        )
        graph_traversal.main(args)

    def test_eval_volume(self, outdir):
        args = eval_vol.add_args(argparse.ArgumentParser()).parse_args(
            [
                os.path.join(outdir, "weights.3.pkl"),
                "--config",
                os.path.join(outdir, "config.yaml"),
                "--zfile",
                os.path.join(outdir, "graph_traversal_zpath.txt"),
                "-o",
                os.path.join(outdir, "eval_vols"),
            ]
        )
        eval_vol.main(args)

    def test_eval_images(self, outdir):
        args = eval_images.add_args(argparse.ArgumentParser()).parse_args(
            [
                self.mrcs_file,
                os.path.join(outdir, "weights.3.pkl"),
                "--config",
                os.path.join(outdir, "config.yaml"),
                "-o",
                os.path.join(outdir, "out_eval_images_losses.pkl"),
                "--out-z",
                os.path.join(outdir, "out_eval_images_z.pkl"),
                "--poses",
                self.poses_file,
                "--log-interval",
                "1",
                "--verbose",
            ]
        )
        eval_images.main(args)


class TestAbinitHetero:

    mrcs_file = os.path.join(DATA_FOLDER, "toy_projections.mrcs")
    ctf_file = os.path.join(DATA_FOLDER, "test_ctf.pkl")

    def test_train_model(self, outdir):
        """Train the initial heterogeneous model."""

        args = abinit_het.add_args(argparse.ArgumentParser()).parse_args(
            [
                self.mrcs_file,
                "--ctf",
                self.ctf_file,
                "-o",
                outdir,
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

    def test_analyze(self, outdir):
        """Produce standard analyses for a particular epoch."""
        args = analyze.add_args(argparse.ArgumentParser()).parse_args(
            [
                outdir,
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
        assert os.path.exists(os.path.join(outdir, "analyze.1"))

    @pytest.mark.parametrize("nb_lbl", ["cryoDRGN_figures"])
    def test_notebooks(self, outdir, nb_lbl):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        os.chdir(os.path.join(outdir, "analyze.1"))
        assert os.path.exists(f"{nb_lbl}.ipynb"), "Upstream tests have failed!"

        with open(f"{nb_lbl}.ipynb") as ff:
            nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

        ExecutePreprocessor(timeout=600, kernel_name="python3").preprocess(nb_in)
        os.chdir(os.path.join("..", ".."))


@pytest.mark.parametrize(
    "star_particles",
    [os.path.join(DATA_FOLDER, "sta_testing_bin8.star")],
    ids=("sta-bin8",),
)
@pytest.mark.parametrize(
    "indices_file",
    [
        None,
        os.path.join(DATA_FOLDER, "ind3.pkl"),
        os.path.join(DATA_FOLDER, "ind3-numpy.pkl"),
    ],
    ids=("no-ind", "ind3", "ind3-numpy"),
)
class TestStarFixedHomo:
    """Run reconstructions using particles from a .star file as input."""

    poses_file = os.path.join(DATA_FOLDER, "sta_pose.pkl")
    ctf_file = os.path.join(DATA_FOLDER, "sta_ctf.pkl")

    def test_train_model(self, outdir, star_particles, indices_file):
        args = [
            star_particles,
            "--datadir",
            DATA_FOLDER,
            "--poses",
            self.poses_file,
            "--ctf",
            self.ctf_file,
            "-o",
            outdir,
            "--dim",
            "32",
        ]
        if indices_file is not None:
            args += ["--ind", indices_file]

        args = train_nn.add_args(argparse.ArgumentParser()).parse_args(args)
        train_nn.main(args)


@pytest.mark.parametrize(
    "star_particles",
    [os.path.join(DATA_FOLDER, "sta_testing_bin8.star")],
    ids=("sta-bin8",),
)
class TestStarFixedHetero:
    """Run reconstructions using particles from a .star file as input."""

    poses_file = os.path.join(DATA_FOLDER, "sta_pose.pkl")
    ctf_file = os.path.join(DATA_FOLDER, "sta_ctf.pkl")

    def test_train_model(self, outdir, star_particles):
        args = [
            star_particles,
            "--datadir",
            DATA_FOLDER,
            "--encode-mode",
            "tilt",
            "--poses",
            self.poses_file,
            "--ctf",
            self.ctf_file,
            "--num-epochs",
            "5",
            "--zdim",
            "4",
            "-o",
            outdir,
            "--tdim",
            "16",
            "--enc-dim",
            "16",
            "--dec-dim",
            "16",
        ]
        args = train_vae.add_args(argparse.ArgumentParser()).parse_args(args)
        train_vae.main(args)

    def test_analyze(self, outdir, star_particles):
        """Produce standard analyses for a particular epoch."""
        args = analyze.add_args(argparse.ArgumentParser()).parse_args(
            [
                outdir,
                "4",  # Epoch number to analyze - 0-indexed
                "--pc",
                "3",  # Number of principal component traversals to generate
                "--ksample",
                "2",  # Number of kmeans samples to generate
                "--vol-start-index",
                "1",
            ]
        )
        analyze.main(args)
        assert os.path.exists(os.path.join(outdir, "analyze.4"))

    @pytest.mark.parametrize("nb_lbl", ["cryoDRGN_figures", "cryoDRGN_ET_viz"])
    def test_notebooks(self, outdir, nb_lbl, star_particles):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        os.chdir(os.path.join(outdir, "analyze.4"))
        assert os.path.exists(f"{nb_lbl}.ipynb"), "Upstream tests have failed!"

        with open(f"{nb_lbl}.ipynb") as ff:
            nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

        ExecutePreprocessor(timeout=600, kernel_name="python3").preprocess(nb_in)
        os.chdir(os.path.join("..", ".."))


@pytest.mark.parametrize(
    "star_particles",
    [os.path.join(DATA_FOLDER, "sta_testing_bin8.star")],
    ids=("sta-bin8",),
)
@pytest.mark.parametrize(
    "indices_file",
    [
        None,
        os.path.join(DATA_FOLDER, "ind3.pkl"),
        os.path.join(DATA_FOLDER, "ind3-numpy.pkl"),
    ],
    ids=("no-ind", "ind3", "ind3-numpy"),
)
class TestStarAbinitHomo:

    ctf_file = os.path.join(DATA_FOLDER, "sta_ctf.pkl")

    def test_train_model(self, outdir, star_particles, indices_file):
        args = [
            star_particles,
            "--datadir",
            DATA_FOLDER,
            "--ctf",
            self.ctf_file,
            "-o",
            outdir,
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
        if indices_file is not None:
            args += ["--ind", indices_file]

        args = abinit_homo.add_args(argparse.ArgumentParser()).parse_args(args)
        abinit_homo.main(args)


@pytest.mark.parametrize(
    "star_particles",
    [os.path.join(DATA_FOLDER, "sta_testing_bin8.star")],
    ids=("sta-bin8",),
)
@pytest.mark.parametrize(
    "indices_file",
    [
        None,
        os.path.join(DATA_FOLDER, "ind3.pkl"),
        os.path.join(DATA_FOLDER, "ind3-numpy.pkl"),
    ],
    ids=("no-ind", "ind3", "ind3-numpy"),
)
class TestStarAbinitHetero:

    ctf_file = os.path.join(DATA_FOLDER, "sta_ctf.pkl")

    def test_train_model(self, outdir, star_particles, indices_file):
        args = [
            star_particles,
            "--datadir",
            DATA_FOLDER,
            "--ctf",
            self.ctf_file,
            "--zdim",
            "8",
            "-o",
            outdir,
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
        if indices_file is not None:
            args += ["--ind", indices_file]

        args = abinit_het.add_args(argparse.ArgumentParser()).parse_args(args)
        abinit_het.main(args)


@pytest.mark.parametrize(
    "indices_file",
    [
        os.path.join(DATA_FOLDER, "ind100.pkl"),
        os.path.join(DATA_FOLDER, "ind100-rand.pkl"),
    ],
    ids=("first-100", "random-100"),
)
class TestIterativeFiltering:

    mrcs_file = os.path.join(DATA_FOLDER, "toy_projections.mrcs")
    poses_file = os.path.join(DATA_FOLDER, "toy_rot_trans.pkl")
    ctf_file = os.path.join(DATA_FOLDER, "test_ctf.pkl")

    def test_train_model(self, outdir, indices_file):
        """Train the initial heterogeneous model."""
        args = train_vae.add_args(argparse.ArgumentParser()).parse_args(
            [
                self.mrcs_file,
                "-o",
                outdir,
                "--ctf",
                self.ctf_file,
                "--ind",
                indices_file,
                "--num-epochs",
                "3",
                "--poses",
                self.poses_file,
                "--zdim",
                "4",
                "--pe-type",
                "gaussian",
                "--no-amp",
            ]
        )
        train_vae.main(args)

    def test_analyze(self, outdir, indices_file):
        """Produce standard analyses for a particular epoch."""
        assert os.path.exists(
            os.path.join(outdir, "weights.2.pkl")
        ), "Upstream tests have failed!"

        args = analyze.add_args(argparse.ArgumentParser()).parse_args(
            [
                outdir,
                "2",  # Epoch number to analyze - 0-indexed
                "--pc",
                "3",  # Number of principal component traversals to generate
                "--ksample",
                "8",  # Number of kmeans samples to generate
                "--vol-start-index",
                "1",
            ]
        )
        analyze.main(args)
        assert os.path.exists(os.path.join(outdir, "analyze.2"))

    @pytest.mark.parametrize(
        "nb_lbl", ["cryoDRGN_figures", "cryoDRGN_filtering", "cryoDRGN_viz"]
    )
    def test_notebooks(self, outdir, nb_lbl, indices_file):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        os.chdir(os.path.join(outdir, "analyze.2"))
        assert os.path.exists(f"{nb_lbl}.ipynb"), "Upstream tests have failed!"

        with open(f"{nb_lbl}.ipynb") as ff:
            nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

        ExecutePreprocessor(timeout=600, kernel_name="python3").preprocess(nb_in)
        os.chdir(os.path.join("..", ".."))
