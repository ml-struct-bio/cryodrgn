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

    mrcs_file = os.path.join(DATA_FOLDER, "hand.mrcs")
    poses_file = os.path.join(DATA_FOLDER, "hand_rot_trans.pkl")

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

        shutil.rmtree("output")


class TestAbinitHetero:

    mrcs_file = os.path.join(DATA_FOLDER, "toy_projections.mrcs")
    ctf_file = os.path.join(DATA_FOLDER, "test_ctf.pkl")

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
        shutil.rmtree("output")


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
    ],
    ids=("no-ind", "ind3"),
)
class TestSta:
    """Run reconstructions using particles from a .star file as input."""

    poses_file = os.path.join(DATA_FOLDER, "sta_pose.pkl")
    ctf_file = os.path.join(DATA_FOLDER, "sta_ctf.pkl")

    def test_train_nn(self, star_particles, indices_file):
        if indices_file is None:
            outdir = os.path.join("output", "sta")
        else:
            outdir = os.path.join(
                "output", os.path.splitext(os.path.basename(indices_file))[0]
            )

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
            "256",
        ]

        if indices_file is not None:
            args += ["--ind", indices_file]

        args = train_nn.add_args(argparse.ArgumentParser()).parse_args(args)
        train_nn.main(args)

    def test_train_vae(self, star_particles, indices_file):
        if indices_file is None:
            outdir = os.path.join("output", "sta")
        else:
            outdir = os.path.join(
                "output", os.path.splitext(os.path.basename(indices_file))[0]
            )

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
            "--zdim",
            "8",
            "-o",
            outdir,
            "--tdim",
            "256",
            "--enc-dim",
            "256",
            "--dec-dim",
            "256",
        ]
        if indices_file is not None:
            args += ["--ind", indices_file]

        args = train_vae.add_args(argparse.ArgumentParser()).parse_args(args)
        train_vae.main(args)

    def test_abinit_homo(self, star_particles, indices_file):
        if indices_file is None:
            outdir = os.path.join("output", "sta")
        else:
            outdir = os.path.join(
                "output", os.path.splitext(os.path.basename(indices_file))[0]
            )

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

    def test_abinit_het(self, star_particles, indices_file):
        if indices_file is None:
            outdir = os.path.join("output", "sta")
        else:
            outdir = os.path.join(
                "output", os.path.splitext(os.path.basename(indices_file))[0]
            )

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

    def test_train_model(self, indices_file):
        """Train the initial heterogeneous model."""
        outdir = "output_rand" if "rand" in indices_file else "output_first"
        shutil.rmtree(outdir, ignore_errors=True)

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

    def test_analyze(self, indices_file):
        """Produce standard analyses for a particular epoch."""
        outdir = "output_rand" if "rand" in indices_file else "output_first"

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
                "20",  # Number of kmeans samples to generate
                "--vol-start-index",
                "1",
            ]
        )
        analyze.main(args)
        assert os.path.exists(os.path.join(outdir, "analyze.2"))

    def test_notebooks(self, indices_file):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        outdir = "output_rand" if "rand" in indices_file else "output_first"

        assert os.path.exists(
            os.path.join(outdir, "analyze.2")
        ), "Upstream tests have failed!"
        os.chdir(os.path.join(outdir, "analyze.2"))

        for ipynb in ["cryoDRGN_figures", "cryoDRGN_filtering", "cryoDRGN_viz"]:
            assert os.path.exists(f"{ipynb}.ipynb")
            with open(f"{ipynb}.ipynb") as ff:
                nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

            ep = ExecutePreprocessor(timeout=600, kernel_name="python3")
            ep.preprocess(nb_in)

        os.chdir(os.path.join("..", ".."))
        shutil.rmtree(outdir)
