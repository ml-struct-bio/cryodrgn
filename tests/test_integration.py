"""Integration tests of ab initio reconstruction and downstream analyses.

Note that the training done here has unrealistically low parameter values to allow the
tests to run in any environment in a reasonable amount of time with or without GPUs.

"""
import pytest
import os
import argparse
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
from cryodrgn.commands import analyze, train_vae
import logging


@pytest.fixture(autouse=True)
def no_logging(caplog):
    caplog.set_level(logging.CRITICAL)


@pytest.mark.parametrize("particles", ["toy.mrcs"], indirect=True)
@pytest.mark.parametrize("poses", ["toy-poses"], indirect=True)
@pytest.mark.parametrize("ctf", ["CTF-Test"], indirect=True)
class TestNotebookFiltering:
    def test_train_model(self, outdir, particles, poses, ctf):
        """Train the initial heterogeneous model."""
        args = train_vae.add_args(argparse.ArgumentParser()).parse_args(
            [
                particles,
                "-o",
                outdir,
                "--ctf",
                ctf,
                "--num-epochs",
                "10",
                "--poses",
                poses,
                "--zdim",
                "4",
                "--tdim",
                "64",
                "--enc-dim",
                "64",
                "--dec-dim",
                "64",
                "--pe-type",
                "gaussian",
                "--no-amp",
            ]
        )
        train_vae.main(args)

    def test_analyze(self, outdir, particles, poses, ctf):
        """Produce standard analyses for a particular epoch."""
        assert os.path.exists(
            os.path.join(outdir, "weights.9.pkl")
        ), "Upstream tests have failed!"

        args = analyze.add_args(argparse.ArgumentParser()).parse_args(
            [
                outdir,
                "9",  # Epoch number to analyze - 0-indexed
                "--pc",
                "2",  # Number of principal component traversals to generate
                "--ksample",
                "10",  # Number of kmeans samples to generate
                "--vol-start-index",
                "1",
            ]
        )
        analyze.main(args)
        assert os.path.exists(os.path.join(outdir, "analyze.9"))

    @pytest.mark.parametrize(
        "nb_lbl", ["cryoDRGN_figures", "cryoDRGN_filtering", "cryoDRGN_viz"]
    )
    def test_notebooks(self, outdir, particles, poses, ctf, nb_lbl):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        orig_cwd = os.path.abspath(os.getcwd())
        os.chdir(os.path.join(outdir, "analyze.9"))
        assert os.path.exists(f"{nb_lbl}.ipynb"), "Upstream tests have failed!"

        with open(f"{nb_lbl}.ipynb") as ff:
            nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

        ExecutePreprocessor(timeout=3, kernel_name="python3").preprocess(nb_in)
        os.chdir(orig_cwd)

    def test_refiltering(self, outdir, particles, poses, ctf):
        ind_keep_fl = [fl for fl in os.listdir(outdir) if fl[:9] == "ind_keep."]
        assert len(ind_keep_fl) == 1
        assert int(ind_keep_fl[0].split(".")[1].split("_particles")[0]) >= 10
        ind_keep_fl = os.path.join(outdir, ind_keep_fl[0])

        args = train_vae.add_args(argparse.ArgumentParser()).parse_args(
            [
                particles,
                "-o",
                outdir,
                "--ctf",
                ctf,
                "--ind",
                ind_keep_fl,
                "--num-epochs",
                "5",
                "--poses",
                poses,
                "--zdim",
                "4",
                "--tdim",
                "64",
                "--enc-dim",
                "64",
                "--dec-dim",
                "64",
                "--pe-type",
                "gaussian",
                "--no-amp",
            ]
        )
        train_vae.main(args)
