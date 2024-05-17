"""Integration tests of ab initio reconstruction and downstream analyses.

Note that the training done here has unrealistically low parameter values to allow the
tests to run in any environment in a reasonable amount of time with or without GPUs.

"""
import pytest
import os
import shutil
import argparse
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
from cryodrgn.commands import analyze, train_vae


@pytest.mark.parametrize("particles", ["toy.mrcs"], indirect=True)
@pytest.mark.parametrize("poses", ["toy-poses"], indirect=True)
@pytest.mark.parametrize("ctf", ["CTF-Test"], indirect=True)
@pytest.mark.parametrize("indices", [None, "first-100", "random-100"], indirect=True)
class TestIterativeFiltering:
    def get_outdir(self, tmpdir_factory, particles, indices, poses, ctf):
        dirname = os.path.join(
            "IterativeFiltering", particles.label, indices.label, poses.label, ctf.label
        )
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)

        return odir

    def test_train_model(self, tmpdir_factory, particles, poses, ctf, indices):
        """Train the initial heterogeneous model without any manual filters."""
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        args = train_vae.add_args(argparse.ArgumentParser()).parse_args(
            [
                particles.file,
                "-o",
                outdir,
                "--ctf",
                ctf.file,
                "--num-epochs",
                "10",
                "--poses",
                poses.file,
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

    def test_analyze(self, tmpdir_factory, particles, poses, ctf, indices):
        """Produce standard analyses for the final epoch."""
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
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
    def test_notebooks(self, tmpdir_factory, particles, poses, ctf, indices, nb_lbl):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        orig_cwd = os.path.abspath(os.getcwd())
        os.chdir(os.path.join(outdir, "analyze.9"))
        assert os.path.exists(f"{nb_lbl}.ipynb"), "Upstream tests have failed!"

        with open(f"{nb_lbl}.ipynb") as ff:
            nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

        ExecutePreprocessor(timeout=60, kernel_name="python3").preprocess(nb_in)
        os.chdir(orig_cwd)

    def test_refiltering(self, tmpdir_factory, particles, poses, ctf, indices):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        ind_keep_fl = [fl for fl in os.listdir(outdir) if fl[:9] == "ind_keep."][0]
        if int(ind_keep_fl.split(".")[1].split("_particles")[0]) < 50:
            ind_keep_fl = [fl for fl in os.listdir(outdir) if fl[:8] == "ind_bad."][0]

        ind_keep_fl = os.path.join(outdir, ind_keep_fl)
        args = train_vae.add_args(argparse.ArgumentParser()).parse_args(
            [
                particles.file,
                "-o",
                outdir,
                "--ctf",
                ctf.file,
                "--ind",
                ind_keep_fl,
                "--num-epochs",
                "5",
                "--poses",
                poses.file,
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

        shutil.rmtree(outdir)
