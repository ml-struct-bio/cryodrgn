"""Running an experiment of training followed by downstream analyses."""

import pytest
import argparse
import os.path
import shutil
import pickle
import numpy as np
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

from cryodrgn.commands import (
    analyze,
    backproject_voxel,
    train_vae,
    abinit_homo,
    abinit_het,
    filter,
)
from cryodrgn.commands_utils import filter_star
from cryodrgn.dataset import TiltSeriesData
from cryodrgn import utils


@pytest.mark.parametrize("particles", ["tilts.star"], indirect=True)
@pytest.mark.parametrize("indices", ["just-4"], indirect=True)
@pytest.mark.parametrize("poses", ["tilt-poses"], indirect=True)
@pytest.mark.parametrize("ctf", ["CTF-Tilt"], indirect=True)
@pytest.mark.parametrize("datadir", ["default-datadir"], indirect=True)
@pytest.mark.parametrize("ntilts", [None, 30])
class TestTiltFixedHetero:
    """Run heterogeneous reconstruction using tilt series from a .star file and poses.

    We use two sets of indices, one that produces a tilt series with all particles
    having the same number of tilts and another that produces a ragged tilt-series.
    """

    def get_outdir(
        self, tmpdir_factory, particles, poses, ctf, indices, datadir, ntilts
    ):
        dirname = os.path.join(
            "TiltFixedHetero",
            particles.label,
            poses.label,
            ctf.label,
            indices.label,
            datadir.label,
            f"ntilts.{ntilts}" if ntilts is not None else "ntilts.def",
        )
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)

        return odir

    def test_train_model(
        self, tmpdir_factory, particles, indices, poses, ctf, datadir, ntilts
    ):
        outdir = self.get_outdir(
            tmpdir_factory, particles, poses, ctf, indices, datadir, ntilts
        )
        args = [
            particles.path,
            "--datadir",
            datadir.path,
            "--encode-mode",
            "tilt",
            "--poses",
            poses.path,
            "--ctf",
            ctf.path,
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
            "--no-analysis",
        ]
        if indices.path is not None:
            args += ["--ind", indices.path]
        if ntilts is not None:
            args += ["--ntilts", str(ntilts)]

        parser = argparse.ArgumentParser()
        train_vae.add_args(parser)
        args = parser.parse_args(args)
        train_vae.main(args)

    def test_filter_command(
        self, tmpdir_factory, particles, indices, poses, ctf, datadir, ntilts
    ):
        outdir = self.get_outdir(
            tmpdir_factory, particles, poses, ctf, indices, datadir, ntilts
        )

        # filter the tilt-series particles
        args = [
            particles.path,
            "--et",
            "-o",
            os.path.join(outdir, "filtered_sta_testing_bin8.star"),
        ]
        if indices.path is not None:
            args += ["--ind", indices.path]
        parser = argparse.ArgumentParser()
        filter_star.add_args(parser)
        filter_star.main(parser.parse_args(args))

        # need to filter poses and CTFs manually due to tilt indices
        pt, tp = TiltSeriesData.parse_particle_tilt(particles.path)
        ind = utils.load_pkl(indices.path)
        new_ind = ind[:3]
        tilt_ind = TiltSeriesData.particles_to_tilts(pt, ind)

        rot, trans = utils.load_pkl(poses.path)
        rot, trans = rot[tilt_ind], trans[tilt_ind]
        utils.save_pkl((rot, trans), os.path.join(outdir, "filtered_sta_pose.pkl"))
        ctf_mat = utils.load_pkl(ctf.path)[tilt_ind]
        utils.save_pkl(ctf_mat, os.path.join(outdir, "filtered_sta_ctf.pkl"))
        utils.save_pkl(new_ind, os.path.join(outdir, "filtered_ind.pkl"))
        args = [
            os.path.join(outdir, "filtered_sta_testing_bin8.star"),
            "--datadir",
            datadir.path,
            "--encode-mode",
            "tilt",
            "--ntilts",
            "5",
            "--poses",
            os.path.join(outdir, "filtered_sta_pose.pkl"),
            "--ctf",
            os.path.join(outdir, "filtered_sta_ctf.pkl"),
            "--ind",
            os.path.join(outdir, "filtered_ind.pkl"),
            "--num-epochs",
            "5",
            "--zdim",
            "4",
            "-o",
            os.path.join(outdir, "filtered"),
            "--tdim",
            "16",
            "--enc-dim",
            "16",
            "--dec-dim",
            "16",
            "--no-analysis",
        ]
        parser = argparse.ArgumentParser()
        train_vae.add_args(parser)
        train_vae.main(parser.parse_args(args))

    def test_analyze(
        self, tmpdir_factory, particles, indices, poses, ctf, datadir, ntilts
    ):
        """Produce standard analyses for a particular epoch."""
        outdir = self.get_outdir(
            tmpdir_factory, particles, poses, ctf, indices, datadir, ntilts
        )

        parser = argparse.ArgumentParser()
        analyze.add_args(parser)
        analyze.main(
            parser.parse_args(
                [
                    outdir,
                    "5",  # Epoch number to analyze - 1-indexed
                    "--pc",
                    "3",  # Number of principal component traversals to generate
                    "--ksample",
                    "2",  # Number of kmeans samples to generate
                ]
            )
        )
        assert os.path.exists(os.path.join(outdir, "analyze.5"))

    @pytest.mark.parametrize(
        "new_indices_file",
        [None, "filtered_ind.pkl"],
        ids=("no-new-ind", "new-ind"),
    )
    def test_backproject(
        self,
        tmpdir_factory,
        particles,
        indices,
        poses,
        ctf,
        datadir,
        ntilts,
        new_indices_file,
    ):
        """Run backprojection using the given particles."""
        outdir = self.get_outdir(
            tmpdir_factory, particles, poses, ctf, indices, datadir, ntilts
        )
        args = [
            os.path.join(outdir, "filtered_sta_testing_bin8.star"),
            "--datadir",
            datadir.path,
            "--tilt",
            "--poses",
            os.path.join(outdir, "filtered_sta_pose.pkl"),
            "--ctf",
            os.path.join(outdir, "filtered_sta_ctf.pkl"),
            "-o",
            os.path.join(outdir, "filtered"),
            "-d",
            "2.93",
            "--no-half-maps",
        ]
        if new_indices_file is not None:
            args += ["--ind", os.path.join(outdir, new_indices_file)]

        parser = argparse.ArgumentParser()
        backproject_voxel.add_args(parser)
        backproject_voxel.main(parser.parse_args(args))
        assert os.path.exists(os.path.join(outdir, "filtered", "backproject.mrc"))
        shutil.rmtree(os.path.join(outdir, "filtered"))

    @pytest.mark.parametrize("nb_lbl", ["cryoDRGN_figures", "cryoDRGN_ET_viz"])
    def test_notebooks(
        self, tmpdir_factory, particles, indices, poses, ctf, datadir, ntilts, nb_lbl
    ):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        outdir = self.get_outdir(
            tmpdir_factory, particles, poses, ctf, indices, datadir, ntilts
        )
        orig_cwd = os.path.abspath(os.getcwd())
        os.chdir(os.path.join(outdir, "analyze.5"))
        assert os.path.exists(f"{nb_lbl}.ipynb"), "Upstream tests have failed!"

        with open(f"{nb_lbl}.ipynb") as ff:
            nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

        ExecutePreprocessor(timeout=600, kernel_name="python3").preprocess(nb_in)
        os.chdir(orig_cwd)

    @pytest.mark.parametrize("plot_ind", [False, True])
    def test_interactive_filtering(
        self, tmpdir_factory, particles, indices, poses, ctf, datadir, ntilts, plot_ind
    ):
        """Launch interface for filtering particles using model covariates."""
        outdir = self.get_outdir(
            tmpdir_factory, particles, poses, ctf, indices, datadir, ntilts
        )
        parser = argparse.ArgumentParser()
        filter.add_args(parser)
        args = [outdir] + ["--epoch", "5", "--force"]
        if plot_ind:
            ind_fl = os.path.join(outdir, "analyze.5", "tmp_ind_test.pkl")
            with open(ind_fl, "wb") as f:
                pickle.dump(np.array([1, 2]), f)
            args += ["--plot-inds", ind_fl]

        filter.main(parser.parse_args(args))
        assert os.path.exists("indices.pkl")
        assert os.path.exists("indices_inverse.pkl")

    @pytest.mark.parametrize("newinds", ["indices.pkl", "indices_inverse.pkl"])
    def test_refiltering(
        self, tmpdir_factory, particles, indices, poses, ctf, datadir, ntilts, newinds
    ):
        """Use particle index creating during analysis."""
        outdir = self.get_outdir(
            tmpdir_factory, particles, poses, ctf, indices, datadir, ntilts
        )
        shutil.rmtree(outdir)

        args = [
            particles.path,
            "--datadir",
            datadir.path,
            "--encode-mode",
            "tilt",
            "--poses",
            poses.path,
            "--ctf",
            ctf.path,
            "--ind",
            newinds,
            "--num-epochs",
            "3",
            "--zdim",
            "8",
            "-o",
            outdir,
            "--tdim",
            "16",
            "--enc-dim",
            "16",
            "--dec-dim",
            "16",
            "--no-analysis",
            "--checkpoint",
            "1",
        ]
        parser = argparse.ArgumentParser()
        train_vae.add_args(parser)
        train_vae.main(parser.parse_args(args))
        assert os.path.exists(os.path.join(outdir, "weights.3.pkl"))
        assert not os.path.exists(os.path.join(outdir, "weights.4.pkl"))

        shutil.rmtree(outdir)


@pytest.mark.parametrize("particles", ["tilts.star"], indirect=True)
@pytest.mark.parametrize("indices", [None, "just-4"], indirect=True)
@pytest.mark.parametrize("ctf", ["CTF-Tilt"], indirect=True)
@pytest.mark.parametrize("datadir", ["default-datadir"], indirect=True)
class TestTiltAbinitHomo:
    def test_train_model(self, tmpdir, particles, indices, ctf, datadir):
        args = [
            particles.path,
            "--datadir",
            datadir.path,
            "--ctf",
            ctf.path,
            "-o",
            str(tmpdir),
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
        if indices.path is not None:
            args += ["--ind", indices.path]

        parser = argparse.ArgumentParser()
        abinit_homo.add_args(parser)
        args = parser.parse_args(args)
        abinit_homo.main(args)


@pytest.mark.parametrize("particles", ["tilts.star"], indirect=True)
@pytest.mark.parametrize("indices", [None, "just-4"], indirect=True)
@pytest.mark.parametrize("ctf", ["CTF-Tilt"], indirect=True)
@pytest.mark.parametrize("datadir", ["default-datadir"], indirect=True)
class TestTiltAbinitHetero:
    def test_train_model(self, tmpdir, particles, indices, ctf, datadir):
        args = [
            particles.path,
            "--datadir",
            datadir.path,
            "--ctf",
            ctf.path,
            "--zdim",
            "8",
            "-o",
            str(tmpdir),
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
            "--no-analysis",
        ]
        if indices.path is not None:
            args += ["--ind", indices.path]

        parser = argparse.ArgumentParser()
        abinit_het.add_args(parser)
        args = parser.parse_args(args)
        abinit_het.main(args)
