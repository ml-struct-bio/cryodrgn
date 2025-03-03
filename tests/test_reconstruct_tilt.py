"""Running an experiment of training followed by downstream analyses."""

import pytest
import argparse
import os.path
import shutil
import pickle
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

from cryodrgn.commands import (
    analyze,
    backproject_voxel,
    train_vae,
    abinit_homo,
    abinit_het,
)
from cryodrgn.commands_utils import filter_star
from cryodrgn.dataset import TiltSeriesData
from cryodrgn import utils


@pytest.mark.parametrize("particles", ["tilts.star"], indirect=True)
@pytest.mark.parametrize("indices", ["just-4", "just-5"], indirect=True)
@pytest.mark.parametrize("poses", ["tilt-poses"], indirect=True)
@pytest.mark.parametrize("ctf", ["CTF-Tilt"], indirect=True)
@pytest.mark.parametrize("datadir", ["default-datadir"], indirect=True)
class TestTiltFixedHetero:
    """Run heterogeneous reconstruction using tilt series from a .star file and poses.

    We use two sets of indices, one that produces a tilt series with all particles
    having the same number of tilts and another that produces a ragged tilt-series.
    """

    def get_outdir(self, tmpdir_factory, particles, poses, ctf, indices, datadir):
        dirname = os.path.join(
            "TiltFixedHetero",
            particles.label,
            poses.label,
            ctf.label,
            indices.label,
            datadir.label,
        )
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)

        return odir

    def test_train_model(self, tmpdir_factory, particles, indices, poses, ctf, datadir):
        outdir = self.get_outdir(
            tmpdir_factory, particles, poses, ctf, indices, datadir
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
        ]
        if indices.path is not None:
            args += ["--ind", indices.path]

        parser = argparse.ArgumentParser()
        train_vae.add_args(parser)
        train_vae.main(parser.parse_args(args))

    def test_filter_command(
        self, tmpdir_factory, particles, indices, poses, ctf, datadir
    ):
        outdir = self.get_outdir(
            tmpdir_factory, particles, poses, ctf, indices, datadir
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
        ]
        parser = argparse.ArgumentParser()
        train_vae.add_args(parser)
        train_vae.main(parser.parse_args(args))

    def test_analyze(self, tmpdir_factory, particles, indices, poses, ctf, datadir):
        """Produce standard analyses for a particular epoch."""
        outdir = self.get_outdir(
            tmpdir_factory, particles, poses, ctf, indices, datadir
        )

        parser = argparse.ArgumentParser()
        analyze.add_args(parser)
        analyze.main(
            parser.parse_args([outdir, "--epoch", "5", "--pc", "3", "--ksample", "2"])
        )
        assert os.path.exists(os.path.join(outdir, "analyze.5"))

    @pytest.mark.parametrize(
        "new_indices_file",
        [None, "filtered_ind.pkl"],
        ids=("no-new-ind", "new-ind"),
    )
    def test_backproject(
        self, tmpdir_factory, particles, indices, poses, ctf, datadir, new_indices_file
    ):
        """Run backprojection using the given particles."""
        outdir = self.get_outdir(
            tmpdir_factory, particles, poses, ctf, indices, datadir
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

    @pytest.mark.parametrize("nb_lbl", ["cryoDRGN_figures", "ET-viz"])
    def test_notebooks(
        self, tmpdir_factory, particles, indices, poses, ctf, datadir, nb_lbl
    ):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        outdir = self.get_outdir(
            tmpdir_factory, particles, poses, ctf, indices, datadir
        )
        orig_cwd = os.path.abspath(os.getcwd())
        os.chdir(os.path.join(outdir, "analyze.5"))
        assert os.path.exists(f"{nb_lbl}.ipynb"), "Upstream tests have failed!"

        with open(f"{nb_lbl}.ipynb") as ff:
            nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

        ExecutePreprocessor(timeout=600, kernel_name="python3").preprocess(nb_in)
        os.chdir(orig_cwd)

    def test_refiltering(self, tmpdir_factory, particles, indices, poses, ctf, datadir):
        """Use particle index creating during analysis."""
        outdir = self.get_outdir(
            tmpdir_factory, particles, poses, ctf, indices, datadir
        )
        orig_cwd = os.path.abspath(os.getcwd())
        os.chdir(os.path.join(outdir, "analyze.5"))
        assert os.path.exists("tmp_ind_selected.pkl"), "Upstream tests have failed!"

        with open("tmp_ind_selected.pkl", "rb") as f:
            indices = pickle.load(f)

        new_indices = indices[:3]
        with open("tmp_ind_selected.pkl", "wb") as f:
            pickle.dump(new_indices, f)

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
            "tmp_ind_selected.pkl",
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

        parser = argparse.ArgumentParser()
        train_vae.add_args(parser)
        train_vae.main(parser.parse_args(args))
        os.chdir(orig_cwd)

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
        abinit_homo.main(parser.parse_args(args))


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
        ]
        if indices.path is not None:
            args += ["--ind", indices.path]

        parser = argparse.ArgumentParser()
        abinit_het.add_args(parser)
        abinit_het.main(parser.parse_args(args))
