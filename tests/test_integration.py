"""Running an experiment of training followed by downstream analyses."""

import pytest
import argparse
import os.path
import pickle
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

from cryodrgn.commands import (
    analyze,
    analyze_landscape,
    analyze_landscape_full,
    backproject_voxel,
    eval_images,
    eval_vol,
    graph_traversal,
    train_nn,
    train_vae,
    abinit_homo,
    abinit_het,
)
from cryodrgn.commands_utils import filter_star
from cryodrgn.dataset import TiltSeriesData
import cryodrgn.utils


@pytest.mark.parametrize("particles", ["hand"], indirect=True)
@pytest.mark.parametrize("poses", ["hand-poses"], indirect=True)
class TestFixedHetero:
    def test_train_model(self, outdir, particles, poses):
        """Train the initial heterogeneous model."""

        args = train_vae.add_args(argparse.ArgumentParser()).parse_args(
            [
                particles,
                "-o",
                outdir,
                "--lr",
                ".0001",
                "--num-epochs",
                "3",
                "--seed",
                "0",
                "--poses",
                poses,
                "--zdim",
                "10",
                "--pe-type",
                "gaussian",
                "--multigpu",
            ]
        )
        train_vae.main(args)

    def test_train_from_checkpoint(self, outdir, particles, poses):
        """Load a cached model and run for another epoch, now without --multigpu."""

        train_vae.main(
            train_vae.add_args(argparse.ArgumentParser()).parse_args(
                [
                    particles,
                    "-o",
                    outdir,
                    "--lr",
                    ".0001",
                    "--num-epochs",
                    "4",
                    "--seed",
                    "0",
                    "--poses",
                    poses,
                    "--zdim",
                    "10",
                    "--pe-type",
                    "gaussian",
                    "--load",
                    os.path.join(outdir, "weights.2.pkl"),
                ]
            )
        )

    def test_analyze(self, outdir, particles, poses):
        """Produce standard analyses for a particular epoch."""
        args = analyze.add_args(argparse.ArgumentParser()).parse_args([outdir, "2"])
        analyze.main(args)
        assert os.path.exists(os.path.join(outdir, "analyze.2"))

    @pytest.mark.parametrize("nb_lbl", ["cryoDRGN_figures"])
    def test_notebooks(self, outdir, particles, poses, nb_lbl):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        orig_cwd = os.path.abspath(os.getcwd())
        os.chdir(os.path.join(outdir, "analyze.2"))
        assert os.path.exists(f"{nb_lbl}.ipynb")

        with open(f"{nb_lbl}.ipynb") as ff:
            nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

        ExecutePreprocessor(timeout=600, kernel_name="python3").preprocess(nb_in)
        os.chdir(orig_cwd)

    def test_landscape(self, outdir, particles, poses):
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

    def test_landscape_full(self, outdir, particles, poses):
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

    def test_graph_traversal(self, outdir, particles, poses):
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

    def test_eval_volume(self, outdir, particles, poses):
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

    def test_eval_images(self, outdir, particles, poses):
        args = eval_images.add_args(argparse.ArgumentParser()).parse_args(
            [
                particles,
                os.path.join(outdir, "weights.3.pkl"),
                "--config",
                os.path.join(outdir, "config.yaml"),
                "-o",
                os.path.join(outdir, "out_eval_images_losses.pkl"),
                "--out-z",
                os.path.join(outdir, "out_eval_images_z.pkl"),
                "--poses",
                poses,
                "--log-interval",
                "1",
                "--verbose",
            ]
        )
        eval_images.main(args)


@pytest.mark.parametrize("particles", ["toy.mrcs"], indirect=True)
@pytest.mark.parametrize("ctf", ["CTF-Test"], indirect=True)
@pytest.mark.parametrize("indices", ["random-100"], indirect=True)
class TestAbinitHetero:
    def test_train_model(self, outdir, particles, ctf, indices):
        """Train the initial heterogeneous model."""

        args = abinit_het.add_args(argparse.ArgumentParser()).parse_args(
            [
                particles,
                "--ctf",
                ctf,
                "--ind",
                indices,
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

    def test_analyze(self, outdir, particles, ctf, indices):
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
    def test_notebooks(self, outdir, particles, ctf, indices, nb_lbl):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        orig_cwd = os.path.abspath(os.getcwd())
        os.chdir(os.path.join(outdir, "analyze.1"))
        assert os.path.exists(f"{nb_lbl}.ipynb"), "Upstream tests have failed!"

        with open(f"{nb_lbl}.ipynb") as ff:
            nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

        ExecutePreprocessor(timeout=600, kernel_name="python3").preprocess(nb_in)
        os.chdir(orig_cwd)


@pytest.mark.parametrize("particles", ["tilts.star"], indirect=True)
@pytest.mark.parametrize("indices", [None, "just-4"], indirect=True)
@pytest.mark.parametrize("poses", ["tilt-poses"], indirect=True)
@pytest.mark.parametrize("ctf", ["CTF-Tilt"], indirect=True)
@pytest.mark.parametrize("datadir", ["default-datadir"], indirect=True)
class TestStarFixedHomo:
    """Run reconstructions using particles from a .star file as input."""

    def test_train_model(self, outdir, particles, indices, poses, ctf, datadir):
        args = [
            particles,
            "--datadir",
            datadir,
            "--poses",
            poses,
            "--ctf",
            ctf,
            "-o",
            outdir,
            "--dim",
            "32",
        ]
        if indices is not None:
            args += ["--ind", indices]

        args = train_nn.add_args(argparse.ArgumentParser()).parse_args(args)
        train_nn.main(args)


@pytest.mark.parametrize("particles", ["tilts.star"], indirect=True)
@pytest.mark.parametrize("indices", ["just-4", "just-5"], indirect=True)
@pytest.mark.parametrize("poses", ["tilt-poses"], indirect=True)
@pytest.mark.parametrize("ctf", ["CTF-Tilt"], indirect=True)
@pytest.mark.parametrize("datadir", ["default-datadir"], indirect=True)
class TestStarFixedHetero:
    """Run heterogeneous reconstruction using tilt series from a .star file and poses.

    We use two sets of indices, one that produces a tilt series with all particles
    having the same number of tilts and another that produces a ragged tilt-series.
    """

    def test_train_model(self, outdir, particles, indices, poses, ctf, datadir):
        args = [
            particles,
            "--datadir",
            datadir,
            "--encode-mode",
            "tilt",
            "--poses",
            poses,
            "--ctf",
            ctf,
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
        if indices is not None:
            args += ["--ind", indices]

        args = train_vae.add_args(argparse.ArgumentParser()).parse_args(args)
        train_vae.main(args)

    def test_filter_command(self, outdir, particles, indices, poses, ctf, datadir):
        # filter the tilt-series particles
        args = [
            particles,
            "--ind",
            indices,
            "--et",
            "-o",
            os.path.join(outdir, "filtered_sta_testing_bin8.star"),
        ]
        parser = argparse.ArgumentParser()
        filter_star.add_args(parser)
        filter_star.main(parser.parse_args(args))

        # need to filter poses and CTFs manually due to tilt indices
        pt, tp = TiltSeriesData.parse_particle_tilt(particles)
        ind = cryodrgn.utils.load_pkl(indices)
        new_ind = ind[:3]
        tilt_ind = TiltSeriesData.particles_to_tilts(pt, ind)

        rot, trans = cryodrgn.utils.load_pkl(poses)
        rot, trans = rot[tilt_ind], trans[tilt_ind]
        cryodrgn.utils.save_pkl(
            (rot, trans), os.path.join(outdir, "filtered_sta_pose.pkl")
        )

        ctf_mat = cryodrgn.utils.load_pkl(ctf)[tilt_ind]
        cryodrgn.utils.save_pkl(ctf_mat, os.path.join(outdir, "filtered_sta_ctf.pkl"))
        cryodrgn.utils.save_pkl(new_ind, os.path.join(outdir, "filtered_ind.pkl"))

        args = [
            os.path.join(outdir, "filtered_sta_testing_bin8.star"),
            "--datadir",
            datadir,
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
        train_vae.main(train_vae.add_args(argparse.ArgumentParser()).parse_args(args))

    def test_analyze(self, outdir, particles, indices, poses, ctf, datadir):
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

    @pytest.mark.parametrize(
        "new_indices_file",
        [None, "filtered_ind.pkl"],
        ids=("no-new-ind", "new-ind"),
    )
    def test_backproject(
        self, outdir, particles, indices, poses, ctf, datadir, new_indices_file
    ):
        """Run backprojection using the given particles."""
        args = [
            os.path.join(outdir, "filtered_sta_testing_bin8.star"),
            "--datadir",
            datadir,
            "--tilt",
            "--poses",
            os.path.join(outdir, "filtered_sta_pose.pkl"),
            "--ctf",
            os.path.join(outdir, "filtered_sta_ctf.pkl"),
            "-o",
            os.path.join(outdir, "filtered.mrc"),
            "-d",
            "2.93",
        ]

        if new_indices_file is not None:
            args += ["--ind", os.path.join(outdir, new_indices_file)]
        parser = argparse.ArgumentParser()
        backproject_voxel.add_args(parser)

        backproject_voxel.main(parser.parse_args(args))
        assert os.path.exists(os.path.join(outdir, "filtered.mrc"))
        os.remove(os.path.join(outdir, "filtered.mrc"))

    @pytest.mark.parametrize("nb_lbl", ["cryoDRGN_figures", "cryoDRGN_ET_viz"])
    def test_notebooks(self, outdir, particles, indices, poses, ctf, datadir, nb_lbl):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        orig_cwd = os.path.abspath(os.getcwd())
        os.chdir(os.path.join(outdir, "analyze.4"))
        assert os.path.exists(f"{nb_lbl}.ipynb"), "Upstream tests have failed!"

        with open(f"{nb_lbl}.ipynb") as ff:
            nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

        ExecutePreprocessor(timeout=600, kernel_name="python3").preprocess(nb_in)
        os.chdir(orig_cwd)

    def test_refiltering(self, outdir, particles, indices, poses, ctf, datadir):
        """Use particle index creating during analysis."""
        orig_cwd = os.path.abspath(os.getcwd())
        os.chdir(os.path.join(outdir, "analyze.4"))
        assert os.path.exists("tmp_ind_selected.pkl"), "Upstream tests have failed!"

        with open("tmp_ind_selected.pkl", "rb") as f:
            indices = pickle.load(f)

        new_indices = indices[:3]
        with open("tmp_ind_selected.pkl", "wb") as f:
            pickle.dump(new_indices, f)

        args = [
            particles,
            "--datadir",
            datadir,
            "--encode-mode",
            "tilt",
            "--poses",
            poses,
            "--ctf",
            ctf,
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

        args = train_vae.add_args(argparse.ArgumentParser()).parse_args(args)
        train_vae.main(args)
        os.chdir(orig_cwd)


@pytest.mark.parametrize("particles", ["tilts.star"], indirect=True)
@pytest.mark.parametrize("indices", [None, "just-4"], indirect=True)
@pytest.mark.parametrize("ctf", ["CTF-Tilt"], indirect=True)
@pytest.mark.parametrize("datadir", ["default-datadir"], indirect=True)
class TestStarAbinitHomo:
    def test_train_model(self, outdir, particles, indices, ctf, datadir):
        args = [
            particles,
            "--datadir",
            datadir,
            "--ctf",
            ctf,
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
        if indices is not None:
            args += ["--ind", indices]

        args = abinit_homo.add_args(argparse.ArgumentParser()).parse_args(args)
        abinit_homo.main(args)


@pytest.mark.parametrize("particles", ["tilts.star"], indirect=True)
@pytest.mark.parametrize("indices", [None, "just-4"], indirect=True)
@pytest.mark.parametrize("ctf", ["CTF-Tilt"], indirect=True)
@pytest.mark.parametrize("datadir", ["default-datadir"], indirect=True)
class TestStarAbinitHetero:
    def test_train_model(self, outdir, particles, indices, ctf, datadir):
        args = [
            particles,
            "--datadir",
            datadir,
            "--ctf",
            ctf,
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
        if indices is not None:
            args += ["--ind", indices]

        args = abinit_het.add_args(argparse.ArgumentParser()).parse_args(args)
        abinit_het.main(args)


@pytest.mark.parametrize("particles", ["toy.mrcs"], indirect=True)
@pytest.mark.parametrize("indices", ["first-100", "random-100"], indirect=True)
@pytest.mark.parametrize("poses", ["toy-poses"], indirect=True)
@pytest.mark.parametrize("ctf", ["CTF-Test"], indirect=True)
class TestIterativeFiltering:
    def test_train_model(self, outdir, particles, indices, poses, ctf):
        """Train the initial heterogeneous model."""
        args = train_vae.add_args(argparse.ArgumentParser()).parse_args(
            [
                particles,
                "-o",
                outdir,
                "--ctf",
                ctf,
                "--ind",
                indices,
                "--num-epochs",
                "3",
                "--poses",
                poses,
                "--zdim",
                "4",
                "--pe-type",
                "gaussian",
                "--no-amp",
            ]
        )
        train_vae.main(args)

    def test_analyze(self, outdir, particles, indices, poses, ctf):
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
    def test_notebooks(self, outdir, particles, indices, poses, ctf, nb_lbl):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        orig_cwd = os.path.abspath(os.getcwd())
        os.chdir(os.path.join(outdir, "analyze.2"))
        assert os.path.exists(f"{nb_lbl}.ipynb"), "Upstream tests have failed!"

        with open(f"{nb_lbl}.ipynb") as ff:
            nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

        ExecutePreprocessor(timeout=600, kernel_name="python3").preprocess(nb_in)
        os.chdir(orig_cwd)
