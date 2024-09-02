"""Running an experiment of training followed by downstream analyses."""

import pytest
import argparse
import os.path
import shutil
import pickle
import random
import nbformat
from nbclient.exceptions import CellExecutionError
from nbconvert.preprocessors import ExecutePreprocessor
import numpy as np

from cryodrgn.commands import (
    analyze,
    analyze_landscape,
    analyze_landscape_full,
    backproject_voxel,
    direct_traversal,
    eval_images,
    eval_vol,
    graph_traversal,
    train_nn,
    train_vae,
    abinit_homo,
    abinit_het,
)
from cryodrgn.commands_utils import clean, filter_star, plot_classes
from cryodrgn.source import ImageSource
from cryodrgn.dataset import TiltSeriesData
from cryodrgn import utils


@pytest.mark.parametrize(
    "particles", ["toy.mrcs", "toy.txt", "toy.star"], indirect=True
)
@pytest.mark.parametrize("poses", ["toy-poses"], indirect=True)
@pytest.mark.parametrize("indices", [None, "random-100"], indirect=True)
class TestFixedHetero:
    def get_outdir(self, tmpdir_factory, particles, poses, ctf, indices):
        dirname = os.path.join(
            "FixedHetero", particles.label, poses.label, ctf.label, indices.label
        )
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)

        return odir

    @pytest.mark.parametrize("ctf", [None, "CTF-Test"], indirect=True)
    def test_train_model(self, tmpdir_factory, particles, poses, ctf, indices):
        """Train the initial heterogeneous model."""
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        args = [
            particles.path,
            "-o",
            outdir,
            "--lr",
            ".0001",
            "--num-epochs",
            "3",
            "--seed",
            "0",
            "--poses",
            poses.path,
            "--zdim",
            "10",
            "--tdim",
            "8",
            "--enc-dim",
            "8",
            "--dec-dim",
            "8",
            "--pe-type",
            "gaussian",
            "--multigpu",
        ]
        if ctf.path is not None:
            args += ["--ctf", ctf.path]
        if indices.path is not None:
            args += ["--ind", indices.path]

        train_vae.main(train_vae.add_args(argparse.ArgumentParser()).parse_args(args))

    @pytest.mark.parametrize("ctf", [None, "CTF-Test"], indirect=True)
    def test_train_from_checkpoint(
        self, tmpdir_factory, particles, poses, ctf, indices
    ):
        """Load a cached model and run for another epoch, now without --multigpu."""
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        args = [
            particles.path,
            "-o",
            outdir,
            "--lr",
            ".0001",
            "--num-epochs",
            "4",
            "--seed",
            "0",
            "--poses",
            poses.path,
            "--zdim",
            "10",
            "--tdim",
            "8",
            "--enc-dim",
            "8",
            "--dec-dim",
            "8",
            "--pe-type",
            "gaussian",
            "--load",
            os.path.join(outdir, "weights.2.pkl"),
        ]
        if ctf.path is not None:
            args += ["--ctf", ctf.path]
        if indices.path is not None:
            args += ["--ind", indices.path]

        train_vae.main(train_vae.add_args(argparse.ArgumentParser()).parse_args(args))

    @pytest.mark.parametrize(
        "ctf, epoch", [("CTF-Test", 2), (None, 3), ("CTF-Test", 3)], indirect=["ctf"]
    )
    def test_analyze(self, tmpdir_factory, particles, poses, ctf, indices, epoch):
        """Produce standard analyses for a particular epoch."""
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        args = analyze.add_args(argparse.ArgumentParser()).parse_args(
            [outdir, str(epoch)]
        )
        analyze.main(args)
        assert os.path.exists(os.path.join(outdir, f"analyze.{epoch}"))

    @pytest.mark.parametrize(
        "nb_lbl, ctf",
        [
            ("cryoDRGN_filtering", "CTF-Test"),
            ("cryoDRGN_figures", "CTF-Test"),
            ("cryoDRGN_viz", "CTF-Test"),
            pytest.param(
                "cryoDRGN_filtering",
                None,
                marks=pytest.mark.xfail(
                    raises=CellExecutionError, reason="need CTFs for filtering ntbook"
                ),
            ),
            ("cryoDRGN_figures", None),
        ],
        indirect=["ctf"],
    )
    def test_notebooks(self, tmpdir_factory, particles, poses, ctf, indices, nb_lbl):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        orig_cwd = os.path.abspath(os.getcwd())
        os.chdir(os.path.join(outdir, "analyze.3"))
        assert os.path.exists(f"{nb_lbl}.ipynb")

        with open(f"{nb_lbl}.ipynb") as ff:
            nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

        try:
            ExecutePreprocessor(timeout=600, kernel_name="python3").preprocess(nb_in)
        except CellExecutionError as e:
            os.chdir(orig_cwd)
            raise e

        os.chdir(orig_cwd)

    @pytest.mark.parametrize(
        "ctf, downsample_dim",
        [
            (None, "16"),
            ("CTF-Test", "16"),
            pytest.param(
                "CTF-Test",
                "64",
                marks=pytest.mark.xfail(
                    raises=AssertionError, reason="box size > resolution"
                ),
            ),
            pytest.param(
                "CTF-Test",
                None,
                marks=pytest.mark.xfail(
                    raises=AssertionError, reason="box size > resolution"
                ),
            ),
        ],
        indirect=["ctf"],
    )
    def test_landscape(
        self, tmpdir_factory, particles, poses, ctf, indices, downsample_dim
    ):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        args = [
            outdir,
            "3",  # Epoch number to analyze - 0-indexed
            "--sketch-size",
            "10",  # Number of volumes to generate for analysis
            "--pc-dim",
            "5",
            "--vol-start-index",
            "1",
        ]
        if downsample_dim:
            args += ["--downsample", downsample_dim]

        analyze_landscape.main(
            analyze_landscape.add_args(argparse.ArgumentParser()).parse_args(args)
        )

    @pytest.mark.parametrize(
        "ctf, downsample_dim",
        [
            (None, "16"),
            ("CTF-Test", "16"),
            pytest.param(
                "CTF-Test",
                "64",
                marks=pytest.mark.xfail(
                    raises=AssertionError, reason="box size > resolution"
                ),
            ),
            pytest.param(
                "CTF-Test",
                None,
                marks=pytest.mark.xfail(
                    raises=AssertionError, reason="box size > resolution"
                ),
            ),
        ],
        indirect=["ctf"],
    )
    def test_landscape_full(
        self, tmpdir_factory, particles, poses, ctf, indices, downsample_dim
    ):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        parser = argparse.ArgumentParser()
        args = [outdir, "3", "-N", "10"]
        if downsample_dim is not None:
            args += ["--downsample", downsample_dim]

        analyze_landscape_full.main(
            analyze_landscape_full.add_args(parser).parse_args(args)
        )

    @pytest.mark.parametrize(
        "ctf, seed, steps, points",
        [
            (None, 915, 5, None),
            ("CTF-Test", 321, 2, None),
            ("CTF-Test", 701, 3, 1),
            ("CTF-Test", 701, 3, 2),
            ("CTF-Test", 55, 3, None),
        ],
        indirect=["ctf"],
    )
    def test_direct_traversal(
        self, tmpdir_factory, particles, poses, ctf, indices, seed, steps, points
    ):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        random.seed(seed)
        anchors = [str(anchor) for anchor in random.sample(range(100), steps)]

        parser = argparse.ArgumentParser()
        direct_traversal.add_args(parser)
        args = [os.path.join(outdir, "z.3.pkl"), "--anchors"] + anchors
        if points is not None:
            args += ["-n", str(points)]

        direct_traversal.main(parser.parse_args(args))

    @pytest.mark.parametrize(
        "ctf, epoch, seed, steps",
        [
            (None, 3, 915, 5),
            ("CTF-Test", 2, 321, 2),
            ("CTF-Test", 3, 701, 5),
            ("CTF-Test", 3, 102, 10),
        ],
        indirect=["ctf"],
    )
    def test_graph_traversal(
        self, tmpdir_factory, particles, poses, ctf, indices, epoch, seed, steps
    ):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        random.seed(seed)
        anchors = [str(anchor) for anchor in random.sample(range(100), steps)]

        parser = argparse.ArgumentParser()
        graph_traversal.add_args(parser)
        args = parser.parse_args(
            [
                os.path.join(outdir, f"z.{epoch}.pkl"),
                "--anchors",
            ]
            + anchors
            + [
                "--outind",
                os.path.join(outdir, f"graph_traversal_path.{epoch}.txt"),
                "--outtxt",
                os.path.join(outdir, f"graph_traversal_zpath.{epoch}.txt"),
            ]
        )
        graph_traversal.main(args)

    @pytest.mark.parametrize("ctf", ["CTF-Test"], indirect=True)
    @pytest.mark.parametrize("epoch", [2, 3])
    def test_eval_volume(self, tmpdir_factory, particles, poses, ctf, indices, epoch):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        args = eval_vol.add_args(argparse.ArgumentParser()).parse_args(
            [
                os.path.join(outdir, f"weights.{epoch}.pkl"),
                "--config",
                os.path.join(outdir, "config.yaml"),
                "--zfile",
                os.path.join(outdir, f"graph_traversal_zpath.{epoch}.txt"),
                "-o",
                os.path.join(outdir, "eval_vols"),
            ]
        )
        eval_vol.main(args)

    @pytest.mark.parametrize("ctf", ["CTF-Test"], indirect=True)
    @pytest.mark.parametrize("epoch", [2, 3])
    def test_eval_images(self, tmpdir_factory, particles, poses, ctf, indices, epoch):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        args = eval_images.add_args(argparse.ArgumentParser()).parse_args(
            [
                particles.path,
                os.path.join(outdir, f"weights.{epoch}.pkl"),
                "--config",
                os.path.join(outdir, "config.yaml"),
                "-o",
                os.path.join(outdir, f"out_eval_images_losses.{epoch}.pkl"),
                "--out-z",
                os.path.join(outdir, f"out_eval_images_z.{epoch}.pkl"),
                "--poses",
                poses.path,
                "--log-interval",
                "1",
                "--verbose",
            ]
        )
        eval_images.main(args)

    @pytest.mark.parametrize("ctf", ["CTF-Test"], indirect=True)
    @pytest.mark.parametrize(
        "epoch, palette, plot_outdir",
        [
            (-1, "rocket", None),
            pytest.param(
                3,
                "Rocket",
                None,
                marks=pytest.mark.xfail(
                    raises=ValueError, reason="palette not available in seaborn!"
                ),
            ),
            (2, None, None),
            (2, None, "plots"),
        ],
    )
    def test_plot_classes(
        self,
        tmpdir_factory,
        particles,
        poses,
        ctf,
        indices,
        epoch,
        palette,
        plot_outdir,
    ):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        ind = utils.load_pkl(indices.path) if indices.path is not None else None
        particles_n = ImageSource.from_file(particles.path, indices=ind, lazy=True).n
        lbl_file = os.path.join(outdir, "plot-classes.pkl")
        if not os.path.exists(lbl_file):
            labels = np.array([0 if i % 3 == 0 else 1 for i in range(particles_n)])
            utils.save_pkl(labels, lbl_file)

        parser = argparse.ArgumentParser()
        plot_classes.add_args(parser)
        args = [outdir, str(epoch), "--labels", lbl_file]
        if palette is not None:
            args += ["--palette", palette]
        if plot_outdir is not None:
            args += ["--outdir", os.path.join(outdir, plot_outdir)]

        plot_classes.main(parser.parse_args(args))
        if plot_outdir is not None:
            use_outdir = os.path.join(outdir, plot_outdir)
        elif epoch == -1:
            use_outdir = os.path.join(outdir, "analyze")
        else:
            use_outdir = os.path.join(outdir, f"analyze.{epoch}")

        assert os.path.exists(os.path.join(use_outdir, "umap_kde_classes.png"))

    @pytest.mark.parametrize("ctf", [None, "CTF-Test"], indirect=True)
    def test_clean_all(self, tmpdir_factory, particles, poses, ctf, indices):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        parser = argparse.ArgumentParser()
        clean.add_args(parser)

        clean.main(parser.parse_args([os.path.relpath(outdir)]))
        shutil.rmtree(outdir)


@pytest.mark.parametrize(
    "particles, indices, ctf",
    [
        ("hand", None, None),
        ("hand", None, "CTF-Test.100"),
        ("toy.txt", "random-100", "CTF-Test"),
        ("toy.star", "first-100", "CTF-Test"),
        ("toy.star-13", None, None),
    ],
    indirect=True,
)
class TestAbinitHetero:
    def get_outdir(self, tmpdir_factory, particles, ctf, indices):
        dirname = os.path.join("AbinitHet", particles.label, ctf.label, indices.label)
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)

        return odir

    def test_train_model(self, tmpdir_factory, particles, ctf, indices):
        """Train the initial heterogeneous model."""
        outdir = self.get_outdir(tmpdir_factory, particles, indices, ctf)
        args = [
            particles.path,
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
        if ctf.path is not None:
            args += ["--ctf", ctf.path]
        if indices.path is not None:
            args += ["--ind", indices.path]

        abinit_het.main(abinit_het.add_args(argparse.ArgumentParser()).parse_args(args))

    def test_analyze(self, tmpdir_factory, particles, ctf, indices):
        """Produce standard analyses for a particular epoch."""
        outdir = self.get_outdir(tmpdir_factory, particles, indices, ctf)
        args = analyze.add_args(argparse.ArgumentParser()).parse_args(
            [
                outdir,
                "1",  # Epoch number to analyze - 0-indexed
                "--pc",
                "3",  # Number of principal component traversals to generate
                "--ksample",
                "10",  # Number of kmeans samples to generate
                "--vol-start-index",
                "1",
            ]
        )
        analyze.main(args)
        assert os.path.exists(os.path.join(outdir, "analyze.1"))

    @pytest.mark.parametrize("nb_lbl", ["cryoDRGN_figures"])
    def test_notebooks(self, tmpdir_factory, particles, ctf, indices, nb_lbl):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        outdir = self.get_outdir(tmpdir_factory, particles, indices, ctf)
        orig_cwd = os.path.abspath(os.getcwd())
        os.chdir(os.path.join(outdir, "analyze.1"))
        assert os.path.exists(f"{nb_lbl}.ipynb"), "Upstream tests have failed!"

        with open(f"{nb_lbl}.ipynb") as ff:
            nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

        ExecutePreprocessor(timeout=600, kernel_name="python3").preprocess(nb_in)
        os.chdir(orig_cwd)

    @pytest.mark.parametrize("epoch", [1, 2])
    def test_graph_traversal(self, tmpdir_factory, particles, ctf, indices, epoch):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, ctf)
        parser = argparse.ArgumentParser()
        graph_traversal.add_args(parser)
        args = parser.parse_args(
            [
                os.path.join(outdir, f"z.{epoch}.pkl"),
                "--anchors",
                "5",
                "0",
                "2",
                "10",
                "--outind",
                os.path.join(outdir, f"graph_traversal_path.{epoch}.txt"),
                "--outtxt",
                os.path.join(outdir, f"graph_traversal_zpath.{epoch}.txt"),
            ]
        )
        graph_traversal.main(args)

    def test_eval_volume(self, tmpdir_factory, particles, ctf, indices):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, ctf)
        args = eval_vol.add_args(argparse.ArgumentParser()).parse_args(
            [
                os.path.join(outdir, "weights.2.pkl"),
                "--config",
                os.path.join(outdir, "config.yaml"),
                "--zfile",
                os.path.join(outdir, "graph_traversal_zpath.2.txt"),
                "-o",
                os.path.join(outdir, "eval_vols"),
            ]
        )
        eval_vol.main(args)

        shutil.rmtree(outdir)


@pytest.mark.parametrize(
    "particles", ["toy.mrcs", "toy.star", "toy.txt"], indirect=True
)
@pytest.mark.parametrize("poses", ["toy-poses", "toy-angles"], indirect=True)
@pytest.mark.parametrize("batch_size", ["24", "36"], ids=("batch24", "batch36"))
@pytest.mark.parametrize("use_amp", [False, True], ids=("no-amp", "use-amp"))
def test_homogeneous_with_poses(tmpdir, particles, poses, batch_size, use_amp):
    args = [
        particles.path,
        "-o",
        str(tmpdir),
        "-n",
        "10",
        "--poses",
        poses.path,
        "--dim",
        "12",
        "-b",
        batch_size,
    ]
    if not use_amp:
        args += ["--no-amp"]

    train_nn.main(train_nn.add_args(argparse.ArgumentParser()).parse_args(args))
    assert "weights.9.pkl" in os.listdir(tmpdir)


@pytest.mark.parametrize(
    "abinit_dir", [{"zdim": zdim} for zdim in [0, 4, 8]], indirect=True
)
def test_abinit_checkpoint_analysis_and_backproject(abinit_dir):
    abinit_dir.train()
    abinit_dir.train(load_epoch=0)
    abinit_dir.backproject()
    abinit_dir.view_config()


@pytest.mark.parametrize(
    "train_dir",
    [{"train_cmd": "train_nn", "epochs": 5}, {"train_cmd": "train_vae", "epochs": 5}],
    indirect=True,
)
@pytest.mark.parametrize("load_epoch", [0, 2])
@pytest.mark.parametrize("train_epochs", [4, 5, 6])
def test_frompose_train_and_from_checkpoint(trained_dir, load_epoch, train_epochs):
    trained_dir.train_load_epoch(load_epoch, train_epochs)


@pytest.mark.parametrize("particles", ["tilts.star"], indirect=True)
@pytest.mark.parametrize("indices", [None, "just-4"], indirect=True)
@pytest.mark.parametrize("poses", ["tilt-poses"], indirect=True)
@pytest.mark.parametrize("ctf", ["CTF-Tilt"], indirect=True)
@pytest.mark.parametrize("datadir", ["default-datadir"], indirect=True)
class TestStarFixedHomo:
    """Run reconstructions using particles from a .star file as input."""

    def test_train_model(self, tmpdir, particles, indices, poses, ctf, datadir):
        args = [
            particles.path,
            "--datadir",
            datadir.path,
            "--poses",
            poses.path,
            "--ctf",
            ctf.path,
            "-o",
            str(tmpdir),
            "--dim",
            "32",
        ]
        if indices.path is not None:
            args += ["--ind", indices.path]

        args = train_nn.add_args(argparse.ArgumentParser()).parse_args(args)
        train_nn.main(args)


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

        args = train_vae.add_args(argparse.ArgumentParser()).parse_args(args)
        train_vae.main(args)

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
        train_vae.main(train_vae.add_args(argparse.ArgumentParser()).parse_args(args))

    def test_analyze(self, tmpdir_factory, particles, indices, poses, ctf, datadir):
        """Produce standard analyses for a particular epoch."""
        outdir = self.get_outdir(
            tmpdir_factory, particles, poses, ctf, indices, datadir
        )
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

    @pytest.mark.parametrize("nb_lbl", ["cryoDRGN_figures", "cryoDRGN_ET_viz"])
    def test_notebooks(
        self, tmpdir_factory, particles, indices, poses, ctf, datadir, nb_lbl
    ):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        outdir = self.get_outdir(
            tmpdir_factory, particles, poses, ctf, indices, datadir
        )
        orig_cwd = os.path.abspath(os.getcwd())
        os.chdir(os.path.join(outdir, "analyze.4"))
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
        os.chdir(os.path.join(outdir, "analyze.4"))
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

        args = train_vae.add_args(argparse.ArgumentParser()).parse_args(args)
        train_vae.main(args)
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

        args = abinit_homo.add_args(argparse.ArgumentParser()).parse_args(args)
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
        ]
        if indices.path is not None:
            args += ["--ind", indices.path]

        args = abinit_het.add_args(argparse.ArgumentParser()).parse_args(args)
        abinit_het.main(args)
