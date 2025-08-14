"""Running an experiment of training followed by downstream analyses."""

import pytest
import argparse
import os.path
import shutil
import random
import nbformat
from nbclient.exceptions import CellExecutionError
from nbconvert.preprocessors import ExecutePreprocessor
import numpy as np

from cryodrgn.commands import (
    analyze,
    analyze_landscape,
    analyze_landscape_full,
    direct_traversal,
    eval_images,
    eval_vol,
    filter,
    graph_traversal,
    train_nn,
    train_vae,
    abinit_het,
)
from cryodrgn.commands_utils import clean, plot_classes
from cryodrgn.source import ImageSource
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
            "--no-analysis",
        ]
        if ctf.path is not None:
            args += ["--ctf", ctf.path]
        if indices.path is not None:
            args += ["--ind", indices.path]

        train_vae.main(train_vae.add_args(argparse.ArgumentParser()).parse_args(args))
        assert os.path.exists(os.path.join(outdir, "weights.3.pkl"))
        assert not os.path.exists(os.path.join(outdir, "weights.4.pkl"))
        assert not os.path.exists(os.path.join(outdir, "analyze.3"))

    @pytest.mark.parametrize(
        "ctf, load", [(None, "latest"), ("CTF-Test", 2)], indirect=["ctf"]
    )
    def test_train_from_checkpoint(
        self,
        tmpdir_factory,
        particles,
        poses,
        ctf,
        indices,
        load,
    ):
        """Load a cached model and run for another epoch, now without --multigpu."""
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        if isinstance(load, int):
            load = os.path.join(outdir, f"weights.{load}.pkl")

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
            load,
        ]
        if ctf.path is not None:
            args += ["--ctf", ctf.path]
        if indices.path is not None:
            args += ["--ind", indices.path]

        train_vae.main(train_vae.add_args(argparse.ArgumentParser()).parse_args(args))
        assert os.path.exists(os.path.join(outdir, "weights.4.pkl"))
        assert not os.path.exists(os.path.join(outdir, "weights.5.pkl"))

    @pytest.mark.parametrize(
        "ctf, epoch", [("CTF-Test", 3), (None, 4), ("CTF-Test", 4)], indirect=["ctf"]
    )
    def test_analyze(self, tmpdir_factory, particles, poses, ctf, indices, epoch):
        """Produce standard analyses for a particular epoch."""
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        parser = argparse.ArgumentParser()
        analyze.add_args(parser)
        analyze.main(parser.parse_args([outdir, str(epoch)]))

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
        os.chdir(os.path.join(outdir, "analyze.4"))
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
        "ctf, epoch",
        [
            ("CTF-Test", 4),
            ("CTF-Test", None),
            pytest.param(
                None,
                None,
                marks=pytest.mark.xfail(raises=NotImplementedError),
            ),
        ],
        indirect=["ctf"],
    )
    def test_interactive_filtering(
        self, tmpdir_factory, particles, poses, ctf, indices, epoch
    ):
        """Launch interface for filtering particles using model covariates."""
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        parser = argparse.ArgumentParser()
        filter.add_args(parser)
        epoch_args = ["--epoch", str(epoch)] if epoch is not None else list()
        filter.main(parser.parse_args([outdir] + epoch_args))

    @pytest.mark.parametrize(
        "ctf, downsample_dim, flip_vol",
        [
            (None, "16", False),
            ("CTF-Test", "16", True),
            pytest.param(
                "CTF-Test",
                "64",
                False,
                marks=pytest.mark.xfail(
                    raises=ValueError, reason="box size > resolution"
                ),
            ),
            pytest.param(
                "CTF-Test",
                None,
                False,
                marks=pytest.mark.xfail(
                    raises=ValueError, reason="box size > resolution"
                ),
            ),
        ],
        indirect=["ctf"],
    )
    def test_landscape(
        self, tmpdir_factory, particles, poses, ctf, indices, downsample_dim, flip_vol
    ):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        args = [
            outdir,
            "4",  # Epoch number to analyze - 1-indexed
            "--sketch-size",
            "10",  # Number of volumes to generate for analysis
            "--pc-dim",
            "5",
            "--vol-start-index",
            "1",
        ]
        if downsample_dim:
            args += ["--downsample", downsample_dim]
        if flip_vol:
            args += ["--flip"]

        parser = argparse.ArgumentParser()
        analyze_landscape.add_args(parser)
        analyze_landscape.main(parser.parse_args(args))

    @pytest.mark.parametrize(
        "ctf, downsample_dim, flip_vol",
        [
            (None, "16", False),
            ("CTF-Test", "16", True),
            pytest.param(
                "CTF-Test",
                "64",
                False,
                marks=pytest.mark.xfail(
                    raises=AssertionError, reason="box size > resolution"
                ),
            ),
            pytest.param(
                "CTF-Test",
                None,
                False,
                marks=pytest.mark.xfail(
                    raises=AssertionError, reason="box size > resolution"
                ),
            ),
        ],
        indirect=["ctf"],
    )
    def test_landscape_full(
        self, tmpdir_factory, particles, poses, ctf, indices, downsample_dim, flip_vol
    ):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        args = [outdir, "4", "-N", "10"]
        if downsample_dim is not None:
            args += ["--downsample", downsample_dim]
        if flip_vol:
            args += ["--flip"]

        parser = argparse.ArgumentParser()
        analyze_landscape_full.add_args(parser)
        analyze_landscape_full.main(parser.parse_args(args))

    @pytest.mark.parametrize("ctf", ["CTF-Test"], indirect=True)
    def test_landscape_notebook(self, tmpdir_factory, particles, poses, ctf, indices):
        """Execute the demo Jupyter notebooks produced by landscape analysis."""
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        orig_cwd = os.path.abspath(os.getcwd())
        os.chdir(os.path.join(outdir, "landscape.4"))
        notebook_fl = "cryoDRGN_analyze_landscape.ipynb"
        assert os.path.exists(notebook_fl)

        with open(notebook_fl) as ff:
            nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

        try:
            ExecutePreprocessor(timeout=600, kernel_name="python3").preprocess(nb_in)
        except CellExecutionError as e:
            os.chdir(orig_cwd)
            raise e

        os.chdir(orig_cwd)

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
        args = [os.path.join(outdir, "z.4.pkl"), "--anchors"] + anchors
        if points is not None:
            args += ["-n", str(points)]

        direct_traversal.main(parser.parse_args(args))

    @pytest.mark.parametrize(
        "ctf, epoch, seed, steps",
        [(None, 4, 915, 5), ("CTF-Test", 3, 321, 2), ("CTF-Test", 4, 655, 3)],
        indirect=["ctf"],
    )
    def test_graph_traversal(
        self, tmpdir_factory, particles, poses, ctf, indices, epoch, seed, steps
    ):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        random.seed(seed)

        if steps == 2:
            anchors = ["95", "35"]
        else:
            anchors = ["95"]
            anchors += [
                str(anchor)
                for anchor in random.sample(
                    list(set(range(100)) - {"95", "35"}), steps - 2
                )
            ]
            anchors += ["35"]

        parser = argparse.ArgumentParser()
        graph_traversal.add_args(parser)
        args = parser.parse_args(
            [
                os.path.join(outdir, f"z.{epoch}.pkl"),
                "--anchors",
            ]
            + anchors
            + ["--max-neighbors", "20", "--avg-neighbors", "10"]
            + [
                "--outind",
                os.path.join(outdir, f"graph_traversal_path.{epoch}.txt"),
                "--outtxt",
                os.path.join(outdir, f"graph_traversal_zpath.{epoch}.txt"),
            ]
        )
        graph_traversal.main(args)

    @pytest.mark.parametrize("ctf", ["CTF-Test"], indirect=True)
    @pytest.mark.parametrize("epoch", [3, 4])
    def test_eval_volume(self, tmpdir_factory, particles, poses, ctf, indices, epoch):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        parser = argparse.ArgumentParser()
        eval_vol.add_args(parser)
        args = parser.parse_args(
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
    @pytest.mark.parametrize("epoch", [3, 4])
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
                4,
                "Rocket",
                None,
                marks=pytest.mark.xfail(
                    raises=ValueError, reason="palette not available in seaborn!"
                ),
            ),
            (3, None, None),
            (3, None, "plots"),
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
            "--no-analysis",
        ]
        if ctf.path is not None:
            args += ["--ctf", ctf.path]
        if indices.path is not None:
            args += ["--ind", indices.path]

        abinit_het.main(abinit_het.add_args(argparse.ArgumentParser()).parse_args(args))
        assert not os.path.exists(os.path.join(outdir, "analyze.2"))

    def test_analyze(self, tmpdir_factory, particles, ctf, indices):
        """Produce standard analyses for a particular epoch."""
        outdir = self.get_outdir(tmpdir_factory, particles, indices, ctf)

        parser = argparse.ArgumentParser()
        analyze.add_args(parser)
        analyze.main(
            parser.parse_args(
                [
                    outdir,
                    "2",  # Epoch number to analyze - 1-indexed
                    "--pc",
                    "3",  # Number of principal component traversals to generate
                    "--ksample",
                    "10",  # Number of kmeans samples to generate
                    "--vol-start-index",
                    "1",
                ]
            )
        )

        assert os.path.exists(os.path.join(outdir, "analyze.2"))

    @pytest.mark.parametrize("nb_lbl", ["cryoDRGN_figures"])
    def test_notebooks(self, tmpdir_factory, particles, ctf, indices, nb_lbl):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        outdir = self.get_outdir(tmpdir_factory, particles, indices, ctf)
        orig_cwd = os.path.abspath(os.getcwd())
        os.chdir(os.path.join(outdir, "analyze.2"))
        assert os.path.exists(f"{nb_lbl}.ipynb"), "Upstream tests have failed!"

        with open(f"{nb_lbl}.ipynb") as ff:
            nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

        ExecutePreprocessor(timeout=600, kernel_name="python3").preprocess(nb_in)
        os.chdir(orig_cwd)

    @pytest.mark.parametrize(
        "epoch",
        [
            2,
            pytest.param(
                None,
                marks=pytest.mark.xfail(
                    raises=ValueError, reason="missing analysis epoch"
                ),
            ),
        ],
    )
    def test_interactive_filtering(
        self, tmpdir_factory, particles, ctf, indices, epoch
    ):
        """Launch interface for filtering particles using model covariates."""
        outdir = self.get_outdir(tmpdir_factory, particles, indices, ctf)
        parser = argparse.ArgumentParser()
        filter.add_args(parser)
        epoch_args = ["--epoch", str(epoch)] if epoch is not None else list()
        filter.main(parser.parse_args([outdir] + epoch_args))

    @pytest.mark.parametrize("epoch", [2, 3])
    def test_graph_traversal(self, tmpdir_factory, particles, ctf, indices, epoch):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, ctf)
        parser = argparse.ArgumentParser()
        graph_traversal.add_args(parser)
        args = parser.parse_args(
            [
                os.path.join(outdir, f"z.{epoch}.pkl"),
                "--anchors",
                "5",
                "2",
                "10",
                "--max-neighbors",
                "100",
                "--avg-neighbors",
                "100",
                "--outind",
                os.path.join(outdir, f"graph_traversal_path.{epoch}.txt"),
                "--outtxt",
                os.path.join(outdir, f"graph_traversal_zpath.{epoch}.txt"),
            ]
        )
        graph_traversal.main(args)

    def test_eval_volume(self, tmpdir_factory, particles, ctf, indices):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, ctf)
        parser = argparse.ArgumentParser()
        eval_vol.add_args(parser)
        args = parser.parse_args(
            [
                os.path.join(outdir, "weights.3.pkl"),
                "--config",
                os.path.join(outdir, "config.yaml"),
                "--zfile",
                os.path.join(outdir, "graph_traversal_zpath.3.txt"),
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
    assert "weights.10.pkl" in os.listdir(tmpdir)


@pytest.mark.parametrize(
    "abinit_dir", [{"zdim": zdim} for zdim in [0, 4, 8]], indirect=True
)
def test_abinit_checkpoint_analysis_and_backproject(abinit_dir):
    abinit_dir.train()
    abinit_dir.train(load_epoch=1)
    abinit_dir.backproject()
    abinit_dir.view_config()


@pytest.mark.parametrize(
    "train_dir",
    [{"train_cmd": "train_nn", "epochs": 5}, {"train_cmd": "train_vae", "epochs": 5}],
    indirect=True,
)
@pytest.mark.parametrize("load_epoch", [1, 3])
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
