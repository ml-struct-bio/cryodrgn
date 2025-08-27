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
    direct_traversal,
    eval_images,
    eval_vol,
    filter,
    graph_traversal,
    train_nn,
    train_vae,
    train_dec,
)
from cryodrgn.commands_utils import clean, plot_classes
from cryodrgn.source import ImageSource
from cryodrgn import utils


@pytest.mark.parametrize("train_cmd", ["train_vae", "train_dec"])
@pytest.mark.parametrize(
    "particles", ["toy.mrcs", "toy.txt", "toy.star"], indirect=True
)
@pytest.mark.parametrize("poses", ["toy-poses"], indirect=True)
@pytest.mark.parametrize(
    "indices", [None, "random-100"], indirect=True, ids=["no.indices", "with.indices"]
)
class TestFixedHetero:
    def get_outdir(self, tmpdir_factory, train_cmd, particles, poses, ctf, indices):
        dirname = os.path.join(
            "FixHet", train_cmd, particles.label, poses.label, ctf.label, indices.label
        )
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)

        return odir

    @pytest.mark.parametrize(
        "ctf", [None, "CTF-Test"], indirect=True, ids=["no.CTF", "with.CTF"]
    )
    def test_train_model(
        self, tmpdir_factory, train_cmd, particles, poses, ctf, indices
    ):
        """Train the initial heterogeneous model."""

        outdir = self.get_outdir(
            tmpdir_factory, train_cmd, particles, indices, poses, ctf
        )
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
            "--pe-type",
            "gaussian",
            "--multigpu",
            "--no-analysis",
        ]
        if train_cmd == "train_vae":
            args += ["--tdim", "8", "--enc-dim", "8", "--dec-dim", "8"]
        elif train_cmd == "train_dec":
            args += ["--dim", "8", "--layers", "2"]
        if ctf.path is not None:
            args += ["--ctf", ctf.path]
        if indices.path is not None:
            args += ["--ind", indices.path]

        parser = argparse.ArgumentParser()
        train_module = (
            train_vae
            if train_cmd == "train_vae"
            else train_dec
            if train_cmd == "train_dec"
            else None
        )
        train_module.add_args(parser)
        train_module.main(parser.parse_args(args))
        assert os.path.exists(os.path.join(outdir, "weights.3.pkl"))
        assert not os.path.exists(os.path.join(outdir, "weights.4.pkl"))
        assert not os.path.exists(os.path.join(outdir, "analyze.3"))

    @pytest.mark.parametrize(
        "ctf, load",
        [(None, "latest"), ("CTF-Test", 2)],
        indirect=["ctf"],
        ids=["no.CTF,load.latest", "with.CTF,load.epoch.2"],
    )
    def test_train_from_checkpoint(
        self,
        tmpdir_factory,
        train_cmd,
        particles,
        poses,
        ctf,
        indices,
        load,
    ):
        """Load a cached model and run for another epoch, now without --multigpu."""

        outdir = self.get_outdir(
            tmpdir_factory, train_cmd, particles, indices, poses, ctf
        )
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
            "--pe-type",
            "gaussian",
            "--load",
            load,
        ]
        if train_cmd == "train_vae":
            args += ["--tdim", "8", "--enc-dim", "8", "--dec-dim", "8"]
        elif train_cmd == "train_dec":
            args += ["--dim", "8", "--layers", "2"]
        if ctf.path is not None:
            args += ["--ctf", ctf.path]
        if indices.path is not None:
            args += ["--ind", indices.path]

        parser = argparse.ArgumentParser()
        train_module = (
            train_vae
            if train_cmd == "train_vae"
            else train_dec
            if train_cmd == "train_dec"
            else None
        )
        train_module.add_args(parser)
        train_module.main(parser.parse_args(args))
        assert os.path.exists(os.path.join(outdir, "weights.4.pkl"))
        assert not os.path.exists(os.path.join(outdir, "weights.5.pkl"))

    @pytest.mark.parametrize(
        "ctf, epoch", [("CTF-Test", 3), (None, 4), ("CTF-Test", 4)], indirect=["ctf"]
    )
    def test_analyze(
        self, tmpdir_factory, train_cmd, particles, poses, ctf, indices, epoch
    ):
        """Produce standard analyses for a particular epoch."""

        outdir = self.get_outdir(
            tmpdir_factory, train_cmd, particles, indices, poses, ctf
        )
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
    def test_notebooks(
        self, tmpdir_factory, train_cmd, particles, poses, ctf, indices, nb_lbl
    ):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        outdir = self.get_outdir(
            tmpdir_factory, train_cmd, particles, indices, poses, ctf
        )
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

    @pytest.mark.parametrize("plotind", [False, True], ids=["dontsave.ind", "save.ind"])
    @pytest.mark.parametrize(
        "ctf, epoch",
        [("CTF-Test", 3), ("CTF-Test", None), (None, None)],
        indirect=["ctf"],
        ids=["with.CTF,epoch.3", "with.CTF,last.epoch", "no.CTF,last.epoch"],
    )
    def test_interactive_filtering(
        self,
        tmpdir_factory,
        train_cmd,
        particles,
        poses,
        ctf,
        indices,
        epoch,
        plotind,
    ):
        """Launch interface for filtering particles using model covariates."""
        outdir = self.get_outdir(
            tmpdir_factory, train_cmd, particles, indices, poses, ctf
        )
        parser = argparse.ArgumentParser()
        filter.add_args(parser)
        args = [outdir, "--force"]
        if epoch is not None:
            args += ["--epoch", str(epoch)]
            sel_dir = os.path.join(outdir, f"analyze.{epoch}")
        else:
            sel_dir = os.path.join(outdir, "analyze.4")
        args += ["--sel-dir", sel_dir]

        if plotind:
            ind_fl = os.path.join(outdir, "tmp_ind_test.pkl")
            with open(ind_fl, "wb") as f:
                pickle.dump(np.array([1, 2]), f)
            args += ["--plot-inds", ind_fl]

        filter.main(parser.parse_args(args))
        if plotind:
            assert os.path.exists(os.path.join(sel_dir, "indices.pkl"))
            with open(os.path.join(sel_dir, "indices.pkl"), "rb") as f:
                inds = pickle.load(f)
            assert isinstance(inds, np.ndarray)
            assert len(inds) == 2
            assert os.path.exists(os.path.join(sel_dir, "indices_inverse.pkl"))
            with open(os.path.join(sel_dir, "indices_inverse.pkl"), "rb") as f:
                inv_inds = pickle.load(f)
            assert isinstance(inv_inds, np.ndarray)
        else:
            assert not os.path.exists(os.path.join(sel_dir, "indices.pkl"))
            assert not os.path.exists(os.path.join(sel_dir, "indices_inverse.pkl"))

    @pytest.mark.parametrize(
        "ctf, downsample_dim, flip_vol, sketch_size, num_clusters",
        [
            (None, "16", False, 10, 3),
            ("CTF-Test", "16", True, 5, 5),
            pytest.param(
                "CTF-Test",
                None,
                False,
                5,
                5,
                marks=pytest.mark.xfail(
                    raises=ValueError, reason="box size > resolution"
                ),
            ),
        ],
        indirect=["ctf"],
        ids=[
            "no.CTF,downsample.16,flip.False,sketch.10,clusters.3",
            "with.CTF,downsample.16,flip.True,sketch.5,clusters.5",
            "with.CTF,downsample.None,flip.False,sketch.5,clusters.5",
        ],
    )
    def test_landscape(
        self,
        tmpdir_factory,
        train_cmd,
        particles,
        poses,
        ctf,
        indices,
        downsample_dim,
        flip_vol,
        sketch_size,
        num_clusters,
    ):
        outdir = self.get_outdir(
            tmpdir_factory, train_cmd, particles, indices, poses, ctf
        )
        args = [
            outdir,
            "4",  # Epoch number to analyze - 1-indexed
            "--sketch-size",
            str(sketch_size),  # Number of volumes to generate for analysis
            "--pc-dim",
            "5",
            "-M",
            str(num_clusters),
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
        ],
        indirect=["ctf"],
        ids=[
            "no.CTF,downsample.16,flip.False",
            "with.CTF,downsample.16,flip.True",
            "with.CTF,downsample.64,flip.False",
        ],
    )
    def test_landscape_full(
        self,
        tmpdir_factory,
        train_cmd,
        particles,
        poses,
        ctf,
        indices,
        downsample_dim,
        flip_vol,
    ):
        outdir = self.get_outdir(
            tmpdir_factory, train_cmd, particles, indices, poses, ctf
        )
        args = [outdir, "4", "-N", "20"]
        if downsample_dim is not None:
            args += ["--downsample", downsample_dim]
        if flip_vol:
            args += ["--flip"]

        parser = argparse.ArgumentParser()
        analyze_landscape_full.add_args(parser)
        analyze_landscape_full.main(parser.parse_args(args))

    @pytest.mark.parametrize("ctf", ["CTF-Test"], indirect=True, ids=["with.CTF"])
    def test_landscape_notebook(
        self, tmpdir_factory, train_cmd, particles, poses, ctf, indices
    ):
        """Execute the demo Jupyter notebooks produced by landscape analysis."""

        outdir = self.get_outdir(
            tmpdir_factory, train_cmd, particles, indices, poses, ctf
        )
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
            (None, 915, 2, None),
            ("CTF-Test", 321, 3, None),
            ("CTF-Test", 701, 2, 3),
            ("CTF-Test", 544, 3, 10),
        ],
        indirect=["ctf"],
        ids=[
            "no.CTF,two.steps,default.points",
            "with.CTF,three.steps,default.points",
            "with.CTF,two.steps,three.points",
            "with.CTF,three.steps,ten.points",
        ],
    )
    def test_direct_traversal(
        self,
        tmpdir_factory,
        train_cmd,
        particles,
        poses,
        ctf,
        indices,
        seed,
        steps,
        points,
    ):
        outdir = self.get_outdir(
            tmpdir_factory, train_cmd, particles, indices, poses, ctf
        )
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
        self,
        tmpdir_factory,
        train_cmd,
        particles,
        poses,
        ctf,
        indices,
        epoch,
        seed,
        steps,
    ):
        outdir = self.get_outdir(
            tmpdir_factory, train_cmd, particles, indices, poses, ctf
        )
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

    @pytest.mark.parametrize("ctf", ["CTF-Test"], indirect=True, ids=["with.CTF"])
    @pytest.mark.parametrize("epoch", [3, 4])
    def test_eval_volume(
        self, tmpdir_factory, train_cmd, particles, poses, ctf, indices, epoch
    ):
        outdir = self.get_outdir(
            tmpdir_factory, train_cmd, particles, indices, poses, ctf
        )
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

    @pytest.mark.parametrize("ctf", ["CTF-Test"], indirect=True, ids=["with.CTF"])
    @pytest.mark.parametrize("epoch", [3, 4])
    def test_eval_images(
        self, tmpdir_factory, train_cmd, particles, poses, ctf, indices, epoch
    ):
        outdir = self.get_outdir(
            tmpdir_factory, train_cmd, particles, indices, poses, ctf
        )
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
        if train_cmd == "train_dec":
            with pytest.raises(NotImplementedError):
                eval_images.main(args)
        else:
            eval_images.main(args)

    @pytest.mark.parametrize("ctf", ["CTF-Test"], indirect=True, ids=["with.CTF"])
    @pytest.mark.parametrize(
        "epoch, palette, plot_outdir, plot_types",
        [
            (-1, "rocket", None, None),
            pytest.param(
                4,
                "Rocket",
                None,
                None,
                marks=pytest.mark.xfail(
                    raises=ValueError, reason="palette not available in seaborn!"
                ),
            ),
            (3, None, None, ["kde", "scatter"]),
            (4, None, "plots", ["kde"]),
        ],
        ids=[
            "plot.defaults",
            "plot.defaults,bad.palette",
            "all.plots,epoch.3",
            "kdeplot.only,chosen.outdir",
        ],
    )
    def test_plot_classes(
        self,
        tmpdir_factory,
        train_cmd,
        particles,
        poses,
        ctf,
        indices,
        epoch,
        palette,
        plot_outdir,
        plot_types,
    ):
        outdir = self.get_outdir(
            tmpdir_factory, train_cmd, particles, indices, poses, ctf
        )
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
        if plot_types is not None:
            args += ["--plot-types"] + plot_types

        plot_classes.main(parser.parse_args(args))
        if plot_outdir is not None:
            use_outdir = os.path.join(outdir, plot_outdir)
        elif epoch == -1:
            use_outdir = os.path.join(outdir, "analyze")
        else:
            use_outdir = os.path.join(outdir, f"analyze.{epoch}")

        if plot_types is None:
            plot_types = ["scatter"]

        scatter_fl = os.path.join(use_outdir, "umap_scatter_classes.png")
        kde_fl = os.path.join(use_outdir, "umap_kde_classes.png")
        assert os.path.exists(scatter_fl) == ("scatter" in plot_types)
        assert os.path.exists(kde_fl) == ("kde" in plot_types)

    @pytest.mark.parametrize("ctf", [None, "CTF-Test"], indirect=True)
    def test_clean_all(self, tmpdir_factory, train_cmd, particles, poses, ctf, indices):
        outdir = self.get_outdir(
            tmpdir_factory, train_cmd, particles, indices, poses, ctf
        )
        parser = argparse.ArgumentParser()
        clean.add_args(parser)

        clean.main(parser.parse_args([os.path.relpath(outdir)]))
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

    parser = argparse.ArgumentParser()
    train_nn.add_args(parser)
    train_nn.main(parser.parse_args(args))
    assert "weights.10.pkl" in os.listdir(tmpdir)


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

        parser = argparse.ArgumentParser()
        train_nn.add_args(parser)
        train_nn.main(parser.parse_args(args))
