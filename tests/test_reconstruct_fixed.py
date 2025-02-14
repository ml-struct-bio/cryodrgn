"""Volume reconstruction training with fixed poses, as well as downstream analyses."""

import pytest
import argparse
import os.path
import shutil
import random
import numpy as np

import nbformat
from nbclient.exceptions import CellExecutionError
from nbconvert.preprocessors import ExecutePreprocessor

from cryodrgn.commands import (
    analyze,
    analyze_landscape,
    analyze_landscape_full,
    direct_traversal,
    eval_images,
    eval_vol,
    graph_traversal,
)
from cryodrgn.commands_utils import clean, plot_classes
from cryodrgn.source import ImageSource
import cryodrgn.utils
from tests.conftest import TrainCommand


@pytest.mark.parametrize(
    "particles", ["toy.mrcs", "toy.star", "toy.txt"], indirect=True
)
@pytest.mark.parametrize(
    "indices, batch_size, use_amp",
    [(None, "36", False), ("random-100", "24", True)],
    indirect=["indices"],
)
class TestFixedHomo:
    @pytest.fixture
    def traincmd(
        self,
        tmpdir_factory,
        train_type,
        particles,
        poses,
        ctf,
        indices,
        batch_size,
        use_amp,
    ):
        """Run an experiment to generate output; remove this output when finished."""

        dirname = os.path.join(
            "FixedHomo",
            train_type,
            particles.label,
            poses.label,
            ctf.label,
            indices.label,
        )
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)

        cmd_args = [
            particles.path,
            "--num-epochs",
            "4",
            "--poses",
            poses.path,
            "--dim",
            "12",
            "-b",
            batch_size,
            "--checkpoint",
            "1",
        ]
        if not use_amp:
            cmd_args += ["--no-amp"]
        if ctf.path is not None:
            cmd_args += ["--ctf", ctf.path]
        if indices.path is not None:
            cmd_args += ["--ind", indices.path]

        return TrainCommand("train_nn", cmd_args, odir, train_type)

    @pytest.mark.parametrize(
        "train_type, ctf, poses",
        [
            ("drgnai", "CTF-Test", "toy-poses"),
            ("cdrgn", None, "toy-poses"),
            ("cdrgn", "CTF-Test", "toy-poses"),
            ("cdrgn-train", None, "toy-poses"),
            ("cdrgn-train", "CTF-Test", "toy-angles"),
        ],
        indirect=["ctf", "poses"],
    )
    def test_train_model(self, traincmd):
        """Train the initial homogeneous model."""

        traincmd.run()
        out_files = os.listdir(traincmd.outdir)
        assert "weights.4.pkl" in out_files, "Missing output model weights!"

    @pytest.mark.parametrize(
        "train_type, ctf, poses",
        [
            ("drgnai", "CTF-Test", "toy-poses"),
            ("cdrgn", None, "toy-poses"),
            ("cdrgn", "CTF-Test", "toy-poses"),
            ("cdrgn-train", None, "toy-poses"),
            ("cdrgn-train", "CTF-Test", "toy-angles"),
        ],
        indirect=["ctf", "poses"],
    )
    def test_train_from_checkpoint(self, traincmd):
        """Train the initial homogeneous model."""

        traincmd.args += ["--load", os.path.join(traincmd.outdir, "weights.4.pkl")]
        i = traincmd.args.index("--num-epochs")
        traincmd.args[i + 1] = str(int(traincmd.args[i + 1]) + 1)

        traincmd.run()
        out_files = os.listdir(traincmd.outdir)
        assert "weights.5.pkl" in out_files, "Missing output model weights!"


@pytest.mark.parametrize(
    "particles", ["toy.mrcs", "toy.txt", "toy.star"], indirect=True
)
@pytest.mark.parametrize("poses", ["toy-poses"], indirect=True)
@pytest.mark.parametrize("indices", [None, "random-100"], indirect=True)
class TestFixedHetero:
    @pytest.fixture
    def traincmd(self, tmpdir_factory, train_type, particles, poses, ctf, indices):
        """Run an experiment to generate output; remove this output when finished."""

        dirname = os.path.join(
            "FixedHetero",
            train_type,
            particles.label,
            poses.label,
            ctf.label,
            indices.label,
        )
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)

        cmd_args = [
            particles.path,
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
            "--hidden-dim",
            "8",
            "--pe-type",
            "gaussian",
            "--multigpu",
            "--checkpoint",
            "1",
        ]
        if ctf.path is not None:
            cmd_args += ["--ctf", ctf.path]
        if indices.path is not None:
            cmd_args += ["--ind", indices.path]

        return TrainCommand("train_vae", cmd_args, odir, train_type)

    @pytest.mark.parametrize(
        "train_type, ctf",
        [
            ("drgnai", "CTF-Test"),
            ("cdrgn", None),
            ("cdrgn", "CTF-Test"),
            ("cdrgn-train", None),
            ("cdrgn-train", "CTF-Test"),
        ],
        indirect=["ctf"],
    )
    def test_train_model(self, traincmd):
        """Train the initial heterogeneous model."""

        traincmd.run()
        out_files = os.listdir(traincmd.outdir)
        assert "weights.3.pkl" in out_files, "Missing output model weights!"
        assert "conf.3.pkl" in out_files, "Missing output latent conformations!"

    @pytest.mark.parametrize(
        "train_type, ctf",
        [
            ("drgnai", "CTF-Test"),
            ("cdrgn", None),
            ("cdrgn", "CTF-Test"),
            ("cdrgn-train", None),
            ("cdrgn-train", "CTF-Test"),
        ],
        indirect=["ctf"],
    )
    def test_train_from_checkpoint(self, traincmd):
        """Load a cached model and run for another epoch, now without --multigpu."""

        traincmd.args += ["--load", os.path.join(traincmd.outdir, "weights.3.pkl")]
        i = traincmd.args.index("--multigpu")
        traincmd.args = traincmd.args[:i] + traincmd.args[(i + 1) :]
        i = traincmd.args.index("--num-epochs")
        traincmd.args[i + 1] = str(int(traincmd.args[i + 1]) + 1)

        # TODO: should number of epochs when loading be orig_epochs + new_epochs as in
        #       current DRGN-AI?
        traincmd.run()
        out_files = os.listdir(traincmd.outdir)
        assert "weights.4.pkl" in out_files, "Missing output model weights!"
        assert "conf.4.pkl" in out_files, "Missing output latent conformations!"

    @pytest.mark.parametrize(
        "train_type, ctf",
        [
            ("drgnai", "CTF-Test"),
            ("cdrgn", None),
            ("cdrgn", "CTF-Test"),
            ("cdrgn-train", None),
            ("cdrgn-train", "CTF-Test"),
        ],
        indirect=["ctf"],
    )
    @pytest.mark.parametrize("epoch", [1, 3, None])
    def test_analyze(self, traincmd, epoch):
        """Produce standard analyses for a particular epoch."""

        args = [traincmd.outdir]
        if epoch is not None:
            args += ["--epoch", str(epoch)]
        else:
            epoch = 4

        parser = argparse.ArgumentParser()
        analyze.add_args(parser)
        analyze.main(parser.parse_args(args))

        assert os.path.exists(os.path.join(traincmd.outdir, f"analyze.{epoch}"))

    @pytest.mark.parametrize(
        "train_type, ctf",
        [
            ("drgnai", "CTF-Test"),
            ("cdrgn", "CTF-Test"),
            ("cdrgn-train", "CTF-Test"),
        ],
        indirect=["ctf"],
    )
    @pytest.mark.parametrize("nb_lbl", ["analysis", "cryoDRGN_figures", "cryoDRGN_viz"])
    def test_notebooks(self, traincmd, nb_lbl):
        """Execute the demonstration Jupyter notebooks produced by analysis."""

        orig_cwd = os.path.abspath(os.getcwd())
        os.chdir(os.path.join(traincmd.outdir, "analyze.3"))
        assert os.path.exists(f"{nb_lbl}.ipynb"), "Upstream tests have failed!"

        with open(f"{nb_lbl}.ipynb") as ff:
            nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

        try:
            ExecutePreprocessor(timeout=600, kernel_name="python3").preprocess(nb_in)
        except CellExecutionError as e:
            os.chdir(orig_cwd)
            raise e

        os.chdir(orig_cwd)

    @pytest.mark.parametrize(
        "train_type, ctf, downsample_dim, flip_vol",
        [
            ("drgnai", "CTF-Test", "16", False),
            ("drgnai", "CTF-Test", "16", True),
            ("cdrgn", None, "16", False),
            ("cdrgn", "CTF-Test", "16", True),
            pytest.param(
                "cdrgn",
                "CTF-Test",
                "64",
                False,
                marks=pytest.mark.xfail(
                    raises=ValueError, reason="box size > resolution"
                ),
            ),
            pytest.param(
                "cdrgn",
                "CTF-Test",
                None,
                False,
                marks=pytest.mark.xfail(
                    raises=ValueError, reason="box size > resolution"
                ),
            ),
            ("cdrgn-train", "CTF-Test", "16", False),
        ],
        indirect=["ctf"],
    )
    def test_landscape(self, traincmd, downsample_dim, flip_vol):
        ldscp_dir = os.path.join(
            traincmd.outdir, f"landscape.4_{downsample_dim}.{flip_vol}"
        )
        args = [
            traincmd.outdir,
            "4",  # Epoch number to analyze - 0-indexed
            "--sketch-size",
            "10",  # Number of volumes to generate for analysis
            "--pc-dim",
            "5",
            "--vol-start-index",
            "1",
            "-o",
            ldscp_dir,
        ]
        if downsample_dim:
            args += ["--downsample", downsample_dim]
        if flip_vol:
            args += ["--flip"]

        parser = argparse.ArgumentParser()
        analyze_landscape.add_args(parser)
        analyze_landscape.main(parser.parse_args(args))

    @pytest.mark.parametrize(
        "train_type, ctf, downsample_dim, flip_vol",
        [
            ("drgnai", "CTF-Test", "16", False),
            ("drgnai", "CTF-Test", "16", True),
            ("cdrgn", None, "16", False),
            ("cdrgn", "CTF-Test", "16", True),
            ("cdrgn-train", "CTF-Test", "16", False),
        ],
        indirect=["ctf"],
    )
    def test_landscape_full(self, traincmd, downsample_dim, flip_vol):
        ldscp_dir = os.path.join(
            traincmd.outdir, f"landscape.4_{downsample_dim}.{flip_vol}"
        )
        ldscpf_dir = os.path.join(
            traincmd.outdir, f"landscape-full.4_{downsample_dim}.{flip_vol}"
        )
        args = [
            traincmd.outdir,
            "4",
            "-N",
            "10",
            "--landscape-dir",
            ldscp_dir,
            "-o",
            ldscpf_dir,
        ]
        if downsample_dim is not None:
            args += ["--downsample", downsample_dim]
        if flip_vol:
            args += ["--flip"]

        parser = argparse.ArgumentParser()
        analyze_landscape_full.add_args(parser)
        analyze_landscape_full.main(parser.parse_args(args))

    @pytest.mark.parametrize(
        "train_type, ctf",
        [("drgnai", "CTF-Test"), ("cdrgn", "CTF-Test")],
        indirect=["ctf"],
    )
    def test_landscape_notebook(self, traincmd):
        """Execute the demo Jupyter notebooks produced by landscape analysis."""
        orig_cwd = os.path.abspath(os.getcwd())
        outlbl = f"landscape.4_{16}.{False}"
        os.chdir(os.path.join(traincmd.outdir, outlbl))
        nb_outfile = "analyze-landscape.ipynb"
        assert os.path.exists(nb_outfile)

        # edit the notebook with the epoch to analyze
        with open(nb_outfile, "r") as f:
            ntbook = nbformat.read(f, as_version=nbformat.NO_CONVERT)
        for cell in ntbook["cells"]:
            cell["source"] = cell["source"].replace("landscape.{EPOCH}", outlbl)
        for cell in ntbook["cells"]:
            cell["source"] = cell["source"].replace(
                "{landscape_dir}/landscape_full",
                f"{{landscape_dir}}/{outlbl}",
            )

        try:
            ExecutePreprocessor(timeout=600, kernel_name="python3").preprocess(ntbook)
        except CellExecutionError as e:
            os.chdir(orig_cwd)
            raise e

        os.chdir(orig_cwd)

    @pytest.mark.parametrize(
        "train_type, ctf, seed, steps, points",
        [
            ("cdrgn", None, 915, 5, None),
            ("drgnai", "CTF-Test", 321, 2, None),
            ("cdrgn-train", "CTF-Test", 701, 3, 1),
            ("drgnai", "CTF-Test", 701, 3, 2),
            ("cdrgn-train", "CTF-Test", 55, 3, None),
        ],
        indirect=["ctf"],
    )
    def test_direct_traversal(self, traincmd, seed, steps, points):
        random.seed(seed)
        anchors = [str(anchor) for anchor in random.sample(range(100), steps)]

        parser = argparse.ArgumentParser()
        direct_traversal.add_args(parser)
        args = [os.path.join(traincmd.outdir, "conf.3.pkl"), "--anchors"] + anchors
        if points is not None:
            args += ["-n", str(points)]

        direct_traversal.main(parser.parse_args(args))

    @pytest.mark.parametrize(
        "train_type, ctf, epoch, seed, steps",
        [
            ("cdrgn", None, 3, 915, 5),
            ("drgnai", "CTF-Test", 2, 321, 2),
            ("cdrgn-train", "CTF-Test", 3, 701, 5),
            ("drgnai", "CTF-Test", 3, 102, 10),
        ],
        indirect=["ctf"],
    )
    def test_graph_traversal(self, traincmd, epoch, seed, steps):
        random.seed(seed)
        anchors = [str(anchor) for anchor in random.sample(range(100), steps)]

        parser = argparse.ArgumentParser()
        graph_traversal.add_args(parser)
        args = parser.parse_args(
            [
                os.path.join(traincmd.outdir, f"conf.{epoch}.pkl"),
                "--anchors",
            ]
            + anchors
            + [
                "--outind",
                os.path.join(traincmd.outdir, f"graph_traversal_path.{epoch}.txt"),
                "--outtxt",
                os.path.join(traincmd.outdir, f"graph_traversal_zpath.{epoch}.txt"),
            ]
        )
        graph_traversal.main(args)

    @pytest.mark.parametrize(
        "train_type, ctf",
        [("drgnai", "CTF-Test"), ("cdrgn-train", "CTF-Test"), ("cdrgn", "CTF-Test")],
        indirect=True,
    )
    @pytest.mark.parametrize("epoch", [2, 3])
    def test_eval_volume(self, traincmd, epoch):
        parser = argparse.ArgumentParser()
        eval_vol.add_args(parser)
        args = parser.parse_args(
            [
                os.path.join(traincmd.outdir, f"weights.{epoch}.pkl"),
                "--config",
                os.path.join(traincmd.outdir, "train-configs.yaml"),
                "--zfile",
                os.path.join(traincmd.outdir, f"graph_traversal_zpath.{epoch}.txt"),
                "-o",
                os.path.join(traincmd.outdir, "eval_vols"),
            ]
        )
        eval_vol.main(args)

    @pytest.mark.parametrize(
        "train_type, ctf",
        [
            ("drgnai", "CTF-Test"),
            ("cdrgn", None),
            ("cdrgn", "CTF-Test"),
            ("cdrgn-train", None),
            ("cdrgn-train", "CTF-Test"),
        ],
        indirect=["ctf"],
    )
    @pytest.mark.parametrize("epoch", [2, 3])
    def test_eval_images(self, traincmd, epoch):
        parser = argparse.ArgumentParser()
        eval_images.add_args(parser)
        args = parser.parse_args(
            [
                traincmd.configs.particles,
                os.path.join(traincmd.outdir, f"weights.{epoch}.pkl"),
                "--config",
                os.path.join(traincmd.outdir, "train-configs.yaml"),
                "-o",
                os.path.join(traincmd.outdir, f"out_eval_images_losses.{epoch}.pkl"),
                "--out-z",
                os.path.join(traincmd.outdir, f"out_eval_images_z.{epoch}.pkl"),
                "--poses",
                traincmd.configs.poses,
                "--log-interval",
                "1",
                "--verbose",
            ]
        )
        eval_images.main(args)

    @pytest.mark.parametrize(
        "train_type, ctf",
        [
            ("drgnai", "CTF-Test"),
            ("cdrgn", None),
            ("cdrgn", "CTF-Test"),
            ("cdrgn-train", None),
            ("cdrgn-train", "CTF-Test"),
        ],
        indirect=["ctf"],
    )
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
    def test_plot_classes(self, traincmd, epoch, palette, plot_outdir):
        if traincmd.configs.ind is not None:
            ind = cryodrgn.utils.load_pkl(traincmd.configs.ind)
        else:
            ind = None

        particles_n = ImageSource.from_file(
            traincmd.configs.particles, indices=ind, lazy=True
        ).n
        lbl_file = os.path.join(traincmd.outdir, "plot-classes.pkl")
        if not os.path.exists(lbl_file):
            labels = np.array([0 if i % 3 == 0 else 1 for i in range(particles_n)])
            cryodrgn.utils.save_pkl(labels, lbl_file)

        parser = argparse.ArgumentParser()
        plot_classes.add_args(parser)
        args = [traincmd.outdir, str(epoch), "--labels", lbl_file]
        if palette is not None:
            args += ["--palette", palette]
        if plot_outdir is not None:
            args += ["--outdir", os.path.join(traincmd.outdir, plot_outdir)]

        plot_classes.main(parser.parse_args(args))
        if plot_outdir is not None:
            use_outdir = os.path.join(traincmd.outdir, plot_outdir)
        elif epoch == -1:
            use_outdir = os.path.join(traincmd.outdir, "analyze.4")
        else:
            use_outdir = os.path.join(traincmd.outdir, f"analyze.{epoch}")

        assert os.path.exists(os.path.join(use_outdir, "umap_kde_classes.png"))

    @pytest.mark.parametrize(
        "train_type, ctf",
        [
            ("drgnai", "CTF-Test"),
            ("cdrgn", None),
            ("cdrgn", "CTF-Test"),
            ("cdrgn-train", None),
            ("cdrgn-train", "CTF-Test"),
        ],
        indirect=["ctf"],
    )
    def test_clean_all(self, traincmd):
        parser = argparse.ArgumentParser()
        clean.add_args(parser)
        clean.main(parser.parse_args([os.path.relpath(traincmd.outdir)]))

        shutil.rmtree(traincmd.outdir)
