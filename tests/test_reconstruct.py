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

import cryodrgn.commands
import cryodrgn.commands_utils
from cryodrgn.source import ImageSource
from cryodrgn.trainers.reconstruction import ReconstructionModelConfigurations
import cryodrgn.utils
from cryodrgn.models.utils import get_model_configurations


class TrainCommand:

    replace_params = {
        "lr": "learning_rate",
        "wd": "weight_decay",
        "zdim": "z_dim",
        "dim": "hidden_dim",
    }
    bool_params = {"no-amp": {"amp", False}}

    def parse_args(self) -> ReconstructionModelConfigurations:
        cur_k = None
        i = 1
        while i < len(self.args):
            if self.args[i][0] == "-":
                cur_k = (
                    self.args[i][2:] if self.args[i][1] == "-" else self.args[i][1:]
                ).replace("-", "_")

                if cur_k in self.replace_params:
                    cur_k = self.replace_params[cur_k]
                if cur_k in self.bool_params:
                    self.cfgs.update(self.bool_params[cur_k])
            else:
                assert cur_k is not None
                self.cfgs[cur_k] = self.args[i]

            i += 1

        if "z_dim" not in self.cfgs:
            self.cfgs["z_dim"] = 0

        return get_model_configurations(self.cfgs)

    def __init__(
        self,
        train_cmd: str,
        args: list[str],
        outdir: str,
        train_type,
        poses=None,
    ) -> None:
        self.train_cmd = train_cmd
        self.args = args
        self.outdir = outdir
        self.train_type = train_type

        self.cfgs = {
            "model": "amort" if self.train_type == "drgnai" else "hps",
            "particles": self.args[0],
        }

        if self.train_cmd in {"train_nn", "train_vae"}:
            self.cfgs["pose_estimation"] = "fixed"
        elif self.train_cmd in {"abinit_homo", "abinit_het"}:
            self.cfgs["pose_estimation"] = "abinit"

        self.configs = self.parse_args()
        self.poses = poses

    def run(self) -> None:
        self.configs = self.parse_args()

        if self.train_type == "cdrgn":
            if "--hidden-dim" in self.args:
                i = self.args.index("--hidden-dim")
                self.args = (
                    self.args[:i]
                    + ["--enc-dim", self.args[i + 1], "--dec-dim", self.args[i + 1]]
                    + self.args[(i + 2) :]
                )

            cmd_args = self.args + ["-o", self.outdir]
        else:
            cfg_file = os.path.join(self.outdir, "configs.yaml")
            cryodrgn.utils.save_yaml(self.cfgs, cfg_file)
            cmd_args = [cfg_file, "-o", self.outdir, "--no-analysis"]

        use_cmd = self.train_cmd if self.train_type == "cdrgn" else "train"
        parser = argparse.ArgumentParser()
        eval(f"cryodrgn.commands.{use_cmd}").add_args(parser)
        eval(f"cryodrgn.commands.{use_cmd}").main(parser.parse_args(cmd_args))


@pytest.mark.parametrize("train_cmd", ["train_nn", "abinit_homo"])
@pytest.mark.parametrize(
    "train_type, particles, ctf, indices, poses, batch_size, use_amp",
    [
        ("drgnai", "hand", "CTF-Test.100", None, "hand-poses", "8", False),
        ("drgnai", "toy.mrcs", "CTF-Test", "random-100", "toy-poses", "8", False),
        ("drgnai", "toy.star", "CTF-Test", "first-100", "toy-poses", "16", True),
        ("drgnai", "toy.txt", "CTF-Test", "random-100", "toy-poses", "24", False),
        ("cdrgn-train", "hand", "CTF-Test.100", None, "hand-poses", "8", False),
        ("cdrgn-train", "toy.star", "CTF-Test", "first-100", "toy-poses", "16", True),
        ("cdrgn-train", "toy.txt", None, "random-100", "toy-poses", "24", False),
        ("cdrgn", "hand", "CTF-Test.100", None, "hand-poses", "36", False),
        ("cdrgn", "toy.mrcs", "CTF-Test", "first-100", "toy-poses", "24", True),
        ("cdrgn", "toy.mrcs", "CTF-Test", "first-100", "toy-angles", "24", True),
    ],
    indirect=["particles", "ctf", "indices", "poses"],
)
class TestHomogeneous:
    @pytest.fixture
    def traincmd(
        self,
        tmpdir_factory,
        train_type,
        train_cmd,
        particles,
        poses,
        ctf,
        indices,
        batch_size,
        use_amp,
    ):
        """Run an experiment to generate output; remove this output when finished."""
        if train_cmd == "train_nn":
            args = [
                "--pretrain",
                "50",
                "--num-epochs",
                "4",
                "--dim",
                "12",
                "--checkpoint",
                "1",
            ]
        elif train_cmd == "abinit_homo":
            args = [
                "--lr",
                ".01",
                "--dim",
                "8",
                "--layers",
                "2",
                "--pretrain",
                "50",
                "--num-epochs",
                "4",
                "--ps-freq",
                "2",
                "--checkpoint",
                "1",
            ]
        else:
            raise ValueError(
                f"Unrecognized homogeneous training command: `{train_cmd}`!"
            )

        dirname = os.path.join(
            "ReconstructHomo",
            train_type,
            train_cmd,
            particles.label,
            poses.label,
            ctf.label,
            indices.label,
            f"b{batch_size}_amp.{use_amp}",
        )
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)

        cmd_args = [particles.path, "-b", batch_size] + args
        if train_cmd == "train_nn":
            cmd_args += ["--poses", poses.path]
        if ctf.path is not None:
            cmd_args += ["--ctf", ctf.path]
        if indices.path is not None:
            cmd_args += ["--ind", indices.path]
        if not use_amp:
            cmd_args += ["--no-amp"]

        return TrainCommand(train_cmd, cmd_args, odir, train_type)

    def test_train_model(self, traincmd):
        """Train the initial homogeneous model."""

        traincmd.run()
        out_files = os.listdir(traincmd.outdir)
        assert "weights.4.pkl" in out_files, "Missing output model weights!"

    def test_train_from_checkpoint(self, traincmd):
        """Train the initial homogeneous model."""

        traincmd.args += ["--load", os.path.join(traincmd.outdir, "weights.4.pkl")]
        i = traincmd.args.index("--num-epochs")
        traincmd.args[i + 1] = str(int(traincmd.args[i + 1]) + 1)

        traincmd.run()
        out_files = os.listdir(traincmd.outdir)
        assert "weights.5.pkl" in out_files, "Missing output model weights!"


@pytest.mark.parametrize("train_cmd", ["train_vae", "abinit_het"])
class TestHeterogeneous:
    @pytest.fixture
    def traincmd(
        self,
        tmpdir_factory,
        train_type,
        train_cmd,
        particles,
        poses,
        ctf,
        indices,
    ) -> TrainCommand:
        """Run an experiment to generate output; remove this output when finished."""
        if train_cmd == "train_vae":
            args = [
                "--lr",
                ".0001",
                "--num-epochs",
                "3",
                "--seed",
                "0",
                "--zdim",
                "10",
                "--tdim",
                "8",
                "--hidden-dim",
                "8",
                "--pretrain",
                "50",
                "--pe-type",
                "gaussian",
                "--multigpu",
                "--checkpoint",
                "1",
            ]
        elif train_cmd == "abinit_het":
            args = [
                "--zdim",
                "4",
                "--lr",
                ".01",
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
                "50",
                "--num-epochs",
                "3",
                "--ps-freq",
                "2",
                "--checkpoint",
                "1",
            ]
        else:
            raise ValueError(
                f"Unrecognized heterogeneous training command: `{train_cmd}`!"
            )

        dirname = os.path.join(
            "ReconstructHetero",
            train_type,
            train_cmd,
            particles.label,
            poses.label,
            ctf.label,
            indices.label,
        )
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)

        cmd_args = [particles.path] + args
        if train_cmd == "train_vae":
            cmd_args += ["--poses", poses.path]
        if ctf.path is not None:
            cmd_args += ["--ctf", ctf.path]
        if indices.path is not None:
            cmd_args += ["--ind", indices.path]

        return TrainCommand(train_cmd, cmd_args, odir, train_type, poses=poses.path)

    @pytest.mark.parametrize(
        "train_type, particles, ctf, indices, poses",
        [
            ("drgnai", "hand", "CTF-Test.100", None, "hand-poses"),
            ("drgnai", "toy.mrcs", "CTF-Test", "random-100", "toy-poses"),
            ("drgnai", "toy.star", "CTF-Test", "first-100", "toy-poses"),
            ("drgnai", "toy.txt", "CTF-Test", "random-100", "toy-poses"),
            ("cdrgn-train", "hand", "CTF-Test.100", None, "hand-poses"),
            ("cdrgn-train", "toy.mrcs", "CTF-Test", "random-100", "toy-poses"),
            ("cdrgn-train", "toy.txt", None, "first-100", "toy-poses"),
            ("cdrgn", "hand", "CTF-Test.100", None, "hand-poses"),
            ("cdrgn", "toy.mrcs", "CTF-Test", "first-100", "toy-angles"),
        ],
        indirect=["particles", "ctf", "indices", "poses"],
    )
    def test_train_model(self, traincmd):
        """Train the initial heterogeneous model."""

        traincmd.run()
        out_files = os.listdir(traincmd.outdir)
        assert "weights.3.pkl" in out_files, "Missing output model weights!"
        assert "z.3.pkl" in out_files, "Missing output latent conformations!"

    @pytest.mark.parametrize(
        "train_type, particles, ctf, indices, poses",
        [
            ("drgnai", "hand", "CTF-Test.100", None, "hand-poses"),
            ("drgnai", "toy.mrcs", "CTF-Test", "random-100", "toy-poses"),
            ("drgnai", "toy.star", "CTF-Test", "first-100", "toy-poses"),
            ("drgnai", "toy.txt", "CTF-Test", "random-100", "toy-poses"),
            ("cdrgn-train", "hand", "CTF-Test.100", None, "hand-poses"),
            ("cdrgn-train", "toy.mrcs", "CTF-Test", "random-100", "toy-poses"),
            ("cdrgn-train", "toy.txt", None, "first-100", "toy-poses"),
            ("cdrgn", "hand", "CTF-Test.100", None, "hand-poses"),
            ("cdrgn", "toy.mrcs", "CTF-Test", "first-100", "toy-angles"),
        ],
        indirect=["particles", "ctf", "indices", "poses"],
    )
    def test_train_from_checkpoint(self, traincmd):
        """Load a cached model and run for another epoch, now without --multigpu."""

        traincmd.args += ["--load", os.path.join(traincmd.outdir, "weights.3.pkl")]
        if "--multigpu" in traincmd.args:
            i = traincmd.args.index("--multigpu")
            traincmd.args = traincmd.args[:i] + traincmd.args[(i + 1) :]
        i = traincmd.args.index("--num-epochs")
        traincmd.args[i + 1] = str(int(traincmd.args[i + 1]) + 1)

        # TODO: should number of epochs when loading be orig_epochs + new_epochs as in
        #       current DRGN-AI?
        traincmd.run()
        out_files = os.listdir(traincmd.outdir)
        assert "weights.4.pkl" in out_files, "Missing output model weights!"
        assert "z.4.pkl" in out_files, "Missing output latent conformations!"

    @pytest.mark.parametrize(
        "train_type, particles, ctf, indices, poses",
        [
            ("drgnai", "hand", "CTF-Test.100", None, "hand-poses"),
            ("drgnai", "toy.mrcs", "CTF-Test", "random-100", "toy-poses"),
            ("drgnai", "toy.star", "CTF-Test", "first-100", "toy-poses"),
            ("drgnai", "toy.txt", "CTF-Test", "random-100", "toy-poses"),
            ("cdrgn-train", "hand", "CTF-Test.100", None, "hand-poses"),
            ("cdrgn-train", "toy.mrcs", "CTF-Test", "random-100", "toy-poses"),
            ("cdrgn-train", "toy.txt", None, "first-100", "toy-poses"),
            ("cdrgn", "hand", "CTF-Test.100", None, "hand-poses"),
            ("cdrgn", "toy.mrcs", "CTF-Test", "first-100", "toy-angles"),
        ],
        indirect=["particles", "ctf", "indices", "poses"],
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
        cryodrgn.commands.analyze.add_args(parser)
        cryodrgn.commands.analyze.main(parser.parse_args(args))

        assert os.path.exists(os.path.join(traincmd.outdir, f"analyze.{epoch}"))

    @pytest.mark.parametrize(
        "train_type, particles, ctf, indices, poses",
        [
            ("drgnai", "hand", "CTF-Test.100", None, "hand-poses"),
            ("drgnai", "toy.mrcs", "CTF-Test", "random-100", "toy-poses"),
            ("drgnai", "toy.star", "CTF-Test", "first-100", "toy-poses"),
            ("drgnai", "toy.txt", "CTF-Test", "random-100", "toy-poses"),
            ("cdrgn-train", "hand", "CTF-Test.100", None, "hand-poses"),
            ("cdrgn-train", "toy.mrcs", "CTF-Test", "random-100", "toy-poses"),
            ("cdrgn", "hand", "CTF-Test.100", None, "hand-poses"),
        ],
        indirect=["particles", "ctf", "indices", "poses"],
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
        "train_type, particles, ctf, indices, poses, downsample_dim, flip_vol",
        [
            ("drgnai", "hand", "CTF-Test.100", None, "hand-poses", "8", False),
            ("drgnai", "toy.mrcs", "CTF-Test", "random-100", "toy-poses", "16", False),
            ("drgnai", "toy.star", "CTF-Test", "first-100", "toy-poses", "16", True),
            pytest.param(
                "drgnai",
                "toy.txt",
                "CTF-Test",
                "random-100",
                "toy-poses",
                "64",
                False,
                marks=pytest.mark.xfail(
                    raises=ValueError, reason="box size > resolution"
                ),
            ),
            ("cdrgn-train", "hand", "CTF-Test.100", None, "hand-poses", "16", True),
            pytest.param(
                "cdrgn-train",
                "toy.mrcs",
                "CTF-Test",
                "random-100",
                "toy-poses",
                None,
                False,
                marks=pytest.mark.xfail(
                    raises=ValueError, reason="box size > resolution"
                ),
            ),
            ("cdrgn", "hand", "CTF-Test.100", None, "hand-poses", "8", False),
        ],
        indirect=["particles", "ctf", "indices", "poses"],
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
        cryodrgn.commands.analyze_landscape.add_args(parser)
        cryodrgn.commands.analyze_landscape.main(parser.parse_args(args))

    @pytest.mark.parametrize(
        "train_type, particles, ctf, indices, poses, downsample_dim, flip_vol",
        [
            ("drgnai", "hand", "CTF-Test.100", None, "hand-poses", "8", False),
            ("drgnai", "toy.mrcs", "CTF-Test", "random-100", "toy-poses", "16", False),
            ("drgnai", "toy.star", "CTF-Test", "first-100", "toy-poses", "16", True),
            ("cdrgn-train", "hand", "CTF-Test.100", None, "hand-poses", "16", True),
            ("cdrgn", "hand", "CTF-Test.100", None, "hand-poses", "8", False),
        ],
        indirect=["particles", "ctf", "indices", "poses"],
    )
    def test_landscape_full(self, traincmd, downsample_dim, flip_vol):
        ldscp_dir = os.path.join(
            traincmd.outdir, f"landscape.4_{downsample_dim}.{flip_vol}"
        )
        args = [
            traincmd.outdir,
            "4",
            "-N",
            "10",
            "--landscape-dir",
            ldscp_dir,
        ]
        if downsample_dim is not None:
            args += ["--downsample", downsample_dim]
        if flip_vol:
            args += ["--flip"]

        parser = argparse.ArgumentParser()
        cryodrgn.commands.analyze_landscape_full.add_args(parser)
        cryodrgn.commands.analyze_landscape_full.main(parser.parse_args(args))

    @pytest.mark.parametrize(
        "train_type, particles, ctf, indices, poses, downsample_dim, flip_vol",
        [
            ("drgnai", "hand", "CTF-Test.100", None, "hand-poses", "8", False),
            ("drgnai", "toy.mrcs", "CTF-Test", "random-100", "toy-poses", "16", False),
            ("drgnai", "toy.star", "CTF-Test", "first-100", "toy-poses", "16", True),
            ("cdrgn-train", "hand", "CTF-Test.100", None, "hand-poses", "16", True),
            ("cdrgn", "hand", "CTF-Test.100", None, "hand-poses", "8", False),
        ],
        indirect=["particles", "ctf", "indices", "poses"],
    )
    def test_landscape_notebook(self, traincmd, downsample_dim, flip_vol):
        """Execute the demo Jupyter notebooks produced by landscape analysis."""
        orig_cwd = os.path.abspath(os.getcwd())
        outlbl = f"landscape.4_{downsample_dim}.{flip_vol}"
        os.chdir(os.path.join(traincmd.outdir, outlbl))
        nb_outfile = "analyze-landscape.ipynb"
        assert os.path.exists(nb_outfile)

        # edit the notebook with the epoch to analyze
        with open(nb_outfile, "r") as f:
            ntbook = nbformat.read(f, as_version=nbformat.NO_CONVERT)
        for cell in ntbook["cells"]:
            cell["source"] = cell["source"].replace("landscape.{EPOCH}", outlbl)

        try:
            ExecutePreprocessor(timeout=600, kernel_name="python3").preprocess(ntbook)
        except CellExecutionError as e:
            os.chdir(orig_cwd)
            raise e

        os.chdir(orig_cwd)

    @pytest.mark.parametrize(
        "train_type, particles, ctf, indices, poses, seed, steps, points",
        [
            ("drgnai", "hand", "CTF-Test.100", None, "hand-poses", 915, 5, None),
            ("drgnai", "toy.mrcs", "CTF-Test", "random-100", "toy-poses", 321, 2, None),
            ("drgnai", "toy.star", "CTF-Test", "first-100", "toy-poses", 701, 3, 1),
            ("drgnai", "toy.txt", "CTF-Test", "random-100", "toy-poses", 701, 3, 2),
            ("cdrgn-train", "hand", "CTF-Test.100", None, "hand-poses", 55, 3, None),
            (
                "cdrgn-train",
                "toy.mrcs",
                "CTF-Test",
                "random-100",
                "toy-poses",
                915,
                5,
                None,
            ),
            ("cdrgn-train", "toy.txt", None, "first-100", "toy-poses", 777, 4, None),
            ("cdrgn", "hand", "CTF-Test.100", None, "hand-poses", 404, 2, 1),
            ("cdrgn", "toy.mrcs", "CTF-Test", "first-100", "toy-angles", 555, 3, None),
        ],
        indirect=["particles", "ctf", "indices", "poses"],
    )
    def test_direct_traversal(self, traincmd, seed, steps, points):
        random.seed(seed)
        anchors = [str(anchor) for anchor in random.sample(range(100), steps)]

        parser = argparse.ArgumentParser()
        cryodrgn.commands.direct_traversal.add_args(parser)
        args = [os.path.join(traincmd.outdir, "z.3.pkl"), "--anchors"] + anchors
        if points is not None:
            args += ["-n", str(points)]

        cryodrgn.commands.direct_traversal.main(parser.parse_args(args))

    @pytest.mark.parametrize(
        "train_type, particles, ctf, indices, poses, epoch, seed, steps",
        [
            ("drgnai", "hand", "CTF-Test.100", None, "hand-poses", 4, 707, 4),
            ("drgnai", "toy.mrcs", "CTF-Test", "random-100", "toy-poses", 4, 707, 4),
            ("drgnai", "toy.star", "CTF-Test", "first-100", "toy-poses", 3, 607, 5),
            ("drgnai", "toy.txt", "CTF-Test", "random-100", "toy-poses", 4, 303, 2),
            ("cdrgn-train", "hand", "CTF-Test.100", None, "hand-poses", 4, 101, 3),
            (
                "cdrgn-train",
                "toy.mrcs",
                "CTF-Test",
                "random-100",
                "toy-poses",
                3,
                103,
                6,
            ),
            ("cdrgn-train", "toy.txt", None, "first-100", "toy-poses", 4, 707, 4),
            ("cdrgn", "hand", "CTF-Test.100", None, "hand-poses", 4, 915, 5),
            ("cdrgn", "toy.mrcs", "CTF-Test", "first-100", "toy-angles", 3, 321, 2),
        ],
        indirect=["particles", "ctf", "indices", "poses"],
    )
    def test_graph_traversal(self, traincmd, epoch, seed, steps):
        random.seed(seed)
        anchors = [str(anchor) for anchor in random.sample(range(100), steps)]

        parser = argparse.ArgumentParser()
        cryodrgn.commands.graph_traversal.add_args(parser)
        args = parser.parse_args(
            [
                os.path.join(traincmd.outdir, f"z.{epoch}.pkl"),
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
        cryodrgn.commands.graph_traversal.main(args)

    @pytest.mark.parametrize(
        "train_type, particles, ctf, indices, poses, epoch",
        [
            ("drgnai", "hand", "CTF-Test.100", None, "hand-poses", 4),
            ("drgnai", "toy.mrcs", "CTF-Test", "random-100", "toy-poses", 4),
            ("drgnai", "toy.star", "CTF-Test", "first-100", "toy-poses", 3),
            ("drgnai", "toy.txt", "CTF-Test", "random-100", "toy-poses", 4),
            ("cdrgn-train", "hand", "CTF-Test.100", None, "hand-poses", 4),
            ("cdrgn-train", "toy.mrcs", "CTF-Test", "random-100", "toy-poses", 3),
            ("cdrgn-train", "toy.txt", None, "first-100", "toy-poses", 4),
            ("cdrgn", "hand", "CTF-Test.100", None, "hand-poses", 4),
            ("cdrgn", "toy.mrcs", "CTF-Test", "first-100", "toy-angles", 3),
        ],
        indirect=["particles", "ctf", "indices", "poses"],
    )
    def test_eval_volume(self, traincmd, epoch):
        parser = argparse.ArgumentParser()
        cryodrgn.commands.eval_vol.add_args(parser)
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
        cryodrgn.commands.eval_vol.main(args)

    @pytest.mark.parametrize(
        "train_type, particles, ctf, indices, poses",
        [
            ("drgnai", "hand", "CTF-Test.100", None, "hand-poses"),
            ("drgnai", "toy.mrcs", "CTF-Test", "random-100", "toy-poses"),
            ("drgnai", "toy.star", "CTF-Test", "first-100", "toy-poses"),
            ("drgnai", "toy.txt", "CTF-Test", "random-100", "toy-poses"),
            ("cdrgn-train", "hand", "CTF-Test.100", None, "hand-poses"),
            ("cdrgn-train", "toy.mrcs", "CTF-Test", "random-100", "toy-poses"),
            ("cdrgn-train", "toy.txt", None, "first-100", "toy-poses"),
            ("cdrgn", "hand", "CTF-Test.100", None, "hand-poses"),
            ("cdrgn", "toy.mrcs", "CTF-Test", "first-100", "toy-angles"),
        ],
        indirect=["particles", "ctf", "indices", "poses"],
    )
    @pytest.mark.parametrize("epoch", [2, 3])
    def test_eval_images(self, traincmd, epoch):
        parser = argparse.ArgumentParser()
        cryodrgn.commands.eval_images.add_args(parser)
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
                traincmd.poses,
                "--log-interval",
                "1",
                "--verbose",
            ]
        )
        cryodrgn.commands.eval_images.main(args)

    @pytest.mark.parametrize(
        "train_type, particles, ctf, indices, poses",
        [
            ("drgnai", "hand", "CTF-Test.100", None, "hand-poses"),
            ("drgnai", "toy.star", "CTF-Test", "first-100", "toy-poses"),
            ("cdrgn-train", "toy.txt", None, "first-100", "toy-poses"),
            ("cdrgn", "toy.mrcs", "CTF-Test", "first-100", "toy-angles"),
        ],
        indirect=["particles", "ctf", "indices", "poses"],
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
        cryodrgn.commands_utils.plot_classes.add_args(parser)
        args = [traincmd.outdir, str(epoch), "--labels", lbl_file]
        if palette is not None:
            args += ["--palette", palette]
        if plot_outdir is not None:
            args += ["--outdir", os.path.join(traincmd.outdir, plot_outdir)]

        cryodrgn.commands_utils.plot_classes.main(parser.parse_args(args))
        if plot_outdir is not None:
            use_outdir = os.path.join(traincmd.outdir, plot_outdir)
        elif epoch == -1:
            use_outdir = os.path.join(traincmd.outdir, "analyze.4")
        else:
            use_outdir = os.path.join(traincmd.outdir, f"analyze.{epoch}")

        assert os.path.exists(os.path.join(use_outdir, "umap_kde_classes.png"))

    @pytest.mark.parametrize(
        "train_type, particles, ctf, indices, poses",
        [
            ("drgnai", "hand", "CTF-Test.100", None, "hand-poses"),
            ("drgnai", "toy.mrcs", "CTF-Test", "random-100", "toy-poses"),
            ("drgnai", "toy.star", "CTF-Test", "first-100", "toy-poses"),
            ("drgnai", "toy.txt", "CTF-Test", "random-100", "toy-poses"),
            ("cdrgn-train", "hand", "CTF-Test.100", None, "hand-poses"),
            ("cdrgn-train", "toy.mrcs", "CTF-Test", "random-100", "toy-poses"),
            ("cdrgn-train", "toy.txt", None, "first-100", "toy-poses"),
            ("cdrgn", "hand", "CTF-Test.100", None, "hand-poses"),
            ("cdrgn", "toy.mrcs", "CTF-Test", "first-100", "toy-angles"),
        ],
        indirect=["particles", "ctf", "indices", "poses"],
    )
    def test_clean_all(self, traincmd):
        parser = argparse.ArgumentParser()
        cryodrgn.commands_utils.clean.add_args(parser)
        cryodrgn.commands_utils.clean.main(
            parser.parse_args([os.path.relpath(traincmd.outdir)])
        )

        shutil.rmtree(traincmd.outdir)
