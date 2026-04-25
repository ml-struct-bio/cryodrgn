"""Running ab-initio volume reconstruction followed by downstream analyses."""

import pytest
import argparse
import os.path
import pickle
import shutil
import numpy as np
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
from cryodrgn.commands import (
    analyze,
    abinit,
    analyze_landscape,
    analyze_landscape_full,
    eval_vol,
    filter,
    graph_traversal,
)


@pytest.mark.parametrize(
    "particles, indices, ctf",
    [
        ("hand", None, "CTF-Test.100"),
        ("toy.txt", "random-100", "CTF-Test"),
    ],
    indirect=True,
)
class TestAbinitHomo:

    model_args = [
        "--zdim",
        "0",
        "--lr",
        ".001",
        "--dim",
        "16",
        "--layers",
        "2",
        "--pe-dim",
        "4",
        "--t-extent",
        "4.0",
        "--t-ngrid",
        "2",
        "--nkeptposes",
        "4",
    ]

    def get_outdir(self, tmpdir_factory, particles, ctf, indices):
        dirname = os.path.join("AbinitHomo", particles.label, ctf.label, indices.label)
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)

        return odir

    def test_train_model(self, tmpdir_factory, particles, ctf, indices):
        """Train the initial homogeneous model."""

        outdir = self.get_outdir(tmpdir_factory, particles, indices, ctf)
        args = [
            particles.path,
            "-o",
            outdir,
            *self.model_args,
            "--num-epochs",
            "3",
            "--epochs-pose-search",
            "1",
            "--n-imgs-pretrain",
            "10",
            "--no-analysis",
            "--checkpoint",
            "1",
        ]
        if ctf.path is not None:
            args += ["--ctf", ctf.path]
        if indices.path is not None:
            args += ["--ind", indices.path]

        parser = argparse.ArgumentParser()
        abinit.add_args(parser)
        abinit.main(parser.parse_args(args))

        assert os.path.exists(os.path.join(outdir, "weights.2.pkl"))
        assert os.path.exists(os.path.join(outdir, "weights.3.pkl"))
        assert not os.path.exists(os.path.join(outdir, "analyze.2"))
        assert not os.path.exists(os.path.join(outdir, "analyze.3"))

    def test_load_checkpoint(self, tmpdir_factory, particles, ctf, indices):
        """Load a checkpoint and continue training."""

        outdir = self.get_outdir(tmpdir_factory, particles, indices, ctf)
        new_outdir = os.path.join(outdir, "checkpoint")
        parser = argparse.ArgumentParser()
        abinit.add_args(parser)
        args = [
            particles.path,
            "-o",
            new_outdir,
            "--load",
            os.path.join(outdir, "weights.3.pkl"),
            *self.model_args,
            "--num-epochs",
            "5",
            "--checkpoint",
            "3",
        ]
        if ctf.path is not None:
            args += ["--ctf", ctf.path]
        if indices.path is not None:
            args += ["--ind", indices.path]

        abinit.main(parser.parse_args(args))
        assert not os.path.exists(os.path.join(new_outdir, "weights.4.pkl"))
        assert os.path.exists(os.path.join(new_outdir, "weights.5.pkl"))
        assert not os.path.exists(os.path.join(new_outdir, "analyze.4"))
        assert os.path.exists(os.path.join(new_outdir, "analyze.5"))


@pytest.mark.parametrize(
    "particles, indices, ctf",
    [
        ("hand", None, "CTF-Test.100"),
        ("toy.txt", "random-100", "CTF-Test"),
        ("toy.star", "first-100", "CTF-Test"),
    ],
    indirect=True,
    ids=["hand,no.ind", "toy.txt,ind.rand.100", "toy.star,ind.f100"],
)
class TestAbinitHetero:

    model_args = [
        "--zdim",
        "4",
        "--lr",
        ".001",
        "--dim",
        "16",
        "--layers",
        "2",
        "--pe-dim",
        "4",
        "--t-extent",
        "4.0",
        "--t-ngrid",
        "2",
        "--nkeptposes",
        "4",
    ]

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
            *self.model_args,
            "--num-epochs",
            "3",
            "--epochs-pose-search",
            "1",
            "--n-imgs-pretrain",
            "10",
            "--no-analysis",
            "--checkpoint",
            "1",
        ]
        if ctf.path is not None:
            args += ["--ctf", ctf.path]
        if indices.path is not None:
            args += ["--ind", indices.path]

        parser = argparse.ArgumentParser()
        abinit.add_args(parser)
        abinit.main(parser.parse_args(args))

        assert os.path.exists(os.path.join(outdir, "weights.2.pkl"))
        assert os.path.exists(os.path.join(outdir, "weights.3.pkl"))
        assert not os.path.exists(os.path.join(outdir, "analyze.2"))
        assert not os.path.exists(os.path.join(outdir, "analyze.3"))

    @pytest.mark.parametrize(
        "epoch, vol_start_index",
        [(3, 1), (2, 1)],
        ids=["epoch.3,volstart.1", "epoch.2,volstart.1"],
    )
    def test_analyze(
        self, tmpdir_factory, particles, ctf, indices, epoch, vol_start_index
    ):
        """Produce standard analyses for a particular epoch."""

        outdir = self.get_outdir(tmpdir_factory, particles, indices, ctf)
        parser = argparse.ArgumentParser()
        analyze.add_args(parser)
        analyze.main(
            parser.parse_args(
                [
                    outdir,
                    str(epoch),  # Epoch number to analyze - 1-indexed
                    "--pc",
                    "3",  # Number of principal component traversals to generate
                    "--ksample",
                    "10",  # Number of kmeans samples to generate
                    "--vol-start-index",
                    str(vol_start_index),
                ]
            )
        )

        kmeans_dir = os.path.join(outdir, f"analyze.{epoch}", "kmeans10")
        for i in range(vol_start_index, 10 + vol_start_index):
            assert os.path.exists(os.path.join(kmeans_dir, f"vol_{i:03d}.mrc"))
        assert not os.path.exists(
            os.path.join(kmeans_dir, f"vol_{(10 + vol_start_index):03d}.mrc")
        )

    def test_load_checkpoint(self, tmpdir_factory, particles, ctf, indices):
        """Load a checkpoint and continue training."""

        outdir = self.get_outdir(tmpdir_factory, particles, indices, ctf)
        new_outdir = os.path.join(outdir, "checkpoint")
        parser = argparse.ArgumentParser()
        abinit.add_args(parser)
        args = [
            particles.path,
            "-o",
            new_outdir,
            "--load",
            os.path.join(outdir, "weights.3.pkl"),
            *self.model_args,
            "--num-epochs",
            "5",
            "--checkpoint",
            "3",
        ]
        if ctf.path is not None:
            args += ["--ctf", ctf.path]
        if indices.path is not None:
            args += ["--ind", indices.path]

        abinit.main(parser.parse_args(args))
        assert not os.path.exists(os.path.join(new_outdir, "weights.4.pkl"))
        assert os.path.exists(os.path.join(new_outdir, "weights.5.pkl"))
        assert not os.path.exists(os.path.join(new_outdir, "analyze.4"))
        assert os.path.exists(os.path.join(new_outdir, "analyze.5"))

    def test_load_poses(self, tmpdir_factory, particles, ctf, indices):
        """Load poses and continue training."""

        outdir = self.get_outdir(tmpdir_factory, particles, indices, ctf)
        new_outdir = os.path.join(outdir, "checkpoint-poses")
        parser = argparse.ArgumentParser()
        abinit.add_args(parser)
        args = [
            particles.path,
            "-o",
            new_outdir,
            "--load-poses",
            os.path.join(outdir, "pose.1.pkl"),
            *self.model_args,
            "--num-epochs",
            "4",
            "--epochs-pose-search",
            "2",
            "--checkpoint",
            "2",
            "--no-analysis",
        ]
        if ctf.path is not None:
            args += ["--ctf", ctf.path]
        if indices.path is not None:
            args += ["--ind", indices.path]

        abinit.main(parser.parse_args(args))
        assert not os.path.exists(os.path.join(new_outdir, "weights.3.pkl"))
        assert os.path.exists(os.path.join(new_outdir, "weights.4.pkl"))
        assert not os.path.exists(os.path.join(new_outdir, "analyze.4"))

    @pytest.mark.parametrize("nb_lbl", ["cryoDRGN_figures", "cryoDRGN_filtering"])
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

    @pytest.mark.parametrize("plotind", [False, True], ids=["dontsave.ind", "save.ind"])
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
        ids=["epoch.2", "epoch.None"],
    )
    def test_interactive_filtering(
        self, tmpdir_factory, particles, ctf, indices, epoch, plotind
    ):
        """Launch interface for filtering particles using model covariates."""

        outdir = self.get_outdir(tmpdir_factory, particles, indices, ctf)
        parser = argparse.ArgumentParser()
        filter.add_args(parser)
        args = [outdir, "--force"]
        if epoch is not None:
            args += ["--epoch", str(epoch)]
            sel_dir = os.path.join(outdir, f"analyze.{epoch}")
        else:
            sel_dir = os.path.join(outdir, "analyze.3")
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

    def test_graph_traversal(self, tmpdir_factory, particles, ctf, indices):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, ctf)
        parser = argparse.ArgumentParser()
        graph_traversal.add_args(parser)
        args = parser.parse_args(
            [
                os.path.join(outdir, "z.3.pkl"),
                "--anchors",
                "3",
                "5",
                "8",
                "--max-neighbors",
                "50",
                "--avg-neighbors",
                "50",
                "--outind",
                os.path.join(outdir, "graph_traversal_path.3.txt"),
                "--outtxt",
                os.path.join(outdir, "graph_traversal_zpath.3.txt"),
            ]
        )
        graph_traversal.main(args)

    def test_analyze_landscape(self, tmpdir_factory, particles, ctf, indices):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, ctf)
        parser = argparse.ArgumentParser()
        analyze_landscape.add_args(parser)
        args = parser.parse_args(
            [
                outdir,
                "3",
                "--sketch-size",
                "10",
                "-M",
                "3",
                "--pc-dim",
                "5",
                "--downsample",
                "64",
            ]
        )
        analyze_landscape.main(args)

        assert os.path.exists(os.path.join(outdir, "landscape.3"))
        assert os.path.exists(os.path.join(outdir, "landscape.3", "umap.pkl"))
        assert os.path.exists(os.path.join(outdir, "landscape.3", "vol_pca_obj.pkl"))
        assert os.path.exists(os.path.join(outdir, "landscape.3", "kmeans10"))
        for i in range(1, 11):
            assert os.path.exists(
                os.path.join(outdir, "landscape.3", "kmeans10", f"vol_{i:03d}.mrc")
            )

    def test_analyze_landscape_full(self, tmpdir_factory, particles, ctf, indices):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, ctf)
        parser = argparse.ArgumentParser()
        analyze_landscape_full.add_args(parser)
        args = parser.parse_args([outdir, "3", "-N", "10", "--downsample", "64"])
        analyze_landscape_full.main(args)

        landfull_dir = os.path.join(outdir, "landscape.3", "landscape_full")
        assert os.path.exists(os.path.join(landfull_dir, "vol_pca_sampled.pkl"))
        assert os.path.exists(os.path.join(landfull_dir, "z.sampled.pkl"))

    def test_eval_volume(self, tmpdir_factory, particles, ctf, indices):
        """Parity with legacy abinit: eval_vol along graph-traversal z-path (new abinit config)."""
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
        assert os.path.isdir(os.path.join(outdir, "eval_vols"))
        shutil.rmtree(outdir)
