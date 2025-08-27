"""Running ab-initio volume reconstruction followed by downstream analyses."""

import pytest
import argparse
import os.path
import shutil
import pickle
import numpy as np
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
from cryodrgn.commands import analyze, eval_vol, filter, graph_traversal, abinit_het


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
    ids=[
        "hand,no.ind,no.ctf",
        "hand,no.ind,with.ctf",
        "toy.txt,ind.rand.100,with.ctf",
        "toy.star,ind.f100,with.ctf",
        "toy.star-13,no.ind,no.ctf",
    ],
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

        parser = argparse.ArgumentParser()
        abinit_het.add_args(parser)
        abinit_het.main(parser.parse_args(args))
        assert not os.path.exists(os.path.join(outdir, "analyze.2"))

    @pytest.mark.parametrize(
        "epoch, vol_start_index",
        [(2, 1), (1, 0)],
        ids=["epoch.2,volstart.1", "epoch.1,volstart.0"],
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

    @pytest.mark.parametrize("epoch", [2, 3], ids=["epoch.2", "epoch.3"])
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
    "abinit_dir", [{"zdim": zdim} for zdim in [0, 4, 8]], indirect=True
)
def test_abinit_checkpoint_analysis_and_backproject(abinit_dir):
    abinit_dir.train()
    abinit_dir.train(load_epoch=1)
    abinit_dir.backproject()
