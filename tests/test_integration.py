"""Integration tests of ab initio reconstruction and downstream analyses.

Note that the training done here has unrealistically low parameter values to allow the
tests to run in any environment in a reasonable amount of time with or without GPUs.

"""
import shutil

import pytest
import os
import argparse
import pickle
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
import numpy as np
import torch
from pathlib import Path

from cryodrgn.commands import (
    analyze,
    backproject_voxel,
    downsample,
    parse_ctf_star,
    parse_pose_star,
    train_vae,
)
from cryodrgn.commands_utils import write_star
from cryodrgn.source import ImageSource
import cryodrgn.utils


@pytest.mark.parametrize("particles", ["toy.mrcs"], indirect=True)
@pytest.mark.parametrize("poses", ["toy-poses"], indirect=True)
@pytest.mark.parametrize("ctf", ["CTF-Test"], indirect=True)
@pytest.mark.parametrize("indices", [None, "first-100", "random-100"], indirect=True)
class TestIterativeFiltering:
    def get_outdir(self, tmpdir_factory, particles, indices, poses, ctf):
        dirname = Path(
            "IterativeFiltering",
            "_".join([particles.label, indices.label, poses.label, ctf.label]),
        )
        odir = tmpdir_factory.getbasetemp() / dirname
        os.makedirs(odir, exist_ok=True)

        return odir

    def test_train_model(self, tmpdir_factory, particles, poses, ctf, indices):
        """Train the initial heterogeneous model without any manual filters."""
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        args = [
            particles.path,
            "-o",
            str(outdir),
            "--ctf",
            ctf.path,
            "--num-epochs",
            "10",
            "--poses",
            poses.path,
            "--zdim",
            "4",
            "--tdim",
            "64",
            "--enc-dim",
            "64",
            "--dec-dim",
            "64",
            "--pe-type",
            "gaussian",
            "--no-amp",
        ]
        if indices.path is not None:
            args += ["--ind", indices.path]

        parser = argparse.ArgumentParser()
        train_vae.add_args(parser)
        train_vae.main(parser.parse_args(args))

    def test_analyze(self, tmpdir_factory, particles, poses, ctf, indices):
        """Produce standard analyses for the final epoch."""
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        assert os.path.exists(
            os.path.join(outdir, "weights.10.pkl")
        ), "Upstream tests have failed!"

        parser = argparse.ArgumentParser()
        analyze.add_args(parser)
        analyze.main(parser.parse_args([str(outdir), "--pc", "2", "--ksample", "10"]))

        assert (outdir / "analyze.10").exists(), "Model analyses have failed!"

    @pytest.mark.parametrize("nb_lbl", ["cryoDRGN_figures", "analysis", "cryoDRGN_viz"])
    def test_notebooks(self, tmpdir_factory, particles, poses, ctf, indices, nb_lbl):
        """Execute the demonstration Jupyter notebooks produced by analysis."""
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        orig_cwd = os.path.abspath(os.getcwd())
        os.chdir(outdir / "analyze.10")
        assert os.path.exists(f"{nb_lbl}.ipynb"), "Upstream tests have failed!"

        with open(f"{nb_lbl}.ipynb") as ff:
            nb_in = nbformat.read(ff, nbformat.NO_CONVERT)

        ExecutePreprocessor(timeout=60, kernel_name="python3").preprocess(nb_in)
        os.chdir(orig_cwd)

    def test_refiltering(self, tmpdir_factory, particles, poses, ctf, indices):
        outdir = self.get_outdir(tmpdir_factory, particles, indices, poses, ctf)
        ind_keep_fl = [fl for fl in os.listdir(outdir) if fl[:9] == "ind_keep."]

        if (
            not ind_keep_fl
            or int(ind_keep_fl[0].split(".")[1].split("_particles")[0]) < 20
        ):
            ind_keep_fl = [fl for fl in os.listdir(outdir) if fl[:8] == "ind_bad."][0]
        else:
            ind_keep_fl = ind_keep_fl[0]

        ind_keep_fl = outdir / ind_keep_fl
        parser = argparse.ArgumentParser()
        train_vae.add_args(parser)
        train_vae.main(
            parser.parse_args(
                [
                    particles.path,
                    "-o",
                    str(outdir),
                    "--ctf",
                    ctf.path,
                    "--ind",
                    str(ind_keep_fl),
                    "--num-epochs",
                    "5",
                    "--poses",
                    poses.path,
                    "--zdim",
                    "4",
                    "--tdim",
                    "64",
                    "--enc-dim",
                    "64",
                    "--dec-dim",
                    "64",
                    "--pe-type",
                    "gaussian",
                    "--no-amp",
                ]
            )
        )

        shutil.rmtree(outdir)


@pytest.mark.parametrize("particles", ["toy.star"], indirect=True)
@pytest.mark.parametrize("datadir", ["default-datadir"], indirect=True)
class TestParseWriteStar:
    def get_outdir(self, tmpdir_factory, particles, datadir) -> Path:
        dirname = Path("ParseWriteStar", "_".join([particles.label, datadir.label]))
        odir = Path(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)

        return odir

    def test_parse_ctf_star(self, tmpdir_factory, particles, datadir):
        outdir = self.get_outdir(tmpdir_factory, particles, datadir)
        out_fl = os.path.join(
            outdir, f"ctf_{os.path.splitext(os.path.basename(particles.path))[0]}.pkl"
        )
        parser = argparse.ArgumentParser()
        parse_ctf_star.add_args(parser)
        args = parser.parse_args(
            [particles.path, "-o", str(out_fl), "-D", "30", "--Apix", "1"]
        )
        parse_ctf_star.main(args)

        with open(out_fl, "rb") as f:
            out_ctf = pickle.load(f)

        assert out_ctf.shape == (1000, 9)
        assert np.allclose(out_ctf[:, 0], 30)  # D
        assert np.allclose(out_ctf[:, 1], 1.0)  # Apix

    @pytest.mark.parametrize(
        "indices", [None, "first-100", "random-100"], indirect=True
    )
    @pytest.mark.parametrize("poses", ["toy-poses", None], indirect=True)
    def test_write_star_from_mrcs(
        self, tmpdir_factory, particles, datadir, indices, poses
    ):
        outdir = self.get_outdir(tmpdir_factory, particles, datadir)
        out_fl = Path(outdir, f"written_{indices.label}_{poses.label}.star")
        parsed_ctf = Path(
            outdir, f"ctf_{os.path.splitext(os.path.basename(particles.path))[0]}.pkl"
        )
        assert parsed_ctf.exists(), "Upstream tests have failed!"

        parser = argparse.ArgumentParser()
        write_star.add_args(parser)
        args = [
            os.path.join(pytest.DATADIR, "toy_projections.mrcs"),
            "--ctf",
            str(parsed_ctf),
            "-o",
            str(out_fl),
        ]
        if indices.path is not None:
            args += ["--ind", indices.path]
        if poses.path is not None:
            args += ["--poses", poses.path]

        write_star.main(parser.parse_args(args))
        data = ImageSource.from_file(
            str(out_fl), lazy=False, datadir=datadir.path
        ).images()
        assert isinstance(data, torch.Tensor)
        assert data.shape == (1000 if indices.path is None else 100, 30, 30)

    @pytest.mark.parametrize("indices", [None], indirect=True)
    @pytest.mark.parametrize("poses", ["toy-poses"], indirect=True)
    def test_parse_pose(self, tmpdir_factory, particles, datadir, indices, poses):
        outdir = self.get_outdir(tmpdir_factory, particles, datadir)
        star_in = Path(outdir, f"written_{indices.label}_{poses.label}.star")
        out_fl = Path(outdir, "test_pose.pkl")
        assert star_in.exists(), "Upstream tests have failed!"

        parser = argparse.ArgumentParser()
        parse_pose_star.add_args(parser)
        parse_pose_star.main(
            parser.parse_args([str(star_in), "-D", "30", "-o", str(out_fl)])
        )

        out_poses = cryodrgn.utils.load_pkl(str(out_fl))
        assert isinstance(out_poses, tuple)
        assert len(out_poses) == 2
        assert isinstance(out_poses[0], np.ndarray)
        assert isinstance(out_poses[1], np.ndarray)
        assert out_poses[0].shape == (1000, 3, 3)
        assert out_poses[1].shape == (1000, 2)

        old_poses = cryodrgn.utils.load_pkl(poses.path)
        assert np.allclose(old_poses[0], out_poses[0], atol=1e-5)
        assert np.allclose(old_poses[1], out_poses[1], atol=1e-5)

    @pytest.mark.parametrize(
        "indices", [None, "first-100", "random-100"], indirect=True
    )
    @pytest.mark.parametrize("downsample_dim, chunk_size", [(28, 80), (14, 100)])
    def test_downsample_and_from_txt(
        self,
        tmpdir_factory,
        particles,
        datadir,
        downsample_dim,
        chunk_size,
        indices,
    ):
        outdir = self.get_outdir(tmpdir_factory, particles, datadir)
        out_mrcs = outdir / f"downsampled_{downsample_dim}.{chunk_size}.mrcs"
        out_txt = outdir / f"downsampled_{downsample_dim}.{chunk_size}.txt"

        parser = argparse.ArgumentParser()
        downsample.add_args(parser)
        args = parser.parse_args(
            [
                particles.path,
                "-D",
                str(downsample_dim),
                "--chunk",
                str(chunk_size),
                "--datadir",
                datadir.path,
                "-o",
                str(out_mrcs),
            ]
        )
        downsample.main(args)

        out_star = str(outdir / f"downsampled_{downsample_dim}.{chunk_size}.star")
        if indices.path is not None:
            out_star = out_star[:-5] + f"_{indices.label}.star"
        parsed_ctf = Path(
            outdir, f"ctf_{os.path.splitext(os.path.basename(particles.path))[0]}.pkl"
        )

        parser = argparse.ArgumentParser()
        write_star.add_args(parser)
        args = [str(out_txt), "--ctf", str(parsed_ctf), "-o", out_star]
        if indices.path is not None:
            args += ["--ind", indices.path]

        write_star.main(parser.parse_args(args))
        data = ImageSource.from_file(
            str(out_star), lazy=False, datadir=str(outdir)
        ).images()
        assert isinstance(data, torch.Tensor)
        assert data.shape == (
            1000 if indices.path is None else 100,
            downsample_dim,
            downsample_dim,
        )

    # NOTE: these must be a subset of the parameters
    #       used in `test_downsample_and_from_txt()` above to get input .txt particles!
    @pytest.mark.parametrize("downsample_dim, chunk_size", [(28, 80), (14, 100)])
    @pytest.mark.parametrize(
        "poses, ctf", [("toy-poses", "CTF-Test"), ("toy-poses", None)], indirect=True
    )
    @pytest.mark.parametrize("indices", [None, "random-100"], indirect=True)
    def test_backproject_from_downsample_txt(
        self,
        tmpdir_factory,
        particles,
        datadir,
        downsample_dim,
        chunk_size,
        poses,
        ctf,
        indices,
    ):
        """Use chunked downsampled .txt particle stack as input for backprojection."""
        outdir = self.get_outdir(tmpdir_factory, particles, datadir)
        out_txt = outdir / f"downsampled_{downsample_dim}.{chunk_size}.txt"
        outpath = outdir / f"backproj_{downsample_dim}.{chunk_size}"

        args = [str(out_txt), "--poses", poses.path]
        if indices.path is not None:
            outpath_parts = list(outpath.parts)
            outpath_parts[-1] = "_".join([outpath_parts[-1], indices.label])
            outpath = Path(*outpath_parts)
            args += ["--ind", indices.path]

        args += ["-o", str(outpath)]
        if ctf.path is not None:
            args += ["--ctf", ctf.path]
        if datadir.path is not None:
            args += ["--datadir", datadir.path]

        parser = argparse.ArgumentParser()
        backproject_voxel.add_args(parser)
        backproject_voxel.main(parser.parse_args(args))

        assert (outpath / "backproject.mrc").exists(), "Missing backprojected volume!"
        assert (outpath / "half_map_a.mrc").exists(), "Missing half-map volume!"
        assert (outpath / "half_map_b.mrc").exists(), "Missing half-map volume!"
        assert (outpath / "fsc-plot.png").exists(), "Missing plot of FSC curve!"
        assert (outpath / "fsc-vals.txt").exists(), "Missing FSC curve values!"

    # NOTE: these must be a subset of the parameters
    #       used in `test_downsample_and_from_txt()` above to get input .txt particles!
    @pytest.mark.parametrize("downsample_dim, chunk_size", [(28, 80), (14, 100)])
    @pytest.mark.parametrize(
        "poses, ctf", [("toy-poses", "CTF-Test"), ("toy-poses", None)], indirect=True
    )
    @pytest.mark.parametrize("indices", [None, "random-100"], indirect=True)
    def test_backproject_from_downsample_star(
        self,
        tmpdir_factory,
        particles,
        datadir,
        downsample_dim,
        chunk_size,
        poses,
        ctf,
        indices,
    ):
        """Use chunked downsampled .txt particle stack as input for backprojection."""
        outdir = self.get_outdir(tmpdir_factory, particles, datadir)
        outpath = outdir / f"backproj_{downsample_dim}.{chunk_size}"
        out_star = outdir / f"downsampled_{downsample_dim}.{chunk_size}.star"

        if indices.path is not None:
            outpath_parts = list(outpath.parts)
            outpath_parts[-1] = "_".join([outpath_parts[-1], indices.label])
            outpath = Path(*outpath_parts)

        args = [str(out_star)]
        args += ["--poses", poses.path, "-o", str(outpath), "--datadir", str(outdir)]
        if indices.path is not None:
            args += ["--ind", indices.path]
        if ctf.path is not None:
            args += ["--ctf", ctf.path]

        parser = argparse.ArgumentParser()
        backproject_voxel.add_args(parser)
        backproject_voxel.main(parser.parse_args(args))

        assert (outpath / "backproject.mrc").exists(), "Missing backprojected volume!"
        assert (outpath / "half_map_a.mrc").exists(), "Missing half-map volume!"
        assert (outpath / "half_map_b.mrc").exists(), "Missing half-map volume!"
        assert (outpath / "fsc-plot.png").exists(), "Missing plot of FSC curve!"
        assert (outpath / "fsc-vals.txt").exists(), "Missing FSC curve values!"

    def test_clean(self, tmpdir_factory, particles, datadir):
        outdir = self.get_outdir(tmpdir_factory, particles, datadir)
        shutil.rmtree(outdir)


@pytest.mark.parametrize(
    "particles, poses, ctf", [("toy.mrcs", "toy-poses", "CTF-Test")], indirect=True
)
@pytest.mark.parametrize("dsample_dim", [8, 10])
@pytest.mark.parametrize("chunksize", [4, 7])
class TestBackprojectFromChunkedDownsampled:
    def get_outpaths(
        self, tmpdir_factory, particles, poses, ctf, dsample_dim, chunksize
    ):
        dirname = Path(
            "BackprojectChunked",
            "_".join(
                [
                    particles.label,
                    poses.label,
                    ctf.label,
                    str(dsample_dim),
                    str(chunksize),
                ]
            ),
        )
        odir = tmpdir_factory.getbasetemp() / dirname
        os.makedirs(odir, exist_ok=True)
        out_mrcs, out_txt = odir / "downsampled.mrcs", odir / "downsampled.txt"

        return odir, out_mrcs, out_txt

    def test_downsample_with_chunks(
        self, tmpdir_factory, particles, poses, ctf, dsample_dim, chunksize
    ):
        outdir, out_mrcs, out_txt = self.get_outpaths(
            tmpdir_factory, particles, poses, ctf, dsample_dim, chunksize
        )

        parser = argparse.ArgumentParser()
        downsample.add_args(parser)
        args = parser.parse_args(
            [
                particles.path,
                "-D",
                str(dsample_dim),
                "--chunk",
                str(chunksize),
                "-o",
                str(out_mrcs),
            ]
        )
        downsample.main(args)

        orig_imgs = ImageSource.from_file(particles.path, lazy=False)
        assert ImageSource.from_file(str(out_txt), lazy=False).shape == (
            orig_imgs.shape[0],
            dsample_dim,
            dsample_dim,
        )

    @pytest.mark.parametrize("indices", [None, "random-100"], indirect=True)
    def test_backprojection_from_chunks(
        self, tmpdir_factory, particles, poses, ctf, indices, dsample_dim, chunksize
    ):
        outdir, out_mrcs, out_txt = self.get_outpaths(
            tmpdir_factory, particles, poses, ctf, dsample_dim, chunksize
        )
        outpath = outdir / "bproj"

        args = [str(out_txt), "--poses", poses.path]
        if ctf.path is not None:
            args += ["--ctf", ctf.path]
        if indices.path is not None:
            args += ["--ind", indices.path]
            outpath_parts = list(outpath.parts())
            outpath_parts[-1] = "_".join([str(outpath_parts[-1]), indices.label])
            outpath = Path(*outpath_parts)
        args += ["-o", str(outpath)]

        parser = argparse.ArgumentParser()
        backproject_voxel.add_args(parser)
        backproject_voxel.main(parser.parse_args(args))

        assert (outpath / "backproject.mrc").exists(), "Missing backprojected volume!"
        assert (outpath / "half_map_a.mrc").exists(), "Missing half-map volume!"
        assert (outpath / "half_map_b.mrc").exists(), "Missing half-map volume!"
        assert (outpath / "fsc-plot.png").exists(), "Missing plot of FSC curve!"
        assert (outpath / "fsc-vals.txt").exists(), "Missing FSC curve values!"

    def test_clean(self, tmpdir_factory, particles, poses, ctf, dsample_dim, chunksize):
        outdir, _, _ = self.get_outpaths(
            tmpdir_factory, particles, poses, ctf, dsample_dim, chunksize
        )
        shutil.rmtree(outdir)
