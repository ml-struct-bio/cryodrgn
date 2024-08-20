import pytest
import os
import argparse
import torch
import numpy as np
from cryodrgn.source import ImageSource
from cryodrgn.commands import downsample
from cryodrgn.utils import load_pkl


@pytest.mark.parametrize(
    "particles, datadir, indices",
    [
        ("toy.mrcs", None, None),
        ("toy.mrcs", None, "just-5"),
        ("toy.txt", None, None),
        ("toy.txt", None, "just-5"),
        ("toy.star-13", "default-datadir", None),
        ("toy.star-13", "default-datadir", "just-4"),
        ("toydatadir.star", "toy", None),
        ("relion31.v2.star", None, None),
    ],
    indirect=True,
)
class TestDownsampleToMRCS:
    @pytest.mark.parametrize(
        "downsample_dim",
        [
            28,
            pytest.param(
                15,
                marks=pytest.mark.xfail(raises=ValueError, reason="odd box size"),
            ),
            10,
        ],
    )
    def test_downsample(self, tmpdir, particles, datadir, indices, downsample_dim):
        use_ind = load_pkl(indices.path) if indices.path is not None else None
        in_imgs = ImageSource.from_file(
            particles.path, datadir=datadir.path, lazy=False, indices=use_ind
        ).images()
        out_mrcs = os.path.join(tmpdir, "downsampled.mrcs")

        parser = argparse.ArgumentParser()
        downsample.add_args(parser)
        args = [
            particles.path,
            "-D",
            str(downsample_dim),  # downsampled from 30x30
            "--datadir",
            datadir.path,  # If specified, prefixed to each _rlnImageName in starfile
            "-o",
            out_mrcs,
        ]
        if indices.path is not None:
            args += ["--ind", indices.path]

        downsample.main(parser.parse_args(args))
        out_imgs = ImageSource.from_file(out_mrcs, lazy=False).images()
        assert isinstance(out_imgs, torch.Tensor)
        assert out_imgs.shape == (in_imgs.shape[0], downsample_dim, downsample_dim)
        assert np.isclose(in_imgs.sum(), out_imgs.sum())

    @pytest.mark.parametrize("downsample_dim, chunk_size", [(8, 5), (12, 6), (16, 8)])
    def test_downsample_with_chunks(
        self, tmpdir, particles, datadir, downsample_dim, chunk_size, indices
    ):
        out_mrcs = os.path.join(tmpdir, "downsampled.mrcs")
        use_ind = load_pkl(indices.path) if indices.path is not None else None
        in_imgs = ImageSource.from_file(
            particles.path, datadir=datadir.path, lazy=False, indices=use_ind
        ).images()

        parser = argparse.ArgumentParser()
        downsample.add_args(parser)
        args = [
            particles.path,
            "-D",
            str(downsample_dim),
            "--datadir",
            datadir.path,  # If specified, prefixed to each _rlnImageName in starfile
            "--chunk",
            str(chunk_size),
            "-o",
            out_mrcs,
        ]
        if indices.path is not None:
            args += ["--ind", indices.path]

        downsample.main(parser.parse_args(args))
        assert ImageSource.from_file(
            out_mrcs.replace(".mrcs", ".txt"), lazy=False
        ).shape == (in_imgs.shape[0], downsample_dim, downsample_dim)

        for i in range((in_imgs.shape[0] // chunk_size) + 4):
            if i < (in_imgs.shape[0] // chunk_size):
                assert ImageSource.from_file(
                    out_mrcs.replace(".mrcs", f".{i}.mrcs")
                ).shape == (
                    chunk_size,
                    downsample_dim,
                    downsample_dim,
                )
            elif (
                i == (in_imgs.shape[0] // chunk_size) and in_imgs.shape[0] % chunk_size
            ):
                assert ImageSource.from_file(
                    out_mrcs.replace(".mrcs", f".{i}.mrcs")
                ).shape == (
                    in_imgs.shape[0] % chunk_size,
                    downsample_dim,
                    downsample_dim,
                )
            else:
                assert not os.path.exists(out_mrcs.replace(".mrcs", f".{i}.mrcs"))


@pytest.mark.parametrize(
    "particles, datadir", [("toydatadir.star", "toy")], indirect=True
)
@pytest.mark.parametrize(
    "outdir", ["downsampled", None], ids=["with-outdir", "wo-outdir"]
)
@pytest.mark.parametrize("downsample_dim", [16, 8])
def test_downsample_starout(tmpdir, particles, datadir, outdir, downsample_dim):
    in_imgs = ImageSource.from_file(
        particles.path, datadir=datadir.path, lazy=False
    ).images()
    out_star = os.path.join(tmpdir, "downsampled.star")
    outdir = os.path.join(tmpdir, f"downsampled-toy.{downsample_dim}")

    parser = argparse.ArgumentParser()
    downsample.add_args(parser)
    args = [
        particles.path,
        "-D",
        str(downsample_dim),
        "--datadir",
        datadir.path,  # If specified, prefixed to each _rlnImageName in starfile
        "-o",
        out_star,
    ]
    if outdir is not None:
        outpath = f"{outdir}-{datadir.label}.{downsample_dim}"
        args += ["--outdir", outpath]

    downsample.main(parser.parse_args(args))
    out_imgs = ImageSource.from_file(out_star, datadir=outpath, lazy=False).images()
    assert out_imgs.shape == (in_imgs.shape[0], downsample_dim, downsample_dim)
    assert np.isclose(in_imgs.sum(), out_imgs.sum())


@pytest.mark.parametrize("particles", ["toy.txt"], indirect=True)
@pytest.mark.parametrize(
    "outdir", ["downsampled", None], ids=["with-outdir", "wo-outdir"]
)
@pytest.mark.parametrize("downsample_dim", [16, 8])
def test_downsample_txtout(tmpdir, particles, outdir, downsample_dim):
    in_imgs = ImageSource.from_file(particles.path, lazy=False).images()
    out_txt = os.path.join(tmpdir, "downsampled.txt")

    parser = argparse.ArgumentParser()
    downsample.add_args(parser)
    args = [
        particles.path,
        "-D",
        str(downsample_dim),
        "-o",
        out_txt,
    ]
    if outdir is not None:
        outpath = f"{outdir}-{particles.label}.{downsample_dim}"
        args += ["--outdir", outpath]

    downsample.main(parser.parse_args(args))
    out_imgs = ImageSource.from_file(out_txt, lazy=False).images()
    assert out_imgs.shape == (in_imgs.shape[0], downsample_dim, downsample_dim)
    assert np.isclose(in_imgs.sum(), out_imgs.sum())
