import pytest
import os
import argparse
import torch
import numpy as np
from cryodrgn.source import ImageSource
from cryodrgn.commands import downsample


@pytest.mark.parametrize(
    "particles, datadir",
    [
        ("toy.mrcs", None),
        ("toy.txt", None),
        ("toy.star-13", "default-datadir"),
        ("toydatadir.star", "toy"),
        ("relion31.v2.star", None),
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
    def test_downsample(self, tmpdir, particles, datadir, downsample_dim):
        in_imgs = ImageSource.from_file(
            particles.path, datadir=datadir.path, lazy=False
        ).images()
        out_mrcs = os.path.join(tmpdir, "downsampled.mrcs")

        parser = argparse.ArgumentParser()
        downsample.add_args(parser)
        args = parser.parse_args(
            [
                particles.path,  # 13 particles
                "-D",
                str(downsample_dim),  # downsampled from 30x30
                "--datadir",
                datadir.path,  # If specified, prefixed to each _rlnImageName in starfile
                "-o",
                out_mrcs,
            ]
        )
        # Note - no filtering is possible in downsample currently
        downsample.main(args)

        out_imgs = ImageSource.from_file(out_mrcs, lazy=False).images()
        assert isinstance(out_imgs, torch.Tensor)
        assert out_imgs.shape == (in_imgs.shape[0], downsample_dim, downsample_dim)
        assert np.isclose(in_imgs.sum(), out_imgs.sum())

    @pytest.mark.parametrize("downsample_dim, chunk_size", [(8, 5), (12, 6), (16, 8)])
    def test_downsample_with_chunks(
        self, tmpdir, particles, datadir, downsample_dim, chunk_size
    ):
        out_mrcs = os.path.join(tmpdir, "downsampled.mrcs")
        in_imgs = ImageSource.from_file(
            particles.path, datadir=datadir.path, lazy=False
        ).images()

        parser = argparse.ArgumentParser()
        downsample.add_args(parser)
        args = parser.parse_args(
            [
                particles.path,  # 13 particles
                "-D",
                str(downsample_dim),
                "--datadir",
                datadir.path,  # If specified, prefixed to each _rlnImageName in starfile
                "--chunk",
                str(chunk_size),
                "-o",
                out_mrcs,
            ]
        )
        downsample.main(args)

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
@pytest.mark.parametrize("downsample_dim", [16, 8])
def test_downsample_dir(tmpdir, particles, datadir, downsample_dim):
    in_imgs = ImageSource.from_file(
        particles.path, datadir=datadir.path, lazy=False
    ).images()
    out_star = os.path.join(tmpdir, "downsampled.star")
    outdir = os.path.join(tmpdir, f"downsampled-toy.{downsample_dim}")

    parser = argparse.ArgumentParser()
    downsample.add_args(parser)
    args = parser.parse_args(
        [
            particles.path,  # 13 particles
            "-D",
            str(downsample_dim),
            "--datadir",
            datadir.path,  # If specified, prefixed to each _rlnImageName in starfile
            "-o",
            out_star,
            "--outdir",
            outdir,
        ]
    )
    downsample.main(args)

    out_imgs = ImageSource.from_file(out_star, datadir=outdir, lazy=False).images()
    assert out_imgs.shape == (1000, downsample_dim, downsample_dim)
    assert np.isclose(in_imgs.sum(), out_imgs.sum())
