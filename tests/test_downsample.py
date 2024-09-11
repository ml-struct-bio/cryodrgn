import pytest
import os
import shutil
import argparse
import torch
import numpy as np
from cryodrgn.source import ImageSource
from cryodrgn.commands import downsample
from cryodrgn.utils import load_pkl
from cryodrgn.mrcfile import parse_mrc, write_mrc
from cryodrgn.starfile import parse_star, write_star


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


@pytest.mark.parametrize(
    "particles, datadir", [("toydatadir.star", "toy")], indirect=True
)
@pytest.mark.parametrize(
    "newdatadir", [False, True], ids=["without-newdir", "with-newdir"]
)
def test_difficult_directory(tmpdir, particles, datadir, newdatadir):
    """Create a hard .star file setup with colliding file names between subfolders."""

    newdir = os.path.join(tmpdir, "toy-new") if newdatadir else tmpdir
    newdatadir_path = newdir if newdatadir else None
    xdir = os.path.join(newdir, "particles_x")
    ydir = os.path.join(newdir, "particles_y")
    os.makedirs(xdir, exist_ok=True)
    os.makedirs(ydir, exist_ok=True)

    shutil.copy(os.path.join(datadir.path, "toy_images_a.mrcs"), xdir)
    shutil.copy(os.path.join(datadir.path, "toy_images_b.mrcs"), xdir)
    shutil.copy(os.path.join(datadir.path, "toy_images_a.mrcs"), ydir)

    ymrcs, yheader = parse_mrc(os.path.join(ydir, "toy_images_a.mrcs"))
    write_mrc(os.path.join(ydir, "toy_images_a.mrcs"), ymrcs * -1, yheader)

    newstar, newoptics = parse_star(particles.path)
    new_imgnames = list()
    for i, imgname in enumerate(newstar["_rlnImageName"]):
        idx, filename = imgname.split("@")
        if filename == "toy_images_a.mrcs":
            if np.random.random() < 0.3:
                newfile = os.path.join("particles_y", "toy_images_a.mrcs")
            else:
                newfile = os.path.join("particles_x", "toy_images_a.mrcs")
        elif filename == "toy_images_b.mrcs":
            newfile = os.path.join("particles_x", "toy_images_b.mrcs")
        else:
            raise ValueError(filename)

        new_imgnames.append("@".join([idx, newfile]))

    newstar["_rlnImageName"] = new_imgnames
    newstarfile = os.path.join(tmpdir, "new-stack.star")
    write_star(newstarfile, newstar, newoptics)
    src = ImageSource.from_file(newstarfile, datadir=newdatadir_path, lazy=False)
    assert src.shape == (1000, 30, 30)

    out_star = os.path.join(tmpdir, "downsampled.star")
    parser = argparse.ArgumentParser()
    downsample.add_args(parser)
    newoutdir = os.path.join(tmpdir, "downsampled")
    args = [newstarfile, "-D", "12", "-o", out_star, "--outdir", newoutdir]
    if newdatadir:
        args += ["--datadir", newdir]
    downsample.main(parser.parse_args(args))

    in_imgs = src.images()
    out_imgs = ImageSource.from_file(out_star, datadir=newoutdir, lazy=False).images()
    assert out_imgs.shape == (src.n, 12, 12)
    assert np.isclose(in_imgs.sum(), out_imgs.sum())
