import pytest
import os
import argparse
import torch
from cryodrgn.source import ImageSource
from cryodrgn.commands import downsample


@pytest.mark.parametrize("particles", ["toy.star-13"], indirect=True)
@pytest.mark.parametrize("datadir", ["default-datadir"], indirect=True)
@pytest.mark.parametrize(
    "downsample_dim",
    [
        28,
        pytest.param(
            15,
            marks=pytest.mark.xfail(
                raises=AssertionError, reason="odd downsampled box size"
            ),
        ),
        10,
    ],
)
def test_downsample(tmpdir, particles, datadir, downsample_dim):
    out_mrcs = os.path.join(tmpdir, "downsampled.mrcs")
    args = downsample.add_args(argparse.ArgumentParser()).parse_args(
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

    output_data = ImageSource.from_file(out_mrcs, lazy=False).images()
    assert isinstance(output_data, torch.Tensor)
    assert output_data.shape == (13, downsample_dim, downsample_dim)


@pytest.mark.parametrize("particles", ["toy.star-13"], indirect=True)
@pytest.mark.parametrize("datadir", ["default-datadir"], indirect=True)
@pytest.mark.parametrize("downsample_dim", [18, 10])
@pytest.mark.parametrize("chunk_size", [5, 6, 8])
def test_downsample_with_chunks(tmpdir, particles, datadir, downsample_dim, chunk_size):
    out_mrcs = os.path.join(tmpdir, "downsampled.mrcs")
    args = downsample.add_args(argparse.ArgumentParser()).parse_args(
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
    ).shape == (13, downsample_dim, downsample_dim)

    for i in range((13 // chunk_size) + 4):
        if i < (13 // chunk_size):
            assert ImageSource.from_file(
                out_mrcs.replace(".mrcs", f".{i}.mrcs")
            ).shape == (
                chunk_size,
                downsample_dim,
                downsample_dim,
            )
        elif i == (13 // chunk_size):
            assert ImageSource.from_file(
                out_mrcs.replace(".mrcs", f".{i}.mrcs")
            ).shape == (
                13 % chunk_size,
                downsample_dim,
                downsample_dim,
            )
        else:
            assert not os.path.exists(out_mrcs.replace(".mrcs", f".{i}.mrcs"))
