import argparse
import os
import os.path
import pickle
import numpy as np
import torch
import pytest
from cryodrgn.source import ImageSource
from cryodrgn.commands import parse_ctf_star
from cryodrgn.commands_utils import filter_star, write_cs, write_star

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture
def input_star_all():
    return f"{DATA_FOLDER}/toy_projections.star"


@pytest.fixture
def input_cs_all():
    # cs file with 2019 particle entries
    return f"{DATA_FOLDER}/cryosparc_P12_J24_001_particles.cs"


@pytest.fixture
def input_cs_proj_dir():
    # All paths stored in the .cs file as 'blob/path' entries are relative w.r.t a 'project directory'
    #   that is implicitly set inside cryoSparc, but needs to be explicitly specified otherwise.
    # TODO: Until the original .mrc files corresponding to input_cs_al() are located and included in testing data,
    #   or tests regenerated using a new .cs file along with new .mrc files, this will return None
    return None


@pytest.fixture
def input_star():
    # starfile refers to 13 random shuffled indices into mrcs between 0-1000
    return f"{DATA_FOLDER}/toy_projections_13.star"


@pytest.fixture
def input_mrcs():
    # 1000 particle mrcs
    return f"{DATA_FOLDER}/toy_projections_0-999.mrcs"


def test_read_mrcs(input_mrcs):
    data = ImageSource.from_file(input_mrcs, lazy=False, datadir=DATA_FOLDER).images()
    assert isinstance(data, torch.Tensor)
    # We have total 1000 particles of size 30x30 to begin with
    assert data.shape == (1000, 30, 30)


def test_read_starfile(input_star):
    data = ImageSource.from_file(input_star, lazy=False, datadir=DATA_FOLDER).images()
    assert isinstance(data, torch.Tensor)
    # We have 13 particles in our starfile, of size 30x30 to begin with
    assert data.shape == (13, 30, 30)


def test_filter(input_star):
    os.makedirs("output", exist_ok=True)

    indices_pkl = "output/indices.pkl"
    with open(indices_pkl, "wb") as f:
        # 0-based indices into the input star file
        # Note that these indices are simply the 0-indexed row numbers in the starfile, and have nothing to do with
        #   the index of the individual particle in the MRCS file (e.g. 00042@mymrcs.mrcs)
        pickle.dump([11, 3, 2, 4], f)

    args = filter_star.add_args(argparse.ArgumentParser()).parse_args(
        [
            input_star,
            "--ind",
            indices_pkl,
            "-o",
            "output/issue150_filtered.star",
        ]  # 13 particles
    )
    filter_star.main(args)

    data = ImageSource.from_file(
        "output/issue150_filtered.star", lazy=False, datadir=DATA_FOLDER
    ).images()
    assert isinstance(data, torch.Tensor)
    assert data.shape == (4, 30, 30)


def test_parse_ctf_star(input_star):
    os.makedirs("output", exist_ok=True)

    args = parse_ctf_star.add_args(argparse.ArgumentParser()).parse_args(
        [
            input_star,  # 13 particles
            "-o",
            "output/ctf.pkl",
            "-D",
            "300",
            "--Apix",
            "1.035",
        ]
    )
    parse_ctf_star.main(args)

    with open("output/ctf.pkl", "rb") as f:
        ctf = pickle.load(f)
        # The ctf pkl file has N rows and 9 columns
        #   D, Apix, _rlnDefocusU, _rlnDefocusV, _rlnDefocusAngle, _rlnVoltage, _rlnSphericalAberration,
        #   _rlnAmplitudeContrast, _rlnPhaseShift
        assert ctf.shape == (13, 9)
        assert np.allclose(ctf[:, 0], 300)  # D
        assert np.allclose(ctf[:, 1], 1.035)  # Apix


def test_write_star(input_mrcs, input_star_all):
    os.makedirs("output", exist_ok=True)

    # The CTF pkl file is required for write_star, for ALL particles in the .mrcs we provide it, so we create it here
    args = parse_ctf_star.add_args(argparse.ArgumentParser()).parse_args(
        [
            input_star_all,  # 1000 particles
            "-o",
            "output/ctf.pkl",
            "-D",
            "300",
            "--Apix",
            "1.035",
        ]
    )
    parse_ctf_star.main(args)
    with open("output/ctf.pkl", "rb") as f:
        ctf = pickle.load(f)
        assert ctf.shape == (1000, 9)

    # write_star can optionally filter the output based on provided indices, so we'll use that here
    indices_pkl = "output/indices.pkl"
    with open(indices_pkl, "wb") as f:
        pickle.dump([11, 3, 2, 4], f)

    args = write_star.add_args(argparse.ArgumentParser()).parse_args(
        [
            input_mrcs,  # 1000 particles
            "--ctf",
            "output/ctf.pkl",  # 1000 ctf params
            "-o",
            "output/issue150_written1.star",
            "--ind",
            "output/indices.pkl",
            "--full-path",
        ]
    )
    write_star.main(args)


def test_write_star2(input_mrcs, input_star_all):
    os.makedirs("output", exist_ok=True)

    # write_star can optionally filter the output based on provided indices, so we'll use that here
    indices_pkl = "output/indices.pkl"
    with open(indices_pkl, "wb") as f:
        pickle.dump([11, 3, 2, 4], f)

    args = write_star.add_args(argparse.ArgumentParser()).parse_args(
        [
            input_star_all,  # 1000 particles
            "-o",
            "output/issue150_written2.star",
            "--ind",
            "output/indices.pkl",
            "--full-path",
        ]
    )
    write_star.main(args)


@pytest.mark.xfail(reason="The source .mrcs file for the .cs file are not available")
def test_write_cs(input_cs_all, input_cs_proj_dir):
    # Test writing out a .cs file from an input .cs file, with filtering
    os.makedirs("output", exist_ok=True)

    # write_cs can optionally filter the output based on provided indices, so we'll use that here
    indices_pkl = "output/indices.pkl"
    with open(indices_pkl, "wb") as f:
        pickle.dump([11, 3, 2, 4], f)

    args = write_cs.add_args(argparse.ArgumentParser()).parse_args(
        [
            input_cs_all,  # 1000 particles
            "--datadir",
            input_cs_proj_dir,
            "-o",
            "output/cs_filtered.cs",
            "--ind",
            "output/indices.pkl",
            "--full-path",
        ]
    )
    write_cs.main(args)
