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


@pytest.fixture
def input_cs_proj_dir():
    # All paths stored in the .cs file as 'blob/path' entries are relative w.r.t a 'project directory'
    #   that is implicitly set inside cryoSparc, but needs to be explicitly specified otherwise.
    # TODO: Until the original .mrc files corresponding to input_cs_al() are located and included in testing data,
    #   or tests regenerated using a new .cs file along with new .mrc files, this will return None
    return None


@pytest.mark.parametrize("particles", ["toy.mrcs-999"], indirect=True)
@pytest.mark.parametrize("datadir", ["default-datadir"], indirect=True)
def test_read_mrcs(particles, datadir):
    data = ImageSource.from_file(particles, lazy=False, datadir=datadir).images()
    assert isinstance(data, torch.Tensor)
    # We have total 1000 particles of size 30x30 to begin with
    assert data.shape == (1000, 30, 30)


@pytest.mark.parametrize("particles", ["toy.star-13"], indirect=True)
@pytest.mark.parametrize("datadir", ["default-datadir"], indirect=True)
def test_read_starfile(particles, datadir):
    data = ImageSource.from_file(particles, lazy=False, datadir=datadir).images()
    assert isinstance(data, torch.Tensor)
    # We have 13 particles in our starfile, of size 30x30 to begin with
    assert data.shape == (13, 30, 30)


@pytest.mark.parametrize("particles", ["toy.star", "toy.star-13"], indirect=True)
@pytest.mark.parametrize("datadir", ["default-datadir"], indirect=True)
def test_filter(outdir, particles, datadir):
    indices_pkl = os.path.join(outdir, "indices.pkl")
    with open(indices_pkl, "wb") as f:
        # 0-based indices into the input star file
        # Note that these indices are simply the 0-indexed row numbers in the starfile,
        # and have nothing to do with the index of the individual particle in the MRCS
        # file (e.g. 00042@mymrcs.mrcs)
        pickle.dump([11, 3, 2, 4], f)

    out_fl = os.path.join(outdir, "issue150_filtered.star")
    parser = argparse.ArgumentParser()
    filter_star.add_args(parser)
    filter_star.main(parser.parse_args([particles, "--ind", indices_pkl, "-o", out_fl]))

    data = ImageSource.from_file(out_fl, lazy=False, datadir=datadir).images()
    assert isinstance(data, torch.Tensor)
    assert data.shape == (4, 30, 30)


@pytest.mark.parametrize("datadir", ["default-datadir"], indirect=True)
class TestParseCTFWriteStar:
    @pytest.mark.parametrize("particles", ["toy.star", "toy.star-13"], indirect=True)
    def test_parse_ctf_star(self, outdir, particles, datadir):
        out_fl = os.path.join(
            outdir, f"ctf_{os.path.splitext(os.path.basename(particles))[0]}.pkl"
        )

        parser = argparse.ArgumentParser()
        parse_ctf_star.add_args(parser)
        args = parser.parse_args(
            [particles, "-o", out_fl, "-D", "300", "--Apix", "1.035"]
        )
        parse_ctf_star.main(args)

        with open(out_fl, "rb") as f:
            out_ctf = pickle.load(f)
            # The ctf pkl file has N rows and 9 columns
            #   D, Apix, _rlnDefocusU, _rlnDefocusV, _rlnDefocusAngle, _rlnVoltage,
            #   _rlnSphericalAberration, _rlnAmplitudeContrast, _rlnPhaseShift
            assert (
                out_ctf.shape
                == {False: (1000, 9), True: (13, 9)}[
                    "13" in os.path.basename(particles)
                ]
            )
            assert np.allclose(out_ctf[:, 0], 300)  # D
            assert np.allclose(out_ctf[:, 1], 1.035)  # Apix

    @pytest.mark.parametrize("particles", ["toy.mrcs-999"], indirect=True)
    def test_write_star_from_mrcs(self, outdir, particles, datadir):
        out_fl = os.path.join(outdir, "written.star")
        parsed_ctf = os.path.join(outdir, "ctf_toy_projections.pkl")
        assert os.path.exists(parsed_ctf), "Upstream tests have failed!"

        parser = argparse.ArgumentParser()
        write_star.add_args(parser)
        args = parser.parse_args(
            [particles, "--ctf", parsed_ctf, "-o", out_fl, "--full-path"]
        )
        write_star.main(args)

        data = ImageSource.from_file(out_fl, lazy=False, datadir=datadir).images()
        assert isinstance(data, torch.Tensor)
        assert data.shape == (1000, 30, 30)
        os.remove(out_fl)

    # TODO: create relion3.1 tests
    @pytest.mark.parametrize("particles", ["toy.star", "toy.star-13"], indirect=True)
    def test_write_star_relion30(self, outdir, particles, datadir):
        indices_pkl = os.path.join(outdir, "indices.pkl")
        out_fl = os.path.join(outdir, "issue150_written_rel30.star")
        with open(indices_pkl, "wb") as f:
            pickle.dump([11, 3, 2, 4], f)

        parser = argparse.ArgumentParser()
        write_star.add_args(parser)
        args = parser.parse_args(
            [
                particles,
                "-o",
                out_fl,
                "--ind",
                indices_pkl,
                "--full-path",
                "--relion30",
            ]
        )
        write_star.main(args)

        data = ImageSource.from_file(out_fl, lazy=False, datadir=datadir).images()
        assert isinstance(data, torch.Tensor)
        assert data.shape == (4, 30, 30)
        os.remove(out_fl)


@pytest.mark.parametrize("particles", ["cryosparc-all"], indirect=True)
@pytest.mark.xfail(reason="The source .mrcs file for the .cs file are not available")
def test_write_cs(outdir, particles, input_cs_proj_dir):
    # Test writing out a .cs file from an input .cs file, with filtering
    # write_cs can optionally filter the output based on provided indices,
    # so we'll use that here
    indices_pkl = os.path.join(outdir, "indices.pkl")
    out_fl = os.path.join(outdir, "cs_filtered.cs")
    with open(indices_pkl, "wb") as f:
        pickle.dump([11, 3, 2, 4], f)

    parser = argparse.ArgumentParser()
    write_cs.add_args(parser)
    args = parser.parse_args(
        [
            particles,
            "--datadir",
            input_cs_proj_dir,
            "-o",
            out_fl,
            "--ind",
            indices_pkl,
            "--full-path",
        ]
    )
    write_cs.main(args)
