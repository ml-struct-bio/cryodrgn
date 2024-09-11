import pandas as pd
import pytest
import argparse
import os
import shutil
import pickle
import numpy as np
import torch
from itertools import product
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
    data = ImageSource.from_file(
        particles.path, lazy=False, datadir=datadir.path
    ).images()
    assert isinstance(data, torch.Tensor)
    # We have total 1000 particles of size 30x30 to begin with
    assert data.shape == (1000, 30, 30)


@pytest.mark.parametrize("particles", ["toy.star-13"], indirect=True)
@pytest.mark.parametrize("datadir", ["default-datadir"], indirect=True)
def test_read_starfile(particles, datadir):
    data = ImageSource.from_file(
        particles.path, lazy=False, datadir=datadir.path
    ).images()
    assert isinstance(data, torch.Tensor)
    # We have 13 particles in our starfile, of size 30x30 to begin with
    assert data.shape == (13, 30, 30)


@pytest.mark.parametrize("particles", ["toy.star", "toy.star-13"], indirect=True)
@pytest.mark.parametrize("datadir", ["default-datadir"], indirect=True)
@pytest.mark.parametrize(
    "index_pair",
    [([11, 3, 2, 4], [1, 2, 3]), ([5, 8, 11], [0, 7, 10])],
    ids=("inds1", "inds2"),
)
def test_filter(tmpdir, particles, datadir, index_pair):
    indices_pkl1 = os.path.join(tmpdir, "indices1.pkl")
    indices_pkl2 = os.path.join(tmpdir, "indices2.pkl")

    with open(indices_pkl1, "wb") as f:
        pickle.dump(index_pair[0], f)
    with open(indices_pkl2, "wb") as f:
        pickle.dump(index_pair[1], f)

    out_fl = os.path.join(tmpdir, "issue150_filtered.star")
    parser = argparse.ArgumentParser()
    filter_star.add_args(parser)
    filter_star.main(
        parser.parse_args([particles.path, "--ind", indices_pkl1, "-o", out_fl])
    )

    data1 = ImageSource.from_file(out_fl, lazy=False, datadir=datadir.path).images()
    assert isinstance(data1, torch.Tensor)
    assert data1.shape == (len(index_pair[0]), 30, 30)
    os.remove(out_fl)

    filter_star.main(
        parser.parse_args([particles.path, "--ind", indices_pkl2, "-o", out_fl])
    )

    data2 = ImageSource.from_file(out_fl, lazy=False, datadir=datadir.path).images()
    assert isinstance(data2, torch.Tensor)
    assert data2.shape == (len(index_pair[1]), 30, 30)
    os.remove(out_fl)

    for (i, ind1), (j, ind2) in product(
        enumerate(index_pair[0]), enumerate(index_pair[1])
    ):
        if ind1 == ind2:
            assert np.allclose(data1[i, ...], data2[j, ...])


@pytest.mark.parametrize("datadir", ["default-datadir"], indirect=True)
@pytest.mark.parametrize("particles", ["tilts.star"], indirect=True)
class TestFilterStar:
    def test_filter_with_indices(self, tmpdir, particles, datadir):
        indices_pkl = os.path.join(tmpdir, "indices.pkl")
        # 0-based indices into the input star file
        # Note that these indices are simply the 0-indexed row numbers in the starfile,
        # and have nothing to do with the index of the individual particle in the MRCS
        # file (e.g. 00042@mymrcs.mrcs)
        with open(indices_pkl, "wb") as f:
            pickle.dump([1, 3, 4], f)

        out_fl = os.path.join(tmpdir, "tilts_filtered.star")
        parser = argparse.ArgumentParser()
        filter_star.add_args(parser)
        filter_star.main(
            parser.parse_args([particles.path, "--ind", indices_pkl, "-o", out_fl])
        )

        data = ImageSource.from_file(out_fl, lazy=False, datadir=datadir.path).images()
        assert isinstance(data, torch.Tensor)
        assert data.shape == (3, 64, 64)
        os.remove(out_fl)

        filter_star.main(
            parser.parse_args(
                [particles.path, "--ind", indices_pkl, "-o", out_fl, "--et"]
            )
        )

        data = ImageSource.from_file(out_fl, lazy=False, datadir=datadir.path).images()
        assert isinstance(data, torch.Tensor)
        assert data.shape == (39 * 3, 64, 64)
        os.remove(out_fl)

        with open(indices_pkl, "wb") as f:
            pickle.dump([1, 3, 5], f)
        filter_star.main(
            parser.parse_args(
                [particles.path, "--ind", indices_pkl, "-o", out_fl, "--et"]
            )
        )

        data = ImageSource.from_file(out_fl, lazy=False, datadir=datadir.path).images()
        assert isinstance(data, torch.Tensor)
        assert data.shape == (39 + 39 + 24, 64, 64)
        os.remove(out_fl)

    @pytest.mark.parametrize("indices", ["just-4", "just-5"], indirect=True)
    def test_filter_with_separate_files(self, tmpdir, particles, indices, datadir):
        out_dir = os.path.join(tmpdir, "tilts_filt-names")
        parser = argparse.ArgumentParser()
        filter_star.add_args(parser)

        out_fl = os.path.join(tmpdir, "tilts_filtered.star")
        args = [particles.path, "--ind", indices.path, "-o", out_fl, "--et"]
        filter_star.main(parser.parse_args(args))

        data = ImageSource.from_file(out_fl, lazy=False, datadir=datadir.path).images()
        assert isinstance(data, torch.Tensor)

        all_n = data.shape[0]
        os.remove(out_fl)

        args = [
            particles.path,
            "--ind",
            indices.path,
            "-o",
            out_dir,
            "--micrograph-files",
            "--et",
        ]
        filter_star.main(parser.parse_args(args))

        split_n = 0
        for out_fl in os.listdir(out_dir):
            data = ImageSource.from_file(
                os.path.join(out_dir, out_fl), lazy=False, datadir=datadir.path
            ).images()
            assert isinstance(data, torch.Tensor)
            split_n += data.shape[0]

        assert split_n == all_n
        shutil.rmtree(out_dir)


@pytest.mark.parametrize("datadir", ["default-datadir"], indirect=True)
class TestParseCTFWriteStar:
    def get_outdir(self, tmpdir_factory, particles, datadir):
        dirname = os.path.join("ParseCTFWriteStar", particles.label, datadir.label)
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)

        return odir

    @pytest.mark.parametrize(
        "particles, resolution, apix",
        [
            ("toy.star", 300, 1.035),
            ("toy.star-13", 160, 2.86),
            pytest.param(
                "toy-star",
                200,
                None,
                marks=pytest.mark.xfail(raises=ValueError, reason="no A/px available!"),
            ),
            pytest.param(
                "toy-star-13",
                None,
                2.86,
                marks=pytest.mark.xfail(
                    raises=ValueError, reason="no resolution available!"
                ),
            ),
        ],
        indirect=["particles"],
    )
    def test_parse_ctf_star(self, tmpdir_factory, particles, datadir, resolution, apix):
        outdir = self.get_outdir(tmpdir_factory, particles, datadir)
        out_lbl = os.path.splitext(os.path.basename(particles.path))[0]
        if resolution is not None:
            out_lbl += "_{}".format(resolution)
        if apix is not None:
            out_lbl += "_{}".format(apix)

        out_fl = os.path.join(outdir, f"ctf_{out_lbl}.pkl")
        in_src = ImageSource.from_file(particles.path, lazy=False, datadir=datadir.path)
        parser = argparse.ArgumentParser()
        parse_ctf_star.add_args(parser)
        args = [particles.path, "-o", out_fl]
        if resolution is not None:
            args += ["-D", str(resolution)]
        if apix is not None:
            args += ["--Apix", str(apix)]

        parse_ctf_star.main(parser.parse_args(args))

        # The ctf pkl file has N rows and 9 columns
        #   D, Apix, _rlnDefocusU, _rlnDefocusV, _rlnDefocusAngle, _rlnVoltage,
        #   _rlnSphericalAberration, _rlnAmplitudeContrast, _rlnPhaseShift
        with open(out_fl, "rb") as f:
            out_ctf = pickle.load(f)

        assert out_ctf.shape == (in_src.n, 9)
        assert np.allclose(out_ctf[:, 0], resolution)  # D
        assert np.allclose(out_ctf[:, 1], apix)  # Apix

    @pytest.mark.parametrize("particles", ["toy.star"], indirect=True)
    def test_write_star_from_mrcs(self, tmpdir_factory, particles, datadir):
        outdir = self.get_outdir(tmpdir_factory, particles, datadir)
        out_fl = os.path.join(outdir, "written.star")
        parsed_ctf = os.path.join(outdir, "ctf_toy_projections_300_1.035.pkl")
        assert os.path.exists(parsed_ctf), "Upstream tests have failed!"

        in_imgs = ImageSource.from_file(
            particles.path, lazy=False, datadir=datadir.path
        ).images()

        args = [
            os.path.join(pytest.DATADIR, "toy_projections_0-999.mrcs"),
            "--ctf",
            parsed_ctf,
            "-o",
            out_fl,
            "--full-path",
        ]

        parser = argparse.ArgumentParser()
        write_star.add_args(parser)
        write_star.main(parser.parse_args(args))
        out_imgs = ImageSource.from_file(
            out_fl, lazy=False, datadir=datadir.path
        ).images()

        assert isinstance(out_imgs, torch.Tensor)
        assert out_imgs.shape == in_imgs.shape
        os.remove(out_fl)

    @pytest.mark.parametrize("particles", ["toy.star", "toy.star-13"], indirect=True)
    @pytest.mark.parametrize(
        "use_relion30", [False, True], ids=("relion3.1", "relion3.0")
    )
    def test_write_filter_star(self, tmpdir_factory, particles, datadir, use_relion30):
        outdir = self.get_outdir(tmpdir_factory, particles, datadir)
        indices_pkl = os.path.join(outdir, "indices.pkl")
        out_fl = os.path.join(outdir, "issue150_written_rel30.star")

        in_src = ImageSource.from_file(particles.path, lazy=False, datadir=datadir.path)
        with open(indices_pkl, "wb") as f:
            pickle.dump([11, 3, 2, 4], f)

        parser = argparse.ArgumentParser()
        write_star.add_args(parser)
        args = [particles.path, "-o", out_fl, "--ind", indices_pkl, "--full-path"]
        if use_relion30:
            args += ["--relion30"]

        write_star.main(parser.parse_args(args))
        out_src = ImageSource.from_file(out_fl, lazy=False, datadir=datadir.path)
        assert isinstance(out_src.df, pd.DataFrame)
        assert out_src.df.shape == (4, 12)
        out_imgs = out_src.images()
        assert isinstance(out_imgs, torch.Tensor)
        assert out_imgs.shape == (4, in_src.D, in_src.D)

        if use_relion30 or in_src.data_optics is None:
            assert out_src.data_optics is None
        else:
            assert out_src.data_optics.shape == in_src.data_optics.shape

        os.remove(out_fl)


@pytest.mark.parametrize("particles", ["cryosparc-all"], indirect=True)
@pytest.mark.xfail(reason="The source .mrcs file for the .cs file are not available")
def test_write_cs(tmpdir, particles, input_cs_proj_dir):
    # Test writing out a .cs file from an input .cs file, with filtering
    # write_cs can optionally filter the output based on provided indices,
    # so we'll use that here
    indices_pkl = os.path.join(tmpdir, "indices.pkl")
    out_fl = os.path.join(tmpdir, "cs_filtered.cs")
    with open(indices_pkl, "wb") as f:
        pickle.dump([11, 3, 2, 4], f)

    parser = argparse.ArgumentParser()
    write_cs.add_args(parser)
    args = parser.parse_args(
        [
            particles.path,
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
