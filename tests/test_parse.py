import pytest
import argparse
import os
import shutil

from cryodrgn.commands import (
    parse_ctf_csparc,
    parse_ctf_star,
    parse_pose_csparc,
    parse_pose_star,
    parse_star,
)
from cryodrgn.commands_utils import write_star
from cryodrgn.utils import assert_pkl_close


@pytest.fixture
def particles_starfile():
    return os.path.join(pytest.DATADIR, "FinalRefinement-OriginalParticles-PfCRT.star")


class TestCtfStar:
    def get_outdir(self, tmpdir_factory, resolution):
        dirname = os.path.join("CTFStar", f"res.{resolution}")
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)

        return odir

    @pytest.mark.parametrize("resolution", ["128", "300"])
    def test_parse(self, tmpdir_factory, particles_starfile, resolution):
        outdir = self.get_outdir(tmpdir_factory, resolution)
        pkl_out = os.path.join(outdir, "ctf.pkl")
        png_out = os.path.join(outdir, "ctf.png")
        parser = argparse.ArgumentParser()
        parse_ctf_star.add_args(parser)
        args = parser.parse_args(
            [
                particles_starfile,
                "-w",
                "0.1",
                "-D",
                resolution,
                "--Apix",
                "1.035",
                "-o",
                pkl_out,
                "--png",
                png_out,
            ]
        )
        parse_ctf_star.main(args)

    @pytest.mark.parametrize("resolution", ["300"])
    def test_fidelity(self, tmpdir_factory, particles_starfile, resolution):
        outdir = self.get_outdir(tmpdir_factory, resolution)
        pkl_out = os.path.join(outdir, "ctf.pkl")
        assert_pkl_close(pkl_out, os.path.join(pytest.DATADIR, "ctf1.pkl"))

    @pytest.mark.parametrize("resolution", ["128", "300"])
    def test_write_star_from_mrcs(self, tmpdir_factory, particles_starfile, resolution):
        outdir = self.get_outdir(tmpdir_factory, resolution)
        mrcs_fl = os.path.join(pytest.DATADIR, "hand.5.mrcs")
        ctf_out = os.path.join(outdir, "ctf.pkl")

        parser = argparse.ArgumentParser()
        write_star.add_args(parser)
        args = parser.parse_args(
            [mrcs_fl, "--ctf", ctf_out, "-o", os.path.join(outdir, "test5.star")]
        )
        write_star.main(args)
        args = parser.parse_args(
            [mrcs_fl, "--ctf", ctf_out, "-o", os.path.join(outdir, "test6.star")]
        )
        write_star.main(args)

        shutil.rmtree(outdir)


@pytest.mark.parametrize("particles", ["csparc_big"], indirect=True)
def test_parse_ctf_cs(tmpdir, particles):
    pkl_out = os.path.join(tmpdir, "ctf.pkl")
    png_out = os.path.join(tmpdir, "ctf.png")
    args = parse_ctf_csparc.add_args(argparse.ArgumentParser()).parse_args(
        [particles.path, "-o", pkl_out, "--png", png_out]
    )
    parse_ctf_csparc.main(args)

    assert_pkl_close(pkl_out, os.path.join(pytest.DATADIR, "ctf2.pkl"))


def test_parse_pose_star(tmpdir, particles_starfile):
    pkl_out = os.path.join(tmpdir, "pose.pkl")
    parser = argparse.ArgumentParser()
    parse_pose_star.add_args(parser)
    args = parser.parse_args([particles_starfile, "-D", "300", "-o", pkl_out])
    parse_pose_star.main(args)

    assert_pkl_close(pkl_out, os.path.join(pytest.DATADIR, "pose.star.pkl"))


@pytest.mark.parametrize("resolution", ["128", "300"])
def test_parse_star(tmpdir, particles_starfile, resolution):
    ctf_out = os.path.join(tmpdir, "ctf.pkl")
    png_out = os.path.join(tmpdir, "ctf.png")
    poses_out = os.path.join(tmpdir, "pose.pkl")

    parser = argparse.ArgumentParser()
    parse_star.add_args(parser)
    args = parser.parse_args(
        [
            particles_starfile,
            "-D",
            resolution,
            "--poses",
            poses_out,
            "--ctf",
            ctf_out,
            "--png",
            png_out,
            "--Apix",
            "1.035",
            "-w",
            "0.1",
        ]
    )
    parse_star.main(args)

    assert os.path.isfile(ctf_out), "Missing CTF output file!"
    assert os.path.isfile(png_out), "Missing CTF plot file!"
    assert os.path.isfile(poses_out), "Missing pose output file!"
    if resolution == "300":
        assert_pkl_close(poses_out, os.path.join(pytest.DATADIR, "pose.star.pkl"))
        assert_pkl_close(ctf_out, os.path.join(pytest.DATADIR, "ctf1.pkl"))


@pytest.mark.parametrize("particles", ["csparc_big"], indirect=True)
def test_parse_pose_cs(tmpdir, particles):
    pkl_out = os.path.join(tmpdir, "pose.pkl")
    parser = argparse.ArgumentParser()
    parse_pose_csparc.add_args(parser)
    args = parser.parse_args([particles.path, "-D", "180", "-o", pkl_out])
    parse_pose_csparc.main(args)

    assert_pkl_close(pkl_out, os.path.join(pytest.DATADIR, "pose.cs.pkl"))
