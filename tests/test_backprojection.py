import pytest
import os
import shutil
import argparse
from cryodrgn.commands import backproject_voxel
from cryodrgn.commands_utils import plot_fsc


class TestBackprojection:
    def get_outdir(self, tmpdir_factory, particles, poses, ctf, indices, datadir):
        dirname = os.path.join(
            "Backprojection",
            particles.label,
            poses.label,
            ctf.label,
            indices.label,
            datadir.label,
        )
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)

        return odir

    # NOTE: these have to be the same as in `test_to_fsc()` below for proper cleanup!
    @pytest.mark.parametrize(
        "particles, poses, ctf, datadir",
        [
            ("toy.mrcs", "toy-poses", None, None),
            ("toy.mrcs", "toy-poses", "CTF-Test", None),
            ("toy.txt", "toy-poses", "CTF-Test", None),
            ("hand", "hand-rot", None, None),
            ("hand", "hand-poses", None, None),
            ("tilts.star", "tilt-poses", "CTF-Tilt", None),
            ("tilts.star", "tilt-poses", "CTF-Tilt", "default-datadir"),
        ],
        indirect=True,
    )
    @pytest.mark.parametrize("indices", [None, "just-5"], indirect=True)
    def test_train(self, tmpdir_factory, particles, poses, ctf, indices, datadir):
        outdir = self.get_outdir(
            tmpdir_factory, particles, poses, ctf, indices, datadir
        )
        args = [
            particles.path,
            "--poses",
            poses.path,
            "-o",
            os.path.join(outdir, "vol.mrc"),
        ]

        if ctf.path is not None:
            args += ["--ctf", ctf.path]
        if indices.path is not None:
            args += ["--ind", indices.path]
        if datadir.path is not None:
            args += ["--datadir", datadir.path]

        parser = argparse.ArgumentParser()
        backproject_voxel.add_args(parser)
        backproject_voxel.main(parser.parse_args(args))

        assert os.path.exists(os.path.join(outdir, "vol.mrc"))
        assert os.path.exists(os.path.join(outdir, "vol_fsc-plot.png"))
        assert os.path.exists(os.path.join(outdir, "vol_fsc-vals.txt"))
        assert os.path.exists(os.path.join(outdir, "vol_half-map1.mrc"))
        assert os.path.exists(os.path.join(outdir, "vol_half-map2.mrc"))

    @pytest.mark.parametrize(
        "particles, poses, ctf, datadir",
        [
            ("hand", "hand-rot", None, None),
            ("hand", "hand-poses", None, None),
        ],
        indirect=True,
    )
    @pytest.mark.parametrize(
        "indices",
        [
            None,
            pytest.param(
                "just-5",
                marks=pytest.mark.xfail(
                    raises=AssertionError, reason="not enough particles!"
                ),
            ),
        ],
        indirect=True,
    )
    def test_fidelity(self, tmpdir_factory, particles, poses, ctf, indices, datadir):
        outdir = self.get_outdir(
            tmpdir_factory, particles, poses, ctf, indices, datadir
        )
        with open(os.path.join(outdir, "vol_fsc-vals.txt"), "r") as f:
            pixres, fsc_val = f.readlines()[-1].strip().split()

        assert round(float(pixres), 3) == 0.484
        assert float(fsc_val) > 0.2

    # NOTE: these have to be the same as in `test_train()` above for proper cleanup!
    @pytest.mark.parametrize(
        "particles, poses, ctf, datadir",
        [
            ("toy.mrcs", "toy-poses", None, None),
            ("toy.mrcs", "toy-poses", "CTF-Test", None),
            ("toy.txt", "toy-poses", "CTF-Test", None),
            ("hand", "hand-rot", None, None),
            ("hand", "hand-poses", None, None),
            ("tilts.star", "tilt-poses", "CTF-Tilt", None),
            ("tilts.star", "tilt-poses", "CTF-Tilt", "default-datadir"),
        ],
        indirect=True,
    )
    @pytest.mark.parametrize("indices", [None, "just-5"], indirect=True)
    def test_to_fsc(self, tmpdir_factory, particles, poses, ctf, indices, datadir):
        """Export calculated FSC scores to plotting function."""
        outdir = self.get_outdir(
            tmpdir_factory, particles, poses, ctf, indices, datadir
        )
        args = [
            os.path.join(outdir, "vol_fsc-vals.txt"),
            "-o",
            os.path.join(outdir, "vol_fsc-plot2.png"),
        ]

        parser = argparse.ArgumentParser()
        plot_fsc.add_args(parser)
        plot_fsc.main(parser.parse_args(args))
        assert os.path.exists(os.path.join(outdir, "vol_fsc-plot2.png"))

        shutil.rmtree(outdir)


class TestTiltBackprojection:
    def get_outdir(
        self, tmpdir_factory, particles, poses, ctf, indices, datadir, ntilts
    ):
        dirname = os.path.join(
            "BackprojectionTilt",
            particles.label,
            poses.label,
            ctf.label,
            indices.label,
            datadir.label,
            f"ntilts.{ntilts}",
        )
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)

        return odir

    # NOTE: these have to be the same as in `test_to_fsc()` below for proper cleanup!
    @pytest.mark.parametrize(
        "particles, poses, ctf, datadir",
        [
            ("tilts.star", "tilt-poses", "CTF-Tilt", None),
            ("tilts.star", "tilt-poses", "CTF-Tilt", "default-datadir"),
        ],
        indirect=True,
    )
    @pytest.mark.parametrize("indices", [None, "just-5"], indirect=True)
    @pytest.mark.parametrize("ntilts", [None, 5, 20])
    def test_train(
        self, tmpdir_factory, particles, poses, ctf, indices, datadir, ntilts
    ):
        outdir = self.get_outdir(
            tmpdir_factory, particles, poses, ctf, indices, datadir, ntilts
        )
        args = [
            particles.path,
            "--poses",
            poses.path,
            "-o",
            os.path.join(outdir, "vol.mrc"),
            "--tilt",
            "--dose-per-tilt",
            "2.93",
        ]

        if ctf.path is not None:
            args += ["--ctf", ctf.path]
        if indices.path is not None:
            args += ["--ind", indices.path]
        if datadir.path is not None:
            args += ["--datadir", datadir.path]
        if ntilts is not None:
            args += ["--ntilts", str(ntilts)]

        parser = argparse.ArgumentParser()
        backproject_voxel.add_args(parser)
        backproject_voxel.main(parser.parse_args(args))

        assert os.path.exists(os.path.join(outdir, "vol.mrc"))
        assert os.path.exists(os.path.join(outdir, "vol_fsc-plot.png"))
        assert os.path.exists(os.path.join(outdir, "vol_fsc-vals.txt"))
        assert os.path.exists(os.path.join(outdir, "vol_half-map1.mrc"))
        assert os.path.exists(os.path.join(outdir, "vol_half-map2.mrc"))

    # NOTE: these have to be the same as in `test_train()` above for proper cleanup!
    @pytest.mark.parametrize(
        "particles, poses, ctf, datadir",
        [
            ("tilts.star", "tilt-poses", "CTF-Tilt", None),
            ("tilts.star", "tilt-poses", "CTF-Tilt", "default-datadir"),
        ],
        indirect=True,
    )
    @pytest.mark.parametrize("indices", [None, "just-5"], indirect=True)
    @pytest.mark.parametrize("ntilts", [None, 5, 20])
    def test_to_fsc(
        self, tmpdir_factory, particles, poses, ctf, indices, datadir, ntilts
    ):
        """Export calculated FSC scores to plotting function."""
        outdir = self.get_outdir(
            tmpdir_factory, particles, poses, ctf, indices, datadir, ntilts
        )
        args = [
            os.path.join(outdir, "vol_fsc-vals.txt"),
            "-o",
            os.path.join(outdir, "vol_fsc-plot2.png"),
        ]

        parser = argparse.ArgumentParser()
        plot_fsc.add_args(parser)
        plot_fsc.main(parser.parse_args(args))
        assert os.path.exists(os.path.join(outdir, "vol_fsc-plot2.png"))

        shutil.rmtree(outdir)
