"""Tests of compatibility with the RELION 3.1 format for .star files (optics groups)."""

import pytest
import argparse
import os
import pickle
import numpy as np
from cryodrgn.commands import parse_ctf_star, parse_pose_star
from cryodrgn.commands_utils import filter_star, select_random
from cryodrgn.starfile import parse_star, write_star, Starfile
from cryodrgn.utils import load_pkl


@pytest.fixture
def relion_starfile(request):
    return os.path.join(pytest.DATADIR, request.param)


pytestmark = pytest.mark.parametrize(
    "relion_starfile",
    ["relion31.star", "relion31.v2.star", "relion31.6opticsgroups.star"],
    indirect=True,
)


@pytest.mark.parametrize("index_fraction, index_seed", [(0.4, 55), (0.3, 101)])
class TestFilterStar:
    def get_outdir(self, tmpdir_factory, relion_starfile, index_seed, index_fraction):
        dirname = os.path.join(
            "r31_FilterStar", relion_starfile, str(index_seed), str(index_fraction)
        )
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)

        return odir

    def test_command(self, tmpdir_factory, relion_starfile, index_seed, index_fraction):
        outlbl = os.path.basename(relion_starfile)
        outdir = self.get_outdir(tmpdir_factory, outlbl, index_seed, index_fraction)
        indata, in_optics = parse_star(relion_starfile)
        sel_file = os.path.join(outdir, "random-index.pkl")

        parser = argparse.ArgumentParser()
        select_random.add_args(parser)
        select_random.main(
            parser.parse_args(
                [
                    str(indata.shape[0]),
                    "-o",
                    sel_file,
                    "--frac",
                    str(index_fraction),
                    "--seed",
                    str(index_seed),
                ]
            )
        )
        selected = load_pkl(sel_file)
        assert len(selected) == int(indata.shape[0] * index_fraction)

        parser = argparse.ArgumentParser()
        filter_star.add_args(parser)
        outfile = os.path.join(
            outdir, f"fltr_{outlbl}_{index_seed}-{index_fraction}.star"
        )
        args = [f"{relion_starfile}", "-o", outfile, "--ind", sel_file]
        filter_star.main(parser.parse_args(args))

        outdata, out_optics = parse_star(outfile)
        assert (out_optics is None) == (in_optics is None)
        assert outdata.shape[0] == len(selected)
        assert out_optics.shape[0] == len(indata["_rlnOpticsGroup"][selected].unique())
        assert (indata.loc[selected].values == outdata.values).all()

    def test_relion30_consistency(
        self, tmpdir_factory, relion_starfile, index_seed, index_fraction
    ):
        outlbl = os.path.basename(relion_starfile)
        outdir = self.get_outdir(tmpdir_factory, outlbl, index_seed, index_fraction)

        starfile = Starfile(relion_starfile)
        sel_file = os.path.join(outdir, "random-index.pkl")
        write_star(os.path.join(outdir, "r30.star"), data=starfile.to_relion30())
        parser = argparse.ArgumentParser()
        filter_star.add_args(parser)
        new_outfile = os.path.join(
            outdir, f"fltr_{outlbl}_{index_seed}-{index_fraction}_r30.star"
        )
        args = [os.path.join(outdir, "r30.star"), "-o", new_outfile, "--ind", sel_file]
        filter_star.main(parser.parse_args(args))

        orig_starfile = Starfile(
            os.path.join(outdir, f"fltr_{outlbl}_{index_seed}-{index_fraction}.star")
        )
        selected = load_pkl(sel_file)
        new_starfile = Starfile(new_outfile)
        assert not new_starfile.relion31
        assert new_starfile.df.shape[0] == len(selected)
        new_data = new_starfile.df.loc[:, orig_starfile.df.columns].values
        assert (orig_starfile.df.values == new_data).all()


@pytest.mark.parametrize(
    "apix, resolution", [(None, None), (1.5, None), (3.0, 256), (None, 128), (2.0, 128)]
)
class TestParsePoseStar:
    def get_outdir(self, tmpdir_factory, relion_starfile, apix, resolution):
        dirname = os.path.join(
            "r31_ParsePoseStar", relion_starfile, str(apix), str(resolution)
        )
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)

        return odir

    def test_command(self, tmpdir_factory, relion_starfile, apix, resolution):
        outlbl = os.path.basename(relion_starfile)
        outdir = self.get_outdir(tmpdir_factory, outlbl, apix, resolution)

        parser = argparse.ArgumentParser()
        parse_pose_star.add_args(parser)
        pose_file = os.path.join(outdir, "orig-poses.pkl")
        starfile = Starfile(relion_starfile)
        orig_apix, orig_D = starfile.apix, starfile.resolution

        args = [f"{relion_starfile}", "-o", pose_file]
        parse_pose_star.main(parser.parse_args(args))
        with open(pose_file, "rb") as f:
            rots, trans = pickle.load(f)

        assert rots.shape == (starfile.df.shape[0], 3, 3)
        assert trans.shape == (starfile.df.shape[0], 2)

        new_posefile = os.path.join(outdir, "parsed-poses.pkl")
        args = [f"{relion_starfile}", "-o", new_posefile]
        if apix is not None:
            args += ["--Apix", str(apix)]
        if resolution is not None:
            args += ["-D", str(resolution)]

        parse_pose_star.main(parser.parse_args(args))
        with open(new_posefile, "rb") as f:
            new_rots, new_trans = pickle.load(f)

        # only translation get modified when we don't use _rlnOrigin[X/Y]Angst
        assert new_rots.shape == (starfile.df.shape[0], 3, 3)
        assert np.allclose(rots, new_rots)
        assert new_trans.shape == (starfile.df.shape[0], 2)

        check_trans = trans.copy()
        if apix is not None:
            check_trans = (check_trans.T * orig_apix / apix).T
        if resolution is not None:
            check_trans = (check_trans.T * orig_D / resolution).T

        assert np.allclose(check_trans, new_trans)

    def test_relion30_consistency(
        self, tmpdir_factory, relion_starfile, apix, resolution
    ):
        outlbl = os.path.basename(relion_starfile)
        outdir = self.get_outdir(tmpdir_factory, outlbl, apix, resolution)

        starfile = Starfile(relion_starfile)
        write_star(os.path.join(outdir, "r30.star"), data=starfile.to_relion30())
        pose_file = os.path.join(outdir, "parsed-poses_r30.pkl")
        parser = argparse.ArgumentParser()
        parse_pose_star.add_args(parser)
        args = [os.path.join(outdir, "r30.star"), "-o", pose_file]
        if apix is not None:
            args += ["--Apix", str(apix)]
        if resolution is not None:
            args += ["-D", str(resolution)]

        parse_pose_star.main(parser.parse_args(args))
        with open(pose_file, "rb") as f:
            rots, trans = pickle.load(f)

        assert rots.shape == (starfile.df.shape[0], 3, 3)
        assert trans.shape == (starfile.df.shape[0], 2)

        old_posefile = os.path.join(outdir, "parsed-poses.pkl")
        with open(old_posefile, "rb") as f:
            old_rots, old_trans = pickle.load(f)

        assert np.allclose(old_rots, rots)
        assert np.allclose(old_trans, trans)


@pytest.mark.parametrize(
    "apix, resolution", [(None, None), (1.5, None), (3.0, 256), (None, 128), (2.0, 128)]
)
@pytest.mark.parametrize("kv", [None, 300])
@pytest.mark.parametrize("cs", [None, 2.7])
@pytest.mark.parametrize("w", [None, 0.15])
@pytest.mark.parametrize("ps", [None, 1.0])
class TestParseCTFStar:
    def get_outdir(
        self, tmpdir_factory, relion_starfile, apix, resolution, kv, cs, w, ps
    ):
        dirname = os.path.join(
            "r31_ParseCTFStar",
            relion_starfile,
            str(apix),
            str(resolution),
            str(kv),
            str(cs),
            str(w),
            str(ps),
        )
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)

        return odir

    def test_command(
        self, tmpdir_factory, relion_starfile, apix, resolution, kv, cs, w, ps
    ):
        outlbl = os.path.basename(relion_starfile)
        outdir = self.get_outdir(
            tmpdir_factory, outlbl, apix, resolution, kv, cs, w, ps
        )

        parser = argparse.ArgumentParser()
        parse_ctf_star.add_args(parser)
        out_fl = os.path.join(outdir, "parsed-ctf.pkl")
        args = [relion_starfile, "-o", out_fl]
        if apix is not None:
            args += ["--Apix", str(apix)]
        if resolution is not None:
            args += ["-D", str(resolution)]
        if kv is not None:
            args += ["--kv", str(kv)]
        if cs is not None:
            args += ["--cs", str(cs)]
        if w is not None:
            args += ["-w", str(w)]
        if ps is not None:
            args += ["--ps", str(ps)]

        parse_ctf_star.main(parser.parse_args(args))
        starfile = Starfile(relion_starfile)
        orig_apix, orig_D = starfile.apix, starfile.resolution
        with open(os.path.join(out_fl), "rb") as f:
            ctf_params = pickle.load(f)

        new_apix = apix or orig_apix
        new_D = resolution or orig_D
        new_kv = kv or starfile.get_optics_values("_rlnVoltage", dtype=float)
        new_cs = cs or starfile.get_optics_values(
            "_rlnSphericalAberration", dtype=float
        )
        new_w = w or starfile.get_optics_values("_rlnAmplitudeContrast", dtype=float)
        new_ps = ps or starfile.get_optics_values("_rlnPhaseShift", dtype=float)

        assert ctf_params.shape == (starfile.df.shape[0], 9)
        assert (ctf_params[:, 1] == new_apix).all()
        assert (ctf_params[:, 0] == new_D).all()

        assert np.allclose(
            starfile.get_optics_values("_rlnDefocusU", dtype=float), ctf_params[:, 2]
        )
        assert np.allclose(
            starfile.get_optics_values("_rlnDefocusV", dtype=float), ctf_params[:, 3]
        )
        assert np.allclose(new_kv, ctf_params[:, 5])
        assert np.allclose(new_cs, ctf_params[:, 6])
        assert np.allclose(new_w, ctf_params[:, 7])
        assert np.allclose(new_ps, ctf_params[:, 8])

    def test_relion30_consistency(
        self, tmpdir_factory, relion_starfile, apix, resolution, kv, cs, w, ps
    ):
        outlbl = os.path.basename(relion_starfile)
        outdir = self.get_outdir(
            tmpdir_factory, outlbl, apix, resolution, kv, cs, w, ps
        )

        starfile = Starfile(relion_starfile)
        write_star(os.path.join(outdir, "r30.star"), data=starfile.to_relion30())
        out_fl = os.path.join(outdir, "parsed-ctf_r30.pkl")
        parser = argparse.ArgumentParser()
        parse_ctf_star.add_args(parser)
        args = [relion_starfile, "-o", out_fl]
        if apix is not None:
            args += ["--Apix", str(apix)]
        if resolution is not None:
            args += ["-D", str(resolution)]
        if kv is not None:
            args += ["--kv", str(kv)]
        if cs is not None:
            args += ["--cs", str(cs)]
        if w is not None:
            args += ["-w", str(w)]
        if ps is not None:
            args += ["--ps", str(ps)]

        parse_ctf_star.main(parser.parse_args(args))
        with open(os.path.join(out_fl), "rb") as f:
            ctf_params = pickle.load(f)
        with open(os.path.join(outdir, "parsed-ctf.pkl"), "rb") as f:
            orig_params = pickle.load(f)

        assert np.allclose(ctf_params, orig_params)


def test_relion50(tmpdir, relion_starfile):
    with open(relion_starfile, "r") as f:
        starlines = f.readlines()

    starlines += [
        "\n",
        "# version 50001\n",
        "\n",
        "data_general\n",
        "\n",
        "_rlnTomoSubTomosAre2DStacks 1\n",
    ]
    with open(os.path.join(tmpdir, "new.star"), "w") as f:
        f.writelines(starlines)

    newfile = Starfile(os.path.join(tmpdir, "new.star"))
    starfile = Starfile(relion_starfile)
    assert newfile == starfile
