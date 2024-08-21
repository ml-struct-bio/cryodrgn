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


# TODO -- convert these starfiles to 3.0 and compare outputs across these tests
@pytest.fixture
def relion_starfile(request):
    return os.path.join(pytest.DATADIR, request.param)


pytestmark = pytest.mark.parametrize(
    "relion_starfile",
    ["relion31.star", "relion31.v2.star", "relion31.6opticsgroups.star"],
    indirect=True,
)


@pytest.mark.parametrize("index_fraction, index_seed", [(0.4, 55), (0.3, 101)])
def test_filter_star(tmpdir, relion_starfile, index_seed, index_fraction):
    indata, in_optics = parse_star(relion_starfile)
    sel_file = os.path.join(tmpdir, "random-index.pkl")

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
    outlbl = os.path.basename(relion_starfile)
    outfile = os.path.join(tmpdir, f"fltr_{outlbl}_{index_seed}-{index_fraction}.star")
    args = [f"{relion_starfile}", "-o", outfile, "--ind", sel_file]
    filter_star.main(parser.parse_args(args))

    outdata, out_optics = parse_star(outfile)
    assert (out_optics is None) == (in_optics is None)
    assert outdata.shape[0] == len(selected)
    assert out_optics.shape[0] == len(indata["_rlnOpticsGroup"][selected].unique())
    assert (indata.loc[selected].values == outdata.values).all()


@pytest.mark.parametrize(
    "apix, resolution", [(None, None), (1.5, None), (3.0, 256), (None, 128), (2.0, 128)]
)
def test_parse_pose_star(tmpdir, relion_starfile, apix, resolution):
    parser = argparse.ArgumentParser()
    parse_pose_star.add_args(parser)
    pose_file = os.path.join(tmpdir, "orig-poses.pkl")
    starfile = Starfile(relion_starfile)
    orig_apix, orig_D = starfile.apix, starfile.resolution

    args = [f"{relion_starfile}", "-o", pose_file]
    parse_pose_star.main(parser.parse_args(args))
    with open(pose_file, "rb") as f:
        rots, trans = pickle.load(f)

    assert rots.shape == (starfile.df.shape[0], 3, 3)
    assert trans.shape == (starfile.df.shape[0], 2)

    new_posefile = os.path.join(tmpdir, "parsed-poses.pkl")
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


@pytest.mark.parametrize(
    "apix, resolution", [(None, None), (1.5, None), (3.0, 256), (None, 128), (2.0, 128)]
)
def test_parse_ctf_star(tmpdir, relion_starfile, apix, resolution):
    parser = argparse.ArgumentParser()
    parse_ctf_star.add_args(parser)
    out_fl = os.path.join(tmpdir, "parsed-ctf.pkl")
    args = [
        relion_starfile,
        "--kv",
        "300",
        "-w",
        ".1",
        "--ps",
        "0",
        "--cs",
        "2.7",
        "-o",
        out_fl,
    ]
    if apix is not None:
        args += ["--Apix", str(apix)]
    if resolution is not None:
        args += ["-D", str(resolution)]

    parse_ctf_star.main(parser.parse_args(args))
    starfile = Starfile(relion_starfile)
    orig_apix, orig_D = starfile.apix, starfile.resolution
    with open(os.path.join(out_fl), "rb") as f:
        ctf_params = pickle.load(f)

    new_apix = apix or orig_apix
    new_D = resolution or orig_D
    assert ctf_params.shape == (starfile.df.shape[0], 9)
    assert (ctf_params[:, 1] == new_apix).all()
    assert (ctf_params[:, 0] == new_D).all()


def test_to_relion30(tmpdir, relion_starfile):
    starfile = Starfile(relion_starfile)
    write_star(os.path.join(tmpdir, "r30.star"), data=starfile.to_relion30())
    r30_starfile = Starfile.load(os.path.join(tmpdir, "r30.star"))
    assert r30_starfile.df.shape[0] == starfile.df.shape[0]
    assert not r30_starfile.relion31
