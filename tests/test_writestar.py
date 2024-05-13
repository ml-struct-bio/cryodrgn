import pytest
import os
import argparse
from cryodrgn.commands_utils import write_star


@pytest.fixture
def particles_starfile():
    return os.path.join(pytest.data_dir, "FinalRefinement-OriginalParticles-PfCRT.star")


@pytest.mark.parametrize(
    "particles, ctf, poses",
    [
        ("hand", "CTF-Test.100", {"rot": "hand-rot", "pose": "hand-poses"}),
        ("toy.mrcs", "CTF-Test", {"pose": "toy-poses", "rot": "toy-angles"}),
    ],
    indirect=True,
    ids=("hand", "toy"),
)
@pytest.mark.parametrize(
    "use_ctf",
    [
        pytest.param(False, marks=pytest.mark.xfail(reason="ctf required for .mrcs")),
        True,
    ],
    ids=("no-ctf", "with-ctf"),
)
@pytest.mark.parametrize(
    "use_poses",
    [
        None,
        pytest.param(
            "rot", marks=pytest.mark.xfail(reason="old pose format not supported")
        ),
        "pose",
    ],
    ids=("no-poses", "with-rot", "with-poses"),
)
@pytest.mark.parametrize(
    "indices", [None, "just-5"], indirect=True, ids=("no-ind", "with-ind")
)
@pytest.mark.parametrize("use_relion30", [False, True], ids=("relion3.1", "relion3.0"))
def test_from_mrcs(
    outdir, particles, ctf, poses, use_ctf, use_poses, indices, use_relion30
):
    parser = argparse.ArgumentParser()
    write_star.add_args(parser)
    args = [particles, "-o", os.path.join(outdir, "out.star")]

    if use_ctf:
        args += ["--ctf", ctf]
    if use_poses:
        args += ["--poses", poses[use_poses]]
    if indices:
        args += ["--ind", indices]
    if use_relion30:
        args += ["--relion30"]

    write_star.main(parser.parse_args(args))
    assert os.path.exists(os.path.join(outdir, "out.star"))
