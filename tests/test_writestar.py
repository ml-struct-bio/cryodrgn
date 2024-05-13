import pytest
import os
import argparse
from cryodrgn.commands_utils import write_star
from cryodrgn.starfile import Starfile
from cryodrgn.source import ImageSource
import cryodrgn.utils


@pytest.fixture
def particles_starfile():
    return os.path.join(pytest.data_dir, "FinalRefinement-OriginalParticles-PfCRT.star")


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
class TestBasic:
    @pytest.mark.parametrize(
        "particles, ctf, poses",
        [
            ("hand", "CTF-Test.100", {"rot": "hand-rot", "pose": "hand-poses"}),
            ("toy.mrcs", "CTF-Test", {"pose": "toy-poses", "rot": "toy-angles"}),
        ],
        indirect=True,
        ids=("hand", "toy"),
    )
    def test_from_mrcs(
        self, outdir, particles, ctf, poses, use_ctf, use_poses, indices, use_relion30
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

        ind = cryodrgn.utils.load_pkl(indices) if indices else None
        particle_data = ImageSource.from_file(particles, indices=ind)
        star_data = Starfile.load(os.path.join(outdir, "out.star"))
        assert star_data.df.shape[0] == particle_data.n

        if use_relion30:
            assert star_data.data_optics is None
        else:
            assert isinstance(star_data.data_optics, Starfile)
            assert star_data.data_optics.df.shape == (1, 6)

        file_names = {
            os.path.basename(name.split("@")[1])
            for name in star_data.df["_rlnImageName"]
        }
        assert file_names == {os.path.basename(particles)}

    @pytest.mark.parametrize(
        "particles, ctf, poses",
        [
            ("toy.txt", "CTF-Test", {"pose": "toy-poses", "rot": "toy-angles"}),
        ],
        indirect=True,
        ids=("hand",),
    )
    def test_from_txt(
        self, outdir, particles, ctf, poses, use_ctf, use_poses, indices, use_relion30
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

        ind = cryodrgn.utils.load_pkl(indices) if indices else None
        particle_data = ImageSource.from_file(particles, indices=ind)
        star_data = Starfile.load(os.path.join(outdir, "out.star"))
        assert star_data.df.shape[0] == particle_data.n

        if use_relion30:
            assert star_data.data_optics is None
        else:
            assert isinstance(star_data.data_optics, Starfile)
            assert star_data.data_optics.df.shape == (1, 6)
