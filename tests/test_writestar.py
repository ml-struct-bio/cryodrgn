import pytest
import os
import argparse
from cryodrgn.commands_utils import write_star, filter_mrcs
from cryodrgn.starfile import Starfile
from cryodrgn.source import ImageSource
import cryodrgn.utils


@pytest.fixture
def particles_starfile():
    return os.path.join(pytest.data_dir, "FinalRefinement-OriginalParticles-PfCRT.star")


@pytest.mark.parametrize(
    "use_ctf",
    [
        pytest.param(False, marks=pytest.mark.xfail(reason="ctf is required input")),
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
        self, tmpdir, particles, ctf, poses, use_ctf, use_poses, indices, use_relion30
    ):
        parser = argparse.ArgumentParser()
        write_star.add_args(parser)
        args = [particles, "-o", os.path.join(tmpdir, "out.star")]

        if use_ctf:
            args += ["--ctf", ctf]
        if use_poses:
            args += ["--poses", poses[use_poses]]
        if indices:
            args += ["--ind", indices]
        if use_relion30:
            args += ["--relion30"]

        write_star.main(parser.parse_args(args))
        assert os.path.exists(os.path.join(tmpdir, "out.star"))

        ind = cryodrgn.utils.load_pkl(indices) if indices else None
        particle_data = ImageSource.from_file(particles, indices=ind)
        ind = list(range(particle_data.n)) if ind is None else ind
        star_data = Starfile.load(os.path.join(tmpdir, "out.star"))
        assert star_data.df.shape[0] == particle_data.n

        if use_relion30:
            assert star_data.data_optics is None
        else:
            assert isinstance(star_data.data_optics, Starfile)
            assert star_data.data_optics.df.shape == (1, 6)

        for i, star_name in enumerate(star_data.df["_rlnImageName"].tolist()):
            file_idx, file_name = star_name.split("@")
            assert int(file_idx) == ind[i] + 1
            assert os.path.basename(file_name) == os.path.basename(particles)

    @pytest.mark.parametrize(
        "particles, ctf, poses",
        [
            ("toy.txt", "CTF-Test", {"pose": "toy-poses", "rot": "toy-angles"}),
        ],
        indirect=True,
        ids=("hand",),
    )
    def test_from_txt(
        self, tmpdir, particles, ctf, poses, use_ctf, use_poses, indices, use_relion30
    ):
        parser = argparse.ArgumentParser()
        write_star.add_args(parser)
        args = [particles, "-o", os.path.join(tmpdir, "out.star")]

        if use_ctf:
            args += ["--ctf", ctf]
        if use_poses:
            args += ["--poses", poses[use_poses]]
        if indices:
            args += ["--ind", indices]
        if use_relion30:
            args += ["--relion30"]

        write_star.main(parser.parse_args(args))
        assert os.path.exists(os.path.join(tmpdir, "out.star"))

        ind = cryodrgn.utils.load_pkl(indices) if indices else None
        particle_data = ImageSource.from_file(particles, indices=ind)
        ind = list(range((particle_data.n))) if ind is None else ind
        star_data = Starfile.load(os.path.join(tmpdir, "out.star"))
        assert star_data.df.shape[0] == particle_data.n

        if use_relion30:
            assert star_data.data_optics is None
        else:
            assert isinstance(star_data.data_optics, Starfile)
            assert star_data.data_optics.df.shape == (1, 6)

        for i, star_name in enumerate(star_data.df["_rlnImageName"].tolist()):
            file_idx, file_name = star_name.split("@")
            assert int(file_idx) == ind[i] + 1
            assert file_name == os.path.basename(
                particle_data.df.loc[ind[i], "__mrc_filename"]
            )


@pytest.mark.parametrize(
    "particles, ctf, poses, indices",
    [("toy.mrcs", "CTF-Test", "toy-poses", "random-100")],
    indirect=True,
    ids=("toy",),
)
@pytest.mark.parametrize("use_relion30", [False, True], ids=("relion3.1", "relion3.0"))
@pytest.mark.parametrize("use_indices", [False, True], ids=("no-ind", "with-ind"))
def test_from_txt_with_two_files(
    tmpdir, particles, ctf, poses, indices, use_relion30, use_indices
):
    orig_dir = os.getcwd()
    os.chdir(tmpdir)
    mrcs_files = [f"particles{i}.mrcs" for i in range(10)]

    for mrcs_file in mrcs_files:
        args = filter_mrcs.add_args(argparse.ArgumentParser()).parse_args(
            [particles, "--ind", indices, "-o", mrcs_file]
        )
        filter_mrcs.main(args)

    txt_file = os.path.join(tmpdir, "particles.txt")
    with open(txt_file, "w") as f:
        f.write("\n".join(mrcs_files))

    parser = argparse.ArgumentParser()
    write_star.add_args(parser)
    args = [
        txt_file,
        "-o",
        os.path.join(tmpdir, "out.star"),
        "--ctf",
        ctf,
        "--poses",
        poses,
    ]
    if use_relion30:
        args += ["--relion30"]
    if use_indices:
        args += ["--ind", indices]
    write_star.main(parser.parse_args(args))
    assert os.path.exists(os.path.join(tmpdir, "out.star"))

    ind = cryodrgn.utils.load_pkl(indices) if use_indices else None
    particle_data = ImageSource.from_file(txt_file, indices=ind)
    ind = list(range(particle_data.n)) if ind is None else ind
    star_data = Starfile.load(os.path.join(tmpdir, "out.star"))
    assert star_data.df.shape[0] == particle_data.n

    if use_relion30:
        assert star_data.data_optics is None
    else:
        assert isinstance(star_data.data_optics, Starfile)
        assert star_data.data_optics.df.shape == (1, 6)

    for i, star_name in enumerate(star_data.df["_rlnImageName"].tolist()):
        file_idx, file_name = star_name.split("@")
        assert int(file_idx) == (ind[i] % 100) + 1
        assert file_name == os.path.basename(
            particle_data.df.loc[ind[i], "__mrc_filename"]
        )

    os.chdir(orig_dir)
