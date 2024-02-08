"""Unit tests of the cryodrgn clean command."""

import os
import shutil
import subprocess
import pytest
from typing import Union


DATA_DIR = os.path.join(
    os.path.abspath(os.path.dirname(__file__)), "..", "testing", "data"
)


def run_command(cmd: str) -> tuple[str, str]:
    cmd_out = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    return cmd_out.stdout, cmd_out.stderr


class TrainDir:
    def __init__(self, data_lbl: str, train_cmd: str) -> None:
        self.data_lbl = data_lbl
        self.outdir = os.path.abspath(data_lbl)
        self.train_cmd = train_cmd

        if self.data_lbl == "toy":
            particles = os.path.join(DATA_DIR, "toy_projections.mrcs")
            poses = os.path.join(DATA_DIR, "toy_angles.pkl")
        else:
            raise ValueError(f"Unrecognized data label `{self.data_lbl}`!")

        cmd = (
            f"cryodrgn {self.train_cmd} {particles} -o {self.outdir} "
            f"--poses {poses} --no-amp -n=10 "
        )

        if train_cmd == "train_vae":
            cmd += "--zdim=4 --tdim=32 --tlayers=2 "
        elif train_cmd == "train_nn":
            cmd += "--dim=32 --layers=2 "

        out, err = run_command(cmd)
        assert ") Finished in 0" in out
        assert self.all_files_present

    @property
    def out_files(self) -> list[str]:
        return os.listdir(self.outdir)

    def epoch_cleaned(self, epoch: Union[int, None]) -> bool:
        cleaned = True
        out_files = self.out_files

        if epoch:
            epoch_lbl = f".{epoch}"
        else:
            epoch_lbl = ""

        if f"weights{epoch_lbl}.pkl" in out_files:
            cleaned = False

        if self.data_lbl == "train_nn":
            if f"reconstruct{epoch_lbl}.mrc" in out_files:
                cleaned = False
        elif self.data_lbl == "train_vae":
            if f"z{epoch_lbl}.pkl" in out_files:
                cleaned = False

        return cleaned

    @property
    def all_files_present(self) -> bool:
        return not any(self.epoch_cleaned(epoch) for epoch in list(range(10)) + [None])


@pytest.fixture
def train_dir(request) -> TrainDir:
    """Set up and tear down an experiment output directory."""
    data_lbl, train_cmd = request.param
    os.makedirs(data_lbl)
    yield TrainDir(data_lbl, train_cmd)
    shutil.rmtree(data_lbl)


@pytest.fixture
def toy_dir(request) -> TrainDir:
    yield TrainDir("toy", request.param)
    shutil.rmtree("toy")


@pytest.mark.parametrize("toy_dir", ["train_nn", "train_vae"], indirect=True)
@pytest.mark.parametrize("every_n", [1, 2, 3, 5])
def test_clean_here(toy_dir, every_n: int) -> None:
    """Test that we can clean the output of the current directory."""

    os.chdir("toy")
    out, err = run_command(f"cryodrgn_utils clean -n {every_n} -d")
    os.chdir("..")
    assert out == f"\tWould remove {2 * (9 - 9 // every_n)} files!\n"
    assert err == ""
    assert toy_dir.all_files_present

    os.chdir("toy")
    out, err = run_command(f"cryodrgn_utils clean -n {every_n}")
    os.chdir("..")

    assert out == f"\tRemoved {2 * (9 - 9 // every_n)} files!\n"
    assert err == ""

    assert not toy_dir.epoch_cleaned(None)
    for epoch in range(10):
        assert toy_dir.epoch_cleaned(epoch) == (epoch % every_n != 0)


@pytest.mark.parametrize("toy_dir", ["train_nn", "train_vae"], indirect=True)
@pytest.mark.parametrize("every_n", [1, 2, 5])
def test_clean_one(toy_dir, every_n: int) -> None:
    """Test that we can clean the output of one directory."""

    out, err = run_command(f"cryodrgn_utils clean toy -n {every_n} -d")
    assert out == f"\tWould remove {2 * (9 - (9 // every_n))} files!\n"
    assert err == ""
    assert toy_dir.all_files_present

    out, err = run_command(f"cryodrgn_utils clean toy -n {every_n}")
    assert out == f"\tRemoved {2 * (9 - (9 // every_n))} files!\n"
    assert err == ""

    assert not toy_dir.epoch_cleaned(None)
    for epoch in range(10):
        assert toy_dir.epoch_cleaned(epoch) == (epoch % every_n != 0)
