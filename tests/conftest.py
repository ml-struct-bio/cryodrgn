"""Fixtures used across many unit test modules."""

import os
import shutil
import pytest
from typing import Optional, Union
from cryodrgn.utils import run_command


DATA_DIR = os.path.join(
    os.path.abspath(os.path.dirname(__file__)), "..", "testing", "data"
)


class TrainDir:
    def __init__(
        self,
        dataset: str,
        train_cmd: str,
        epochs: int = 10,
        seed: Optional[int] = None,
        out_lbl: Optional[str] = None
    ) -> None:
        self.dataset = dataset
        self.train_cmd = train_cmd
        self.epochs = epochs
        self.out_lbl = out_lbl or "_".join([dataset, train_cmd])
        self.outdir = os.path.abspath(self.out_lbl)
        self.orig_cache = os.path.join(self.outdir, "orig_cache")

        if self.dataset == "toy":
            particles = os.path.join(DATA_DIR, "toy_projections.mrcs")
            poses = os.path.join(DATA_DIR, "toy_angles.pkl")
        elif self.dataset == "hand":
            particles = os.path.join(DATA_DIR, "hand.mrcs")
            poses = os.path.join(DATA_DIR, "hand_rot.pkl")
        else:
            raise ValueError(f"Unrecognized dataset label `{self.dataset}`!")

        cmd = (
            f"cryodrgn {self.train_cmd} {particles} -o {self.outdir} "
            f"--poses {poses} --no-amp -n={self.epochs} "
        )

        if train_cmd == "train_vae":
            cmd += "--zdim=8 --tdim=16 --tlayers=1 "
        elif train_cmd == "train_nn":
            cmd += "--dim=16 --layers=2 "

        if seed:
            cmd += f" --seed={seed}"

        out, err = run_command(cmd)
        assert ") Finished in " in out
        assert self.all_files_present

        orig_files = self.out_files
        os.makedirs(self.orig_cache)
        for orig_file in orig_files:
            shutil.copy(
                os.path.join(self.outdir, orig_file),
                os.path.join(self.orig_cache, orig_file),
            )

    @property
    def out_files(self) -> list[str]:
        return os.listdir(self.outdir)

    def epoch_cleaned(self, epoch: Union[int, None]) -> bool:
        if epoch and not 0 <= epoch < self.epochs:
            raise ValueError(
                f"Cannot check if given epoch {epoch} has been cleaned "
                f"for output folder `{self.outdir}` which only contains "
                f"{self.epochs} epochs!"
            )

        cleaned = True
        out_files = self.out_files

        if epoch:
            epoch_lbl = f".{epoch}"
        else:
            epoch_lbl = ""

        if f"weights{epoch_lbl}.pkl" in out_files:
            cleaned = False

        if self.train_cmd == "train_nn":
            if f"reconstruct{epoch_lbl}.mrc" in out_files:
                cleaned = False
        elif self.train_cmd == "train_vae":
            if f"z{epoch_lbl}.pkl" in out_files:
                cleaned = False

        return cleaned

    @property
    def all_files_present(self) -> bool:
        return not any(
            self.epoch_cleaned(epoch) for epoch in list(range(self.epochs)) + [None]
        )

    def replace_files(self) -> None:
        for orig_file in os.listdir(self.orig_cache):
            shutil.copy(
                os.path.join(self.orig_cache, orig_file),
                os.path.join(self.outdir, orig_file),
            )


@pytest.fixture(scope="module")
def train_dir(request) -> TrainDir:
    """Set up and tear down an experiment output directory."""
    dataset = request.param["dataset"] if "dataset" in request.param else "hand"
    cmd = request.param["train_cmd"] if "train_cmd" in request.param else "train_nn"
    epochs = request.param["epochs"] if "epochs" in request.param else 10
    seed = request.param["seed"] if "seed" in request.param else None

    out_lbl = f"{request.module.__name__}__{dataset}_{cmd}"
    yield TrainDir(dataset, cmd, epochs, seed, out_lbl)
    shutil.rmtree(out_lbl)


@pytest.fixture(scope="module")
def train_dirs(request) -> list[TrainDir]:
    """Set up and tear down a series of experiment output directories."""
    out_lbls = {
        (dataset, train_cmd): f"{request.module.__name__}-dirs__{dataset}_{train_cmd}"
        for dataset, train_cmd in request.param
    }

    yield [
        TrainDir(dataset, train_cmd, out_lbl=out_lbl)
        for (dataset, train_cmd), out_lbl in out_lbls.items()
    ]

    for out_lbl in out_lbls.values():
        shutil.rmtree(out_lbl)
