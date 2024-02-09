"""Fixtures used across many unit test modules."""

import os
import shutil
import pytest
from typing import Optional, Union, Any
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
        out_lbl: Optional[str] = None,
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

    @classmethod
    def parse_request(cls, req: dict[str, Any]) -> dict[str, Any]:
        train_args = dict()

        if "dataset" in req:
            train_args["dataset"] = req["dataset"]
        else:
            train_args["dataset"] = "hand"

        if "train_cmd" in req:
            train_args["train_cmd"] = req["train_cmd"]
        else:
            train_args["train_cmd"] = "train_nn"

        if "epochs" in req:
            train_args["epochs"] = req["epochs"]
        else:
            train_args["epochs"] = 10

        if "seed" in req:
            train_args["seed"] = req["seed"]
        else:
            train_args["seed"] = None

        train_args["out_lbl"] = "_".join([str(x) for x in train_args.values()])

        return train_args

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


@pytest.fixture(scope="session")
def train_dir(request) -> TrainDir:
    """Run an experiment to generate output; remove this output when finished."""
    train_args = TrainDir.parse_request(request.param)
    yield TrainDir(**train_args)
    shutil.rmtree(train_args["out_lbl"])


@pytest.fixture(scope="function")
def trained_dir(train_dir) -> TrainDir:
    """Get an experiment that has been run; restore its output when done."""
    yield train_dir
    train_dir.replace_files()


@pytest.fixture(scope="session")
def train_dirs(request) -> list[TrainDir]:
    """Run experiments to generate outputs; remove these outputs when finished."""
    train_args_list = [TrainDir.parse_request(req) for req in request.param]
    yield [TrainDir(**train_args) for train_args in train_args_list]
    for train_args in train_args_list:
        shutil.rmtree(train_args["out_lbl"])


@pytest.fixture(scope="function")
def trained_dirs(train_dirs) -> list[TrainDir]:
    """Get experiments that have been run; restore their output when done."""
    yield train_dirs
    for train_dir in train_dirs:
        train_dir.replace_files()
