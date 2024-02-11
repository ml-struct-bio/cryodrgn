"""Fixtures used across many unit test modules."""

import os
import shutil
import pytest
from typing import Optional, Union, Any
from cryodrgn.utils import run_command

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


@pytest.fixture(scope="session", autouse=True)
def output_dir() -> None:
    """Helper fixture to remove output folder upon completion of all tests."""
    yield None
    shutil.rmtree("output")


def get_testing_datasets(dataset_lbl: str) -> tuple[str, str]:
    """Retrieve the input files corresponding to a given dataset label."""

    if dataset_lbl == "toy":
        particles = os.path.join(DATA_DIR, "toy_projections.mrcs")
        poses = os.path.join(DATA_DIR, "toy_angles.pkl")
    elif dataset_lbl == "hand":
        particles = os.path.join(DATA_DIR, "hand.mrcs")
        poses = os.path.join(DATA_DIR, "hand_rot.pkl")
    else:
        raise ValueError(f"Unrecognized dataset label `{dataset_lbl}`!")

    return particles, poses


class TrainDir:
    """A folder containing trained cryoDRGN reconstruction output for use in tests.

    Note that the reconstruction model training is done during initialization of an
    instance of this class, and that the `replace_files` method can be used to clean
    up after individual tests via the `orig_cache` subfolder.
    """

    def __init__(
        self,
        dataset: str,
        train_cmd: str,
        epochs: int = 10,
        seed: Optional[int] = None,
        out_lbl: Optional[str] = None,
    ) -> None:
        self.dataset = dataset
        self.particles, self.poses = get_testing_datasets(dataset)

        self.train_cmd = train_cmd
        self.epochs = epochs
        self.out_lbl = out_lbl or "_".join([dataset, train_cmd])
        self.outdir = os.path.abspath(self.out_lbl)
        self.orig_cache = os.path.join(self.outdir, "orig_cache")

        cmd = (
            f"cryodrgn {self.train_cmd} {self.particles} -o {self.outdir} "
            f"--poses {self.poses} --no-amp -n={self.epochs} "
        )

        if self.train_cmd == "train_vae":
            cmd += "--zdim=8 --tdim=16 --tlayers=1 "
        elif self.train_cmd == "train_nn":
            cmd += "--dim=16 --layers=2 "

        if seed:
            cmd += f" --seed={seed}"

        out, err = run_command(cmd)
        assert ") Finished in " in out, err
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

    def train_load_epoch(self, load_epoch: int, train_epochs: int) -> None:
        if not 0 <= load_epoch < self.epochs:
            raise ValueError(
                f"Given epoch to load {load_epoch} is not valid for experiment "
                f"with {self.epochs} epochs!"
            )

        cmd = (
            f"cryodrgn {self.train_cmd} {self.particles} -o {self.outdir} "
            f"--poses {self.poses} --no-amp -n={train_epochs} "
            f"--load {os.path.join(self.outdir, f'weights.{load_epoch}.pkl')} "
        )
        if self.train_cmd == "train_vae":
            cmd += "--zdim=8 --tdim=16 --tlayers=1 "
        elif self.train_cmd == "train_nn":
            cmd += "--dim=16 --layers=2 "

        out, err = run_command(cmd)
        assert ") Finished in " in out, err
        self.epochs = train_epochs
        assert self.all_files_present


# Note the use of nested fixtures here with a function->session scope hierarchy that
# allows for separate setup/teardown routines to be defined for the same `train_dir`
# fixture across these scopes. Thus each time we need a directory with trained
# experiment results we get one and replace any moved/deleted output files afterwards,
# but the experiment is only run as often as necessary for each session, with the output
# folder being cleaned up each time.


@pytest.fixture(scope="session")
def train_dir(request) -> TrainDir:
    """Run an experiment to generate output; remove this output when finished."""
    tdir = TrainDir(**TrainDir.parse_request(request.param))
    yield tdir
    shutil.rmtree(tdir.outdir)


@pytest.fixture(scope="function")
def trained_dir(train_dir) -> TrainDir:
    """Get an experiment that has been run; restore its output when done."""
    yield train_dir
    train_dir.replace_files()


@pytest.fixture(scope="session")
def train_dirs(request) -> list[TrainDir]:
    """Run experiments to generate outputs; remove these outputs when finished."""
    tdirs = [TrainDir(**TrainDir.parse_request(req)) for req in request.param]
    yield tdirs
    for tdir in tdirs:
        shutil.rmtree(tdir.outdir)


@pytest.fixture(scope="function")
def trained_dirs(train_dirs) -> list[TrainDir]:
    """Get experiments that have been run; restore their output when done."""
    yield train_dirs
    for train_dir in train_dirs:
        train_dir.replace_files()


class AbInitioDir:
    """A folder for running cryoDRGN ab initio experiment tests."""

    def __init__(
        self,
        zdim: int,
        dataset: str = "hand",
        epochs: int = 2,
        seed: Optional[int] = None,
    ) -> None:
        self.zdim = zdim
        self.dataset = dataset
        self.particles, _ = get_testing_datasets(dataset)
        self.outdir = os.path.abspath(f"test-output_{dataset}")
        os.makedirs(self.outdir)
        self.epochs = epochs
        self.seed = seed

    @classmethod
    def parse_request(cls, req: dict[str, Any]) -> dict[str, Any]:
        train_args = dict()

        if "zdim" not in req:
            raise ValueError("AbinitioDir fixture request must specify a zdim!")
        train_args["zdim"] = req["zdim"]

        if "dataset" in req:
            train_args["dataset"] = req["dataset"]
        else:
            train_args["dataset"] = "hand"

        if "epochs" in req:
            train_args["epochs"] = req["epochs"]
        else:
            train_args["epochs"] = 2

        if "seed" in req:
            train_args["seed"] = req["seed"]
        else:
            train_args["seed"] = None

        return train_args

    def train(self, load_epoch: Optional[int] = None) -> None:
        train_cmd = "abinit_het" if self.zdim > 0 else "abinit_homo"

        cmd = (
            f"cryodrgn {train_cmd} {self.particles} -o {self.outdir} "
            f"--num-epochs {self.epochs} --no-window --pretrain 100 "
        )
        if self.zdim > 0:
            cmd += f"--zdim {self.zdim} "
            cmd += "--enc-dim 8 --enc-layers 2 --dec-dim 8 --pe-dim 8 "
        else:
            cmd += "--dim 16 --pe-dim 8 "

        if self.seed is not None:
            cmd += f"--seed={self.seed} "

        if load_epoch is not None:
            cmd += f"--load {os.path.join(self.outdir, f'weights.{load_epoch}.pkl')} "
            cmd += (
                f"--load-poses {os.path.join(self.outdir, f'pose.{load_epoch}.pkl')} "
            )

        out, err = run_command(cmd)
        assert ") Finished in " in out, err

    def analyze(self, analysis_epoch: int) -> None:
        run_command(f"cryodrgn analyze {self.outdir} {analysis_epoch}")

    def backproject(self) -> None:
        run_command(
            f"cryodrgn backproject_voxel {self.particles} "
            f"-o {os.path.join(self.outdir, 'backproject', 'vol.mrc')} "
            f"--poses {os.path.join(self.outdir, 'pose.pkl')} "
        )

    def view_config(self) -> None:
        run_command(f"cryodrgn view_config {self.outdir}")


@pytest.fixture(scope="session")
def abinit_dir(request) -> AbInitioDir:
    adir = AbInitioDir(**AbInitioDir.parse_request(request.param))
    yield adir
    shutil.rmtree(adir.outdir)
