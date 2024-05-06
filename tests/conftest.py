"""Fixtures used across many unit test modules."""

import pytest
import os
import shutil
from typing import Optional, Union, Any
from cryodrgn.utils import run_command

os.environ["NUMEXPR_NUM_THREADS"] = "1"
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


def pytest_configure():
    pytest.data_dir = DATA_DIR


@pytest.fixture(scope="session", autouse=True)
def default_outdir() -> None:
    """Helper fixture to remove default  output folder upon completion of all tests."""
    yield None

    # we don't always create this folder, e.g. if we are only doing some of the tests
    if os.path.exists("output"):
        shutil.rmtree("output")


@pytest.fixture(scope="class")
def outdir(tmpdir_factory, request) -> str:
    odir = tmpdir_factory.mktemp(f"output_{request.node.__class__.__name__}")
    yield str(odir)
    shutil.rmtree(odir)


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


PARTICLES_FILES = {
    "hand": "hand.mrcs",
    "hand-tilt": "hand_tilt.mrcs",
    "toy.mrcs": "toy_projections.mrcs",
    "toy.mrcs-999": "toy_projections_0-999.mrcs",
    "toy.star": "toy_projections.star",
    "toy.star-13": "toy_projections_13.star",
    "toy.txt": "toy_projections.txt",
    "tilts.star": "sta_testing_bin8.star",
    "cryosparc-all": "cryosparc_P12_J24_001_particles.cs",
}
POSES_FILES = {
    "hand-rot": "hand_rot.pkl",
    "hand-poses": "hand_rot_trans.pkl",
    "toy-poses": "toy_rot_trans.pkl",
    "toy-angles": "toy_angles.pkl",
    "tilt-poses": "sta_pose.pkl",
}
CTF_FILES = {
    "CTF-Test": "test_ctf.pkl",
    "CTF-Tilt": "sta_ctf.pkl",
}
IND_FILES = {
    "first-100": "ind100.pkl",
    "random-100": "ind100-rand.pkl",
    "just-4": "ind4.pkl",
    "just-5": "ind5.pkl",
}
DATA_FOLDERS = {
    "default-datadir": DATA_DIR,
}


@pytest.fixture(scope="function")
def particles(request) -> Union[None, str]:
    if request.param:
        assert (
            request.param in PARTICLES_FILES
        ), f"Unknown testing particles label `{request.param}` !"
        return os.path.join(DATA_DIR, PARTICLES_FILES[request.param])


@pytest.fixture(scope="function")
def poses(request) -> Union[None, str]:
    if request.param:
        assert (
            request.param in POSES_FILES
        ), f"Unknown testing poses label `{request.param}` !"
        return os.path.join(DATA_DIR, POSES_FILES[request.param])


@pytest.fixture(scope="function")
def ctf(request) -> Union[None, str]:
    if request.param:
        assert (
            request.param in CTF_FILES
        ), f"Unknown testing CTF file label `{request.param}` !"
        return os.path.join(DATA_DIR, CTF_FILES[request.param])


@pytest.fixture(scope="function")
def indices(request) -> Union[None, str]:
    if request.param:
        assert (
            request.param in IND_FILES
        ), f"Unknown testing indices label `{request.param}` !"
        return os.path.join(DATA_DIR, IND_FILES[request.param])


@pytest.fixture(scope="function")
def datadir(request) -> Union[None, str]:
    if request.param:
        assert (
            request.param in DATA_FOLDERS
        ), f"Unknown --datadir path `{request.param}` !"
        return DATA_FOLDERS[request.param]


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

        if "out_lbl" not in train_args:
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
def train_dir(request, tmpdir_factory) -> TrainDir:
    """Run an experiment to generate output; remove this output when finished."""
    args = TrainDir.parse_request(request.param)
    out_lbl = f"train-outs_{request.node.__class__.__name__}"
    args.update(dict(out_lbl=tmpdir_factory.mktemp(out_lbl)))

    tdir = TrainDir(**args)
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
        out_lbl: Optional[str] = None,
        seed: Optional[int] = None,
    ) -> None:
        self.zdim = zdim
        self.dataset = dataset
        self.particles, _ = get_testing_datasets(dataset)

        if out_lbl is None:
            self.outdir = os.path.abspath(f"test-output_{dataset}")
        else:
            self.outdir = os.path.abspath(out_lbl)

        shutil.rmtree(self.outdir, ignore_errors=True)
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
        assert os.path.exists(
            os.path.join(self.outdir, f"weights.{self.epochs - 1}.pkl")
        ), err

    def analyze(self, analysis_epoch: int) -> None:
        out, err = run_command(f"cryodrgn analyze {self.outdir} {analysis_epoch}")
        assert ") Finished in " in out, err
        assert os.path.isdir(
            os.path.join(self.outdir, f"analysis.{analysis_epoch}")
        ), err

    def backproject(self) -> None:
        out_path = os.path.join(self.outdir, "backproject")
        out_fl = os.path.join(out_path, "vol.mrc")
        in_poses = os.path.join(self.outdir, "pose.pkl")

        out, err = run_command(
            f"cryodrgn backproject_voxel {self.particles} "
            f"-o {out_fl} --poses {in_poses} "
        )
        assert "Backprojected 100 images in" in out, err

    def view_config(self) -> None:
        out, err = run_command(f"cryodrgn view_config {self.outdir}")
        assert "'cmd'" in out and "'dataset_args'" in out and "'model_args'" in out, out
        assert err == "", err


@pytest.fixture
def abinit_dir(request, tmpdir_factory) -> AbInitioDir:
    args = AbInitioDir.parse_request(request.param)
    out_lbl = f"abinit-outs_{request.function.__name__}"
    args.update(dict(out_lbl=tmpdir_factory.mktemp(out_lbl)))

    adir = AbInitioDir(**args)
    yield adir
    shutil.rmtree(adir.outdir)
