"""Fixtures used across many unit test modules."""

import pytest
import os
import argparse
import resource
import shutil
import tempfile
import threading
import time
import numpy as np
from contextlib import contextmanager, nullcontext
from functools import lru_cache
from pathlib import Path
from typing import Any, Generator, Optional, Union
from dataclasses import dataclass

from cryodrgn.commands import analyze, train_vae
from cryodrgn.dashboard import app as dash_app
from cryodrgn.dashboard.data import DashboardExperiment, load_experiment
from cryodrgn.utils import run_command

os.environ["NUMEXPR_NUM_THREADS"] = "1"
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "tests", "data")


def pytest_configure():
    pytest.DATADIR = DATA_DIR
    # Playwright's headless Chromium (WebGL/SwiftShader) can segfault on teardown
    # and leave multi-hundred-MB core.* files in the pytest cwd. Tests still pass.
    try:
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    except (ValueError, OSError):
        pass


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


# Data fixtures for cryoDRGN inputs
PARTICLES_FILES = {
    "hand": "hand.mrcs",
    "hand-5": "hand.5.mrcs",
    "hand-tilt": "hand_tilt.mrcs",
    "toy.mrcs": "toy_projections.mrcs",
    "toy.mrcs-999": "toy_projections_0-999.mrcs",
    "toy.star": "toy_projections.star",
    "toydatadir.star": "toy_projections_dir.star",
    "toy.star-13": "toy_projections_13.star",
    "toy.txt": "toy_projections.txt",
    "tilts.star": "sta_testing_bin8.star",
    "csparc_small": "cryosparc_J2_particles_exported.cs",
    "csparc_big": "cryosparc_P12_J24_001_particles.cs",
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
    "CTF-Test.100": "test_ctf.100.pkl",
    "CTF-Tilt": "sta_ctf.pkl",
    "CTF1": "ctf1.pkl",
}
IND_FILES = {
    "first-100": "ind100.pkl",
    "random-100": "ind100-rand.pkl",
    "just-4": "ind4.pkl",
    "just-5": "ind5.pkl",
}
WEIGHTS_FILES = {
    "het": "het_weights.pkl",
}
CONFIG_FILES = {
    "het": "het_config.yaml",
}
DATA_FOLDERS = {
    "default-datadir": ".",
    "toy": "toy_datadir",
}
TRANS_FILES = {
    "toy": "toy_trans.pkl",
    "toy-zero": "toy_zero.pkl",
}


# Data fixtures for cryoDRGN outputs (and sometimes inputs)
VOLUME_FILES = {
    "toy": "toy_projections.mrc",
    "toy-small": "toymodel_small_nocenter.mrc",
    "hand": "hand-vol.mrc",
    "spike": "spike-vol.mrc",
    "empiar": "empiar_10076_7.mrc",
}


@dataclass
class DataFixture:
    label: str
    path: Union[None, str]


def produce_data_fixture(
    data_dict: dict[str, str], labels: str
) -> Union[DataFixture, dict[str, DataFixture]]:
    """Retrieves and parses a request for a fixture defined in a data dictionary."""
    if labels is None:
        files = DataFixture(label="None", path=None)

    else:
        lbls = labels if isinstance(labels, dict) else {None: labels}
        files = dict()

        for k, lbl in lbls.items():
            if lbl in data_dict:
                files[k] = DataFixture(
                    label=lbl, path=os.path.join(DATA_DIR, data_dict[lbl])
                )
            else:
                files[k] = DataFixture(label=lbl, path=os.path.join(DATA_DIR, lbl))

        if not isinstance(labels, dict):
            files = files[None]

    return files


@pytest.fixture(scope="function")
def particles(request) -> Union[DataFixture, dict[str, DataFixture]]:
    return produce_data_fixture(PARTICLES_FILES, request.param)


@pytest.fixture(scope="function")
def poses(request) -> Union[DataFixture, dict[str, DataFixture]]:
    return produce_data_fixture(POSES_FILES, request.param)


@pytest.fixture(scope="function")
def ctf(request) -> Union[DataFixture, dict[str, DataFixture]]:
    return produce_data_fixture(CTF_FILES, request.param)


@pytest.fixture(scope="function")
def indices(request) -> Union[DataFixture, dict[str, DataFixture]]:
    return produce_data_fixture(IND_FILES, request.param)


@pytest.fixture(scope="function")
def datadir(request) -> Union[DataFixture, dict[str, DataFixture]]:
    return produce_data_fixture(DATA_FOLDERS, request.param)


@pytest.fixture(scope="function")
def trans(request) -> Union[DataFixture, dict[str, DataFixture]]:
    return produce_data_fixture(TRANS_FILES, request.param)


@pytest.fixture(scope="function")
def volume(request) -> Union[DataFixture, dict[str, DataFixture]]:
    return produce_data_fixture(VOLUME_FILES, request.param)


@pytest.fixture(scope="function")
def weights(request) -> Union[DataFixture, dict[str, DataFixture]]:
    return produce_data_fixture(WEIGHTS_FILES, request.param)


@pytest.fixture(scope="function")
def configs(request) -> Union[DataFixture, dict[str, DataFixture]]:
    return produce_data_fixture(CONFIG_FILES, request.param)


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
            f"--poses {self.poses} -b 8 --no-amp -n={self.epochs} "
        )
        if self.train_cmd == "train_vae":
            cmd += "--zdim=8 --tdim=16 --tlayers=1 --no-analysis "
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
        if epoch and not 1 <= epoch <= self.epochs:
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
            self.epoch_cleaned(epoch)
            for epoch in list(range(1, self.epochs + 1)) + [None]
        )

    def replace_files(self) -> None:
        for orig_file in os.listdir(self.orig_cache):
            shutil.copy(
                os.path.join(self.orig_cache, orig_file),
                os.path.join(self.outdir, orig_file),
            )

    def train_load_epoch(self, load_epoch: int, train_epochs: int) -> None:
        if not 1 <= load_epoch <= self.epochs:
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
def train_dir(request, tmpdir_factory) -> Generator[TrainDir, None, None]:
    """Run an experiment to generate output; remove this output when finished."""
    args = TrainDir.parse_request(request.param)
    out_lbl = f"train-outs_{request.node.__class__.__name__}"
    args.update(dict(out_lbl=tmpdir_factory.mktemp(out_lbl)))

    tdir = TrainDir(**args)
    yield tdir
    shutil.rmtree(tdir.outdir)


@pytest.fixture(scope="function")
def trained_dir(train_dir: TrainDir) -> Generator[TrainDir, None, None]:
    """Get an experiment that has been run; restore its output when done."""
    yield train_dir
    train_dir.replace_files()


@pytest.fixture(scope="session")
def train_dirs(request) -> Generator[list[TrainDir], None, None]:
    """Run experiments to generate outputs; remove these outputs when finished."""
    tdirs = [TrainDir(**TrainDir.parse_request(req)) for req in request.param]
    yield tdirs
    for tdir in tdirs:
        shutil.rmtree(tdir.outdir)


@pytest.fixture(scope="function")
def trained_dirs(train_dirs) -> Generator[list[TrainDir], None, None]:
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
        train_cmd = "abinit_het_old" if self.zdim > 0 else "abinit_homo_old"

        cmd = (
            f"cryodrgn {train_cmd} {self.particles} -o {self.outdir} "
            f"--num-epochs {self.epochs} --no-window --pretrain 100 "
        )
        if self.zdim > 0:
            cmd += f"--zdim {self.zdim} --no-analysis "
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
            os.path.join(self.outdir, f"weights.{self.epochs}.pkl")
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


@pytest.fixture
def abinit_dir(request, tmpdir_factory) -> AbInitioDir:
    args = AbInitioDir.parse_request(request.param)
    out_lbl = f"abinit-outs_{request.function.__name__}"
    args.update(dict(out_lbl=tmpdir_factory.mktemp(out_lbl)))

    adir = AbInitioDir(**args)
    yield adir
    shutil.rmtree(adir.outdir)


# ---------------------------------------------------------------------------
# Dashboard split-suite shared fixtures
# ---------------------------------------------------------------------------

DASHBOARD_TRAIN_EPOCHS = 3
DASHBOARD_ANALYZE_EPOCH = 2
_DASHBOARD_DEFAULT_TEST_CACHE = os.path.join(
    tempfile.gettempdir(), "cryodrgn_smoke_cache"
)
_DASHBOARD_FIXTURE_SUBDIR = "pytest_dashboard_fixture"

_EXPLORER_ELIGIBLE_PATCH_TARGETS = (
    "cryodrgn.dashboard.routes_analysis.explorer_volumes_eligible",
    "cryodrgn.dashboard.routes_explorer.explorer_volumes_eligible",
    "cryodrgn.dashboard.route_helpers.explorer_volumes_eligible",
)


def _monkeypatch_explorer_volumes_eligible(
    monkeypatch: pytest.MonkeyPatch, *, eligible: bool
) -> None:
    for target in _EXPLORER_ELIGIBLE_PATCH_TARGETS:
        monkeypatch.setattr(target, lambda _e, eligible=eligible: eligible)


def _dashboard_data_dir() -> str:
    return os.path.join(os.path.dirname(__file__), "data")


def _dashboard_is_usable_workdir(workdir: str) -> bool:
    required = [
        os.path.join(workdir, "config.yaml"),
        os.path.join(workdir, f"weights.{DASHBOARD_ANALYZE_EPOCH}.pkl"),
        os.path.join(workdir, f"z.{DASHBOARD_ANALYZE_EPOCH}.pkl"),
        os.path.join(workdir, f"analyze.{DASHBOARD_ANALYZE_EPOCH}", "umap.pkl"),
    ]
    return all(os.path.exists(p) for p in required)


@lru_cache(maxsize=1)
def _dashboard_torch_cuda_kernel_usable() -> bool:
    """True when PyTorch can run kernels on the visible CUDA device.

    ``torch.cuda.is_available()`` alone is insufficient: e.g. an H100 (sm_90)
    with an older PyTorch wheel reports CUDA available but kernel launch fails.
    """
    try:
        import torch
    except ImportError:
        return False
    if not torch.cuda.is_available():
        return False
    try:
        probe = torch.zeros(1, device="cuda")
        del probe
        torch.cuda.synchronize()
        return True
    except RuntimeError:
        return False


@contextmanager
def _dashboard_hide_cuda_devices():
    """Hide GPUs from cryoDRGN CLI entry points that auto-select CUDA."""
    prev = os.environ.get("CUDA_VISIBLE_DEVICES")
    os.environ["CUDA_VISIBLE_DEVICES"] = ""
    try:
        yield
    finally:
        if prev is None:
            os.environ.pop("CUDA_VISIBLE_DEVICES", None)
        else:
            os.environ["CUDA_VISIBLE_DEVICES"] = prev


def _dashboard_train_device_context():
    """Use CPU for fixture training when the local GPU is incompatible with torch."""
    if torch_cuda_reports_available_but_broken():
        return _dashboard_hide_cuda_devices()
    return nullcontext()


def torch_cuda_reports_available_but_broken() -> bool:
    """CUDA visible to PyTorch but kernel launch is known to fail on this device."""
    try:
        import torch
    except ImportError:
        return False
    return torch.cuda.is_available() and not _dashboard_torch_cuda_kernel_usable()


def _dashboard_resolve_fixture_workdir(
    tmp_path_factory: pytest.TempPathFactory,
) -> tuple[str, bool]:
    """Pick workdir for ``dashboard_workdir`` and whether it is a shared cache path.

    Order: ``CRYODRGN_DASHBOARD_TEST_OUTDIR`` → complete default cache → fresh tmp.
    """
    cache_root = os.environ.get("CRYODRGN_DASHBOARD_TEST_OUTDIR")
    if cache_root:
        workdir = os.path.join(cache_root, _DASHBOARD_FIXTURE_SUBDIR)
        os.makedirs(workdir, exist_ok=True)
        return workdir, True

    default_workdir = os.path.join(
        _DASHBOARD_DEFAULT_TEST_CACHE, _DASHBOARD_FIXTURE_SUBDIR
    )
    if _dashboard_is_usable_workdir(default_workdir):
        return default_workdir, True

    return str(tmp_path_factory.mktemp("dashboard_vae")), False


def _dashboard_run_train_and_analyze(workdir: str) -> None:
    data_dir = _dashboard_data_dir()
    with _dashboard_train_device_context():
        parser = argparse.ArgumentParser()
        train_vae.add_args(parser)
        train_args = parser.parse_args(
            [
                os.path.join(data_dir, "hand.mrcs"),
                "-o",
                workdir,
                "--poses",
                os.path.join(data_dir, "hand_rot_trans.pkl"),
                "--ctf",
                os.path.join(data_dir, "test_ctf.100.pkl"),
                "-b",
                "8",
                "--no-amp",
                "-n",
                str(DASHBOARD_TRAIN_EPOCHS),
                "--zdim",
                "4",
                "--tdim",
                "16",
                "--tlayers",
                "1",
                "--seed",
                "0",
                "--no-analysis",
            ]
        )
        train_vae.main(train_args)

        parser = argparse.ArgumentParser()
        analyze.add_args(parser)
        analyze_args = parser.parse_args(
            [workdir, str(DASHBOARD_ANALYZE_EPOCH), "--ksample", "5", "--pc", "2"]
        )
        analyze.main(analyze_args)


@pytest.fixture(scope="session")
def dashboard_workdir(tmp_path_factory: pytest.TempPathFactory) -> str:
    """Trained VAE + analyze outputs for dashboard tests.

      Set ``CRYODRGN_DASHBOARD_TEST_OUTDIR`` to a persistent directory to reuse the
    cached ``pytest_dashboard_fixture`` tree and skip re-training between runs.
    """
    workdir, _shared = _dashboard_resolve_fixture_workdir(tmp_path_factory)
    if _dashboard_is_usable_workdir(workdir):
        return workdir
    _dashboard_run_train_and_analyze(workdir)
    assert _dashboard_is_usable_workdir(
        workdir
    ), f"Dashboard fixture incomplete at {workdir!r}"
    return workdir


@pytest.fixture(scope="session")
def dashboard_experiment(dashboard_workdir: str) -> DashboardExperiment:
    return load_experiment(dashboard_workdir)


def _dashboard_copy_workdir(
    dashboard_workdir: str,
    tmp_path_factory: pytest.TempPathFactory,
    parent_name: str,
) -> str:
    """Copy session workdir into ``parent_name/cryo_out`` (writable per xdist worker)."""
    parent = tmp_path_factory.mktemp(parent_name)
    dst = Path(parent) / "cryo_out"
    shutil.copytree(dashboard_workdir, dst)
    return str(dst)


def _dashboard_strip_landscape_dirs(workdir: str) -> None:
    """Remove ``landscape.N`` trees so copies can add minimal analyze_landscape outputs."""
    for name in list(os.listdir(workdir)):
        if name.startswith("landscape.") and os.path.isdir(os.path.join(workdir, name)):
            shutil.rmtree(os.path.join(workdir, name), ignore_errors=True)


def _dashboard_write_minimal_landscape_volpca(
    workdir: str,
    *,
    epoch: int = DASHBOARD_ANALYZE_EPOCH,
    k: int = 3,
) -> str:
    """``landscape.{epoch}/vol_pca_K.pkl`` + ``kmeansK/`` sketch files for vol PCA UI."""
    from cryodrgn import utils

    land = os.path.join(workdir, f"landscape.{epoch}")
    km = os.path.join(land, f"kmeans{k}")
    os.makedirs(km, exist_ok=True)
    pc = np.array(
        [[0.0, 0.1], [1.0, -0.5], [0.2, 0.3]],
        dtype=np.float64,
    )
    utils.save_pkl(pc, os.path.join(land, f"vol_pca_{k}.pkl"))
    for i in range(1, k + 1):
        p = os.path.join(km, f"vol_{i:03d}.mrc")
        with open(p, "wb"):
            pass
    centers_path = os.path.join(km, "centers_ind.txt")
    with open(centers_path, "w", encoding="utf-8") as fh:
        for row in range(k):
            fh.write(f"{row}\n")
    umap_full = np.array(
        [[0.0, 1.0], [2.0, 3.0], [4.0, 5.0], [6.0, 7.0], [8.0, 9.0]],
        dtype=np.float64,
    )
    utils.save_pkl(umap_full, os.path.join(land, "umap.pkl"))
    return land


def _dashboard_write_minimal_landscape_full(
    workdir: str,
    *,
    epoch: int = DASHBOARD_ANALYZE_EPOCH,
    n_sampled: int = 25,
) -> str:
    """``landscape.{epoch}/landscape_full`` outputs for ``analyze_landscape_full`` UI."""
    from cryodrgn import utils

    z_path = os.path.join(workdir, f"z.{epoch}.pkl")
    z_all = np.asarray(utils.load_pkl(z_path), dtype=np.float64)
    if z_all.ndim != 2:
        raise ValueError(f"expected 2-D z array at {z_path}")
    n_particles = int(z_all.shape[0])
    outdir = os.path.join(workdir, f"landscape.{epoch}", "landscape_full")
    os.makedirs(outdir, exist_ok=True)
    n_take = min(int(n_sampled), n_particles)
    ind = np.linspace(0, n_particles - 1, num=n_take, dtype=np.int64)
    z_s = z_all[ind]
    utils.save_pkl(ind, os.path.join(outdir, "ind.sampled.pkl"))
    utils.save_pkl(z_s, os.path.join(outdir, "z.sampled.pkl"))
    rng = np.random.default_rng(0)
    vol_pca = rng.standard_normal((n_particles, 5), dtype=np.float64)
    utils.save_pkl(vol_pca, os.path.join(outdir, "vol_pca_all.pkl"))
    return outdir


@pytest.fixture(scope="module")
def dashboard_workdir_with_landscape_volpca(
    dashboard_workdir: str,
    tmp_path_factory: pytest.TempPathFactory,
) -> str:
    """Copy shared output and add minimal ``analyze_landscape`` sketch files."""
    d = _dashboard_copy_workdir(
        dashboard_workdir, tmp_path_factory, "workdir_landscape"
    )
    _dashboard_strip_landscape_dirs(d)
    _dashboard_write_minimal_landscape_volpca(d)
    return d


@pytest.fixture(scope="module")
def dashboard_workdir_plain_copy(
    dashboard_workdir: str,
    tmp_path_factory: pytest.TempPathFactory,
) -> str:
    """Writable copy with no ``landscape.*`` trees."""
    d = _dashboard_copy_workdir(
        dashboard_workdir, tmp_path_factory, "wd_plain_no_landscape"
    )
    _dashboard_strip_landscape_dirs(d)
    return d


@pytest.fixture(scope="module")
def dashboard_workdir_with_landscape_full(
    dashboard_workdir: str,
    tmp_path_factory: pytest.TempPathFactory,
) -> str:
    """Dashboard tree with mocked ``analyze_landscape`` + ``analyze_landscape_full`` outputs."""
    d = _dashboard_copy_workdir(
        dashboard_workdir, tmp_path_factory, "workdir_landscape_full"
    )
    _dashboard_strip_landscape_dirs(d)
    _dashboard_write_minimal_landscape_volpca(d)
    _dashboard_write_minimal_landscape_full(d)
    return d


@contextmanager
def _patch_explorer_volumes_eligible(eligible: bool = True):
    """Force trajectory / volume-explorer eligibility (no CUDA required in tests)."""
    from unittest.mock import patch

    target = lambda _e, eligible=eligible: eligible  # noqa: E731
    patches = [patch(t, target) for t in _EXPLORER_ELIGIBLE_PATCH_TARGETS]
    for p in patches:
        p.start()
    try:
        yield
    finally:
        for p in patches:
            p.stop()


@pytest.fixture(scope="module")
def dashboard_landscape_volpca_live_url(dashboard_workdir_with_landscape_volpca: str):
    with dashboard_live_server(dashboard_workdir_with_landscape_volpca) as url:
        yield url


@pytest.fixture(scope="module")
def dashboard_landscape_full_live_url(dashboard_workdir_with_landscape_full: str):
    with dashboard_live_server(dashboard_workdir_with_landscape_full) as url:
        yield url


@pytest.fixture(scope="module")
def dashboard_volumes_eligible_live_url(dashboard_workdir: str):
    with _patch_explorer_volumes_eligible(True):
        with dashboard_live_server(dashboard_workdir) as url:
            yield url


@pytest.fixture(scope="module")
def dashboard_volumes_ineligible_live_url(dashboard_workdir: str):
    with _patch_explorer_volumes_eligible(False):
        with dashboard_live_server(dashboard_workdir) as url:
            yield url


@pytest.fixture
def flask_client_landscape_full(dashboard_workdir_with_landscape_full: str):
    app = dash_app.create_app(workdir=dashboard_workdir_with_landscape_full)
    with app.test_client() as client:
        yield client


@pytest.fixture
def flask_client_landscape(dashboard_workdir_with_landscape_volpca: str):
    app = dash_app.create_app(workdir=dashboard_workdir_with_landscape_volpca)
    with app.test_client() as client:
        yield client


@pytest.fixture
def flask_client_no_landscape(dashboard_workdir_plain_copy: str):
    app = dash_app.create_app(workdir=dashboard_workdir_plain_copy)
    with app.test_client() as client:
        yield client


@pytest.fixture
def flask_client_volumes_eligible(
    dashboard_workdir: str, monkeypatch: pytest.MonkeyPatch
) -> Generator:
    _monkeypatch_explorer_volumes_eligible(monkeypatch, eligible=True)
    app = dash_app.create_app(workdir=dashboard_workdir)
    with app.test_client() as client:
        yield client


@pytest.fixture(scope="function")
def flask_client(dashboard_workdir: str):
    app = dash_app.create_app(workdir=dashboard_workdir)
    with app.test_client() as client:
        yield client


# ---------------------------------------------------------------------------
# Shared helpers for dashboard test modules
# ---------------------------------------------------------------------------


def is_plotly_typed_array(value: Any) -> bool:
    return isinstance(value, dict) and "dtype" in value and "bdata" in value


def _plotly_typed_array_shape_dims(shape: Any) -> tuple[int, ...]:
    """Normalize Plotly 6 ``shape`` (comma string or JSON list) for ``numpy.reshape``."""
    if shape is None:
        return ()
    if isinstance(shape, (list, tuple)):
        return tuple(int(d) for d in shape)
    parts = [s.strip() for s in str(shape).split(",") if s.strip()]
    return tuple(int(s) for s in parts)


def decode_plotly_value(value: Any) -> Any:
    """Expand Plotly 6 typed-array blobs (``{dtype, bdata[, shape]}``) to Python lists."""
    if is_plotly_typed_array(value):
        import base64

        import numpy as np

        raw = base64.b64decode(value["bdata"])
        arr = np.frombuffer(raw, dtype=np.dtype(value["dtype"]))
        dims = _plotly_typed_array_shape_dims(value.get("shape"))
        if dims:
            arr = arr.reshape(dims)
        return arr.tolist()
    if isinstance(value, list):
        return [decode_plotly_value(v) for v in value]
    if isinstance(value, dict):
        return {k: decode_plotly_value(v) for k, v in value.items()}
    return value


def decode_plotly_figure(fig: dict[str, Any]) -> dict[str, Any]:
    return decode_plotly_value(fig)


def plotly_trace_array(trace: dict[str, Any], key: str) -> Any:
    return decode_plotly_value(trace[key])


def dashboard_repo_root() -> "Path":
    from pathlib import Path

    return Path(__file__).resolve().parents[1]


def png_b64_rgb(
    w: int = 8, h: int = 8, rgb: tuple[int, int, int] = (200, 30, 40)
) -> str:
    """Solid-colour PNG as standard base64 (GIF API + plot_gif_utils tests)."""
    import base64
    from io import BytesIO

    from PIL import Image

    buf = BytesIO()
    Image.new("RGB", (w, h), rgb).save(buf, format="PNG")
    return base64.standard_b64encode(buf.getvalue()).decode("ascii")


@lru_cache(maxsize=None)
def read_dashboard_template(rel: str) -> str:
    """Cached read of a Jinja template under ``cryodrgn/dashboard/templates/``."""
    path = dashboard_repo_root() / "cryodrgn" / "dashboard" / "templates" / rel
    return path.read_text(encoding="utf-8")


def _read_dashboard_static(rel: str, subdir: str) -> str:
    path = dashboard_repo_root() / "cryodrgn" / "dashboard" / "static" / subdir / rel
    return path.read_text(encoding="utf-8")


@lru_cache(maxsize=None)
def read_dashboard_static_js(rel: str) -> str:
    """Cached read of ``cryodrgn/dashboard/static/js/<rel>``."""
    return _read_dashboard_static(rel, "js")


@lru_cache(maxsize=None)
def read_dashboard_static_css(rel: str) -> str:
    """Cached read of ``cryodrgn/dashboard/static/css/<rel>``."""
    return _read_dashboard_static(rel, "css")


def read_latent_3d_html() -> str:
    """Cached read of ``latent_3d.html``."""
    return read_dashboard_template("latent_3d.html")


def js_function_body(source: str, fn_marker: str, until_marker: str) -> str:
    """Slice a JS source string from ``fn_marker`` up to (but not including) ``until_marker``."""
    return source.split(fn_marker, 1)[1].split(until_marker, 1)[0]


# ---------------------------------------------------------------------------
# Dashboard Plotly browser smoke (Playwright; optional ``dev`` extra)
# ---------------------------------------------------------------------------

DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS = 120_000
DASHBOARD_BROWSER_FAST_TIMEOUT_MS = 30_000
DASHBOARD_BROWSER_SMOKE_CACHE_SIZE = 25


def make_fake_trajectory_volume_pngs(cache_token: str = "cache-tok"):
    """Return a ``generate_trajectory_volume_pngs`` stub accepting ``**kwargs``."""

    def _fake(exp, z_values, **kwargs):
        import io

        from PIL import Image

        blobs = []
        for _ in range(len(z_values)):
            buf = io.BytesIO()
            Image.new("RGB", (4, 4)).save(buf, format="PNG")
            blobs.append(buf.getvalue())
        return blobs, cache_token

    return _fake


def _volume_viewer_stub_payloads():
    """JSON bodies for stubbed ``/api/volume_viewer/*`` responses."""
    import json

    tiny = png_b64_rgb()
    catalog = [
        {
            "id": "kmeans:0",
            "kind": "kmeans",
            "cluster_label": 0,
            "znorm": 0.1,
        },
        {
            "id": "kmeans:1",
            "kind": "kmeans",
            "cluster_label": 1,
            "znorm": 0.2,
        },
        {
            "id": "kmeans:2",
            "kind": "kmeans",
            "cluster_label": 2,
            "znorm": 0.3,
        },
    ]
    markers = [
        {
            "vol_id": f"kmeans:{i}",
            "plot_row": i * 10,
            "label": f"K{i + 1}",
            "kind": "kmeans",
        }
        for i in range(3)
    ]
    return {
        "tiny": tiny,
        "catalog_with_markers": json.dumps(
            {
                "ok": True,
                "catalog": catalog,
                "default_vol_id": "kmeans:0",
                "markers": markers,
                "D": None,
            }
        ),
        "catalog_only": json.dumps(
            {
                "ok": True,
                "catalog": catalog,
                "default_vol_id": "kmeans:0",
                "markers": [],
                "D": None,
            }
        ),
        "markers_only": json.dumps({"ok": True, "markers": markers}),
        "single_volume": json.dumps(
            {"ok": True, "id": "kmeans:0", "volume_b64": tiny, "D": 32}
        ),
    }


def fulfill_volume_viewer_render_route(route) -> bool:
    """Fulfill slow volume-viewer render POSTs; return True if handled."""
    import json

    tiny = _volume_viewer_stub_payloads()["tiny"]
    url = route.request.url
    method = route.request.method
    if "analyze_volumes_chimerax_batch" in url and method == "POST":
        try:
            req = json.loads(route.request.post_data or "{}")
        except json.JSONDecodeError:
            req = {}
        ids = req.get("ids") or ["kmeans:0"]
        images = [tiny] * max(1, len(ids))
        route.fulfill(
            status=200,
            content_type="application/json",
            body=json.dumps({"ok": True, "images": images}),
        )
        return True
    if "analyze_volumes_batch" in url and method == "POST":
        try:
            req = json.loads(route.request.post_data or "{}")
        except json.JSONDecodeError:
            req = {}
        ids = req.get("ids") or ["kmeans:0"]
        volumes = {
            str(vol_id): {"volume_b64": tiny, "D": 32, "id": str(vol_id)}
            for vol_id in ids
        }
        route.fulfill(
            status=200,
            content_type="application/json",
            body=json.dumps({"ok": True, "volumes": volumes}),
        )
        return True
    return False


def playwright_route_volume_viewer_render_stub(page) -> None:
    """Stub only slow ChimeraX/VTK batch renders; catalog GETs hit the live server."""

    def _handler(route):
        if not fulfill_volume_viewer_render_route(route):
            route.continue_()

    page.route("**/api/volume_viewer/**", _handler)


def fulfill_volume_viewer_route(route) -> bool:
    """Fulfill a Playwright route with volume-viewer stubs; return True if handled."""
    if fulfill_volume_viewer_render_route(route):
        return True

    payloads = _volume_viewer_stub_payloads()
    url = route.request.url
    method = route.request.method
    if "analyze_volumes" in url and method == "GET":
        body = (
            payloads["catalog_only"]
            if "include_markers=0" in url
            else payloads["catalog_with_markers"]
        )
        route.fulfill(status=200, content_type="application/json", body=body)
        return True
    if "analyze_markers" in url and method == "GET":
        route.fulfill(
            status=200,
            content_type="application/json",
            body=payloads["markers_only"],
        )
        return True
    if "analyze_volume" in url and method == "GET":
        route.fulfill(
            status=200,
            content_type="application/json",
            body=payloads["single_volume"],
        )
        return True
    return False


def playwright_route_volume_viewer_stub(page) -> None:
    """Intercept ``/api/volume_viewer/*`` with fast JSON stubs for Playwright tests."""

    def _handler(route):
        if not fulfill_volume_viewer_route(route):
            route.continue_()

    page.route("**/api/volume_viewer/**", _handler)


@contextmanager
def dashboard_live_server(workdir: str) -> Generator[str, None, None]:
    """Threaded Werkzeug server for headless browser tests against a real workdir."""
    from werkzeug.serving import make_server

    app = dash_app.create_app(workdir=workdir)
    server = make_server("127.0.0.1", 0, app)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{server.server_port}"
    finally:
        server.shutdown()
        thread.join(timeout=10)


@pytest.fixture(scope="module")
def dashboard_live_url(dashboard_workdir: str) -> Generator[str, None, None]:
    with dashboard_live_server(dashboard_workdir) as url:
        yield url


@pytest.fixture(scope="module")
def dashboard_plain_live_url(dashboard_workdir_plain_copy: str):
    with dashboard_live_server(dashboard_workdir_plain_copy) as url:
        yield url


@pytest.fixture(scope="module")
def playwright_browser():
    pytest.importorskip("playwright")
    from playwright.sync_api import sync_playwright

    with sync_playwright() as p:
        browser = _launch_playwright_chromium(p)
        yield browser
        browser.close()


@pytest.fixture(scope="module")
def playwright_page(playwright_browser):
    context = playwright_browser.new_context(viewport={"width": 1400, "height": 900})
    page = context.new_page()
    yield page
    context.close()


def _playwright_chromium_launch_kwargs() -> dict:
    return {
        "headless": True,
        "args": [
            "--enable-unsafe-swiftshader",
            "--ignore-gpu-blocklist",
            "--enable-webgl",
            "--max-active-webgl-contexts=32",
        ],
    }


def _launch_playwright_chromium(playwright):
    """Launch headless Chromium; SwiftShader helps WebGL on GPU-less CI nodes."""
    from playwright.sync_api import Error as PlaywrightError

    try:
        return playwright.chromium.launch(**_playwright_chromium_launch_kwargs())
    except PlaywrightError as exc:
        pytest.skip(
            f"Playwright Chromium not installed ({exc}). "
            "Run: playwright install chromium"
        )


def _dashboard_smoke_scene_camera(page, plot_id: str) -> dict | None:
    return page.evaluate(
        """(plotId) => {
          var gd = document.getElementById(plotId);
          var P3S = window.CryoPlotlyScatter3dScene;
          var snap = gd && P3S && P3S.snapshot(gd);
          if (!snap || !snap.camera) return null;
          var cam = snap.camera;
          function vec(v) { return v ? [v.x, v.y, v.z] : null; }
          return { eye: vec(cam.eye), center: vec(cam.center), up: vec(cam.up) };
        }""",
        plot_id,
    )


def _dashboard_smoke_wait_scene_camera(page, plot_id: str, *, timeout_ms: int) -> dict:
    page.wait_for_function(
        """(plotId) => {
          var gd = document.getElementById(plotId);
          var P3S = window.CryoPlotlyScatter3dScene;
          var snap = gd && P3S && P3S.snapshot(gd);
          return !!(snap && snap.camera && snap.camera.eye);
        }""",
        arg=plot_id,
        timeout=timeout_ms,
    )
    cam = _dashboard_smoke_scene_camera(page, plot_id)
    if not cam or not cam.get("eye"):
        raise RuntimeError(f"plot #{plot_id} camera missing after wait: {cam!r}")
    return cam


def _dashboard_smoke_cameras_match(
    a: dict | None, b: dict | None, *, tol: float = 1e-4
) -> bool:
    if not a or not b:
        return False
    for key in ("eye", "center", "up"):
        va, vb = a.get(key), b.get(key)
        if not va or not vb or len(va) != 3 or len(vb) != 3:
            return False
        if any(abs(x - y) > tol for x, y in zip(va, vb)):
            return False
    return True


def _dashboard_smoke_wait_latent3d_overlay_hidden(page, *, timeout_ms: int) -> None:
    page.wait_for_function(
        """() => {
          var ov = document.getElementById('latent3d-rendering-overlay');
          if (!ov) return true;
          var st = window.getComputedStyle(ov);
          return st.display === 'none' || st.visibility === 'hidden' || st.opacity === '0';
        }""",
        timeout=timeout_ms,
    )


def _dashboard_smoke_set_select_value(page, select_id: str, value: str) -> None:
    """Set a native ``<select>`` and fire ``change`` (works when options are CSS-hidden)."""
    page.evaluate(
        """([id, val]) => {
          var el = document.getElementById(id);
          if (!el) throw new Error('missing select #' + id);
          el.value = val;
          el.dispatchEvent(new Event('change', { bubbles: true }));
        }""",
        [select_id, value],
    )


def _dashboard_smoke_wait_plot_ready(page, plot_id: str, *, timeout_ms: int) -> dict:
    page.wait_for_selector(
        f"#{plot_id} .js-plotly-plot, #{plot_id} .main-svg, #{plot_id} .gl-canvas",
        timeout=timeout_ms,
    )
    deadline = time.time() + timeout_ms / 1000.0
    last: dict = {}
    while time.time() < deadline:
        last = page.evaluate(
            """(plotId) => {
              var gd = document.getElementById(plotId);
              var P = window.CryoPlotlyArrays;
              if (!gd || !P) return { ok: false, reason: 'missing gd or CryoPlotlyArrays' };
              if (!gd.data || !gd.data[0]) return { ok: false, reason: 'no trace yet' };
              var tr = gd.data[0];
              var n = P.length(tr.x);
              return { ok: n > 0, n: n, type: tr.type || '' };
            }""",
            plot_id,
        )
        if last.get("ok"):
            return last
        page.wait_for_timeout(250)
    raise RuntimeError(f"plot #{plot_id} not ready: {last!r}")


def _dashboard_smoke_preload_overlay_hidden(
    page, *, timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS
) -> None:
    """Wait until the montage preload overlay is dismissed (class + aria-hidden)."""
    page.wait_for_function(
        """() => {
          var ov = document.getElementById('montage-preload-overlay');
          if (!ov) return true;
          return ov.getAttribute('aria-hidden') === 'true'
            && !ov.classList.contains('cryo-plot-rendering-overlay--show');
        }""",
        timeout=timeout_ms,
    )


def _dashboard_smoke_invalidate_server_preload_cache(page) -> None:
    """Drop Flask ``PRELOAD_CACHE`` for this browser session (survives ``page.goto``)."""
    page.evaluate(
        """() => fetch('/api/preload_images', {
          method: 'POST',
          headers: {'Content-Type': 'application/json'},
          body: JSON.stringify({invalidate_cache: true})
        }).then(function(r) { return r.json(); })"""
    )


def _dashboard_smoke_clear_explorer_cache_if_any(page, *, timeout_ms: int) -> None:
    """Drop explorer preload state (server + any client UI cache)."""
    _dashboard_smoke_invalidate_server_preload_cache(page)
    _dashboard_smoke_ensure_panel_open(page, "cache-panel-toggle")
    has_cache = page.evaluate(
        """() => {
          var st = document.getElementById('preload-status');
          return !!(st && /\\bcached\\b/i.test(st.textContent || ''));
        }"""
    )
    if not has_cache:
        return
    clear_btn = page.locator("#btn-clear-image-cache")
    clear_btn.wait_for(state="visible", timeout=timeout_ms)
    if clear_btn.is_disabled():
        _dashboard_smoke_preload_overlay_hidden(page, timeout_ms=timeout_ms)
        if clear_btn.is_disabled():
            return
    clear_btn.click()
    page.wait_for_function(
        """() => {
          var st = document.getElementById('preload-status');
          var txt = st ? st.textContent || '' : '';
          return !/\\bcached\\b/i.test(txt);
        }""",
        timeout=timeout_ms,
    )


def _dashboard_smoke_ensure_panel_open(page, toggle_id: str) -> None:
    expanded = page.evaluate(
        """(id) => {
          var b = document.getElementById(id);
          return b && b.getAttribute('aria-expanded') === 'true';
        }""",
        toggle_id,
    )
    if not expanded:
        page.click(f"#{toggle_id}")


def dashboard_smoke_particle_explorer(
    page,
    base_url: str,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
    cache_size: int | None = DASHBOARD_BROWSER_SMOKE_CACHE_SIZE,
) -> dict:
    base = base_url.rstrip("/")
    page.goto(f"{base}/explorer", wait_until="domcontentloaded", timeout=timeout_ms)
    info = _dashboard_smoke_wait_plot_ready(page, "scatter", timeout_ms=timeout_ms)

    _dashboard_smoke_clear_explorer_cache_if_any(page, timeout_ms=timeout_ms)
    _dashboard_smoke_ensure_panel_open(page, "cache-panel-toggle")
    if cache_size is not None:
        page.fill("#montage-cache-size-input", str(cache_size))
    btn = page.locator("#btn-expand-cache")
    btn.wait_for(state="visible", timeout=30_000)
    if btn.is_disabled():
        meta = page.evaluate(
            """() => {
              var b = document.getElementById('btn-expand-cache');
              return { title: b && b.title, text: b && b.textContent };
            }"""
        )
        raise RuntimeError(f"Build cache disabled: {meta!r}")

    btn.click()
    _dashboard_smoke_preload_overlay_hidden(page, timeout_ms=timeout_ms)
    page.wait_for_function(
        """() => {
          var st = document.getElementById('preload-status');
          return st && /\\bcached\\b/i.test(st.textContent || '');
        }""",
        timeout=timeout_ms,
    )

    preload = page.evaluate(
        """() => {
          var st = document.getElementById('preload-status');
          var txt = st ? st.textContent.trim() : '';
          var m = txt.match(/([\\d,]+)\\s+images\\s+cached/i);
          return { status: txt, cached: m ? parseInt(m[1].replace(/,/g, ''), 10) : 0 };
        }"""
    )
    if not preload.get("cached"):
        raise RuntimeError(f"cache not populated: {preload!r}")

    _dashboard_smoke_ensure_panel_open(page, "image-grid-menu-toggle")
    view_btn = page.locator("#btn-view-images")
    if view_btn.is_disabled():
        raise RuntimeError("Load all images button still disabled after cache build")
    view_btn.click()

    page.wait_for_selector("#montage-grid img[src]:not([src=''])", timeout=timeout_ms)
    img_count = page.locator("#montage-grid img[src]:not([src=''])").count()

    resample = page.locator("#btn-montage-resample-cache")
    if resample.count() and not resample.is_disabled():
        resample.click()
        page.wait_for_function(
            """() => {
              var gd = document.getElementById('scatter');
              return gd && gd.data && gd.data.length > 1 && gd.data[1] && gd.data[1].x;
            }""",
            timeout=30_000,
        )

    letters = page.evaluate(
        """() => {
          var glyphs = document.querySelectorAll('.cryo-explorer-grid-letter-glyph');
          var sample = Array.from(glyphs).slice(0, 5).map(function(el) { return el.textContent; });
          return { count: glyphs.length, sample: sample.join('') };
        }"""
    )
    return {
        "scatter_points": info.get("n"),
        "cached_images": preload.get("cached"),
        "grid_images": img_count,
        "scatter_letters": letters,
    }


def dashboard_smoke_landscape_volpca(
    page,
    base_url: str,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
) -> dict | None:
    base = base_url.rstrip("/")
    page.goto(
        f"{base}/landscape-volpca", wait_until="domcontentloaded", timeout=timeout_ms
    )
    if page.locator("#volsketch").count() == 0:
        body = page.inner_text("body")
        if "landscape" in body.lower() or "analyze_landscape" in body:
            return None
        raise RuntimeError("landscape-volpca page missing #volsketch")
    info = _dashboard_smoke_wait_plot_ready(page, "volsketch", timeout_ms=timeout_ms)

    page.wait_for_function(
        """() => {
          var b = document.getElementById('volsketch-random-sel');
          return b && !b.disabled;
        }""",
        timeout=timeout_ms,
    )
    page.click("#volsketch-random-sel")
    page.wait_for_function(
        """() => {
          var gd = document.getElementById('volsketch');
          if (!gd || !gd.data) return false;
          for (var i = 0; i < gd.data.length; i++) {
            if (gd.data[i] && gd.data[i].name === 'cdrgnVolSelectionOverlay') return true;
          }
          return false;
        }""",
        timeout=30_000,
    )

    overlay = page.evaluate(
        """() => {
          var gd = document.getElementById('volsketch');
          var P = window.CryoPlotlyArrays;
          if (!gd || !P || !gd.data) return { traces: 0 };
          var overlayIdx = -1;
          for (var i = 0; i < gd.data.length; i++) {
            if (gd.data[i] && gd.data[i].name === 'cdrgnVolSelectionOverlay') overlayIdx = i;
          }
          var tr = overlayIdx >= 0 ? gd.data[overlayIdx] : null;
          var texts = tr && tr.text ? P.length(tr.text) : 0;
          return { overlayIdx: overlayIdx, overlayTexts: texts };
        }"""
    )
    if overlay.get("overlayIdx", -1) < 0 or overlay.get("overlayTexts", 0) < 1:
        raise RuntimeError(f"vol selection overlay missing letters: {overlay!r}")
    return {"scatter_points": info.get("n"), **overlay}


def dashboard_smoke_landscape_full_3d(
    page,
    base_url: str,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
) -> dict | None:
    base = base_url.rstrip("/")
    page.goto(
        f"{base}/landscape-full-3d", wait_until="domcontentloaded", timeout=timeout_ms
    )
    if page.locator("#latent3d").count() == 0:
        body = page.inner_text("body")
        if "analyze_landscape_full" in body or "need" in body.lower():
            return None
        raise RuntimeError("unexpected landscape-full-3d page without #latent3d")

    info = _dashboard_smoke_wait_plot_ready(page, "latent3d", timeout_ms=timeout_ms)
    page.wait_for_function(
        """() => {
          var b = document.getElementById('l3dva-random-sel');
          return b && !b.disabled;
        }""",
        timeout=timeout_ms,
    )
    page.click("#l3dva-random-sel")
    page.wait_for_function(
        """() => {
          var gd = document.getElementById('latent3d');
          var anns = gd && gd.layout && gd.layout.scene && gd.layout.scene.annotations;
          return anns && anns.length > 0;
        }""",
        timeout=timeout_ms,
    )
    ann = page.evaluate(
        """() => {
          var gd = document.getElementById('latent3d');
          var anns = (gd && gd.layout && gd.layout.scene && gd.layout.scene.annotations) || [];
          return {
            count: anns.length,
            sample: anns.slice(0, 3).map(function(a) { return String(a.text || ''); }),
          };
        }"""
    )
    if ann.get("count", 0) < 1:
        raise RuntimeError(f"no scene annotations after random vol selection: {ann!r}")
    return {"scatter_points": info.get("n"), **ann}


def _dashboard_smoke_eval_plot(page, plot_id: str) -> dict:
    return page.evaluate(
        """(plotId) => {
          var gd = document.getElementById(plotId);
          var P = window.CryoPlotlyArrays;
          if (!gd || !P || !gd.data || !gd.data[0]) {
            return { traces: 0, points: 0, trace0type: '' };
          }
          return {
            traces: gd.data.length,
            points: P.length(gd.data[0].x),
            trace0type: String(gd.data[0].type || ''),
          };
        }""",
        plot_id,
    )


def _dashboard_smoke_wait_pairplot_image(page, *, timeout_ms: int) -> str:
    page.wait_for_function(
        """() => {
          var img = document.getElementById('pairplot');
          var vp = document.getElementById('pairplot-viewport');
          return img && img.src && img.src.length > 32
            && vp && vp.getAttribute('aria-busy') === 'false';
        }""",
        timeout=timeout_ms,
    )
    src = page.evaluate(
        """() => {
      var img = document.getElementById('pairplot');
      return img ? String(img.src || '') : '';
    }"""
    )
    if not src:
        raise RuntimeError("pairplot image missing src")
    return src


def dashboard_smoke_index(
    page, base_url: str, *, timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS
) -> dict:
    base = base_url.rstrip("/")
    page.goto(f"{base}/", wait_until="domcontentloaded", timeout=timeout_ms)
    page.wait_for_selector(".landing-cards", timeout=timeout_ms)
    links = page.evaluate(
        """() => Array.from(document.querySelectorAll('a.landing-card-link'))
          .map(function(a) { return a.getAttribute('href') || ''; })"""
    )
    required_cards = {"/explorer", "/pairplot", "/latent-3d", "/command-builder"}
    missing_cards = sorted(required_cards - set(links))
    if missing_cards:
        raise RuntimeError(f"landing page missing card links: {missing_cards}")
    body_text = page.inner_text("body")
    for link in required_cards:
        if link not in body_text and link not in links:
            raise RuntimeError(f"landing page missing reference to {link!r}")
    return {
        "landing_links": len(links),
        "card_links": links,
        "has_trajectory": "/trajectory" in links,
        "has_landscape_volpca": "/landscape-volpca" in links,
        "mentions_volume_landscapes": "3D volume landscapes" in body_text,
        "trajectory_ineligible_note": "CUDA-enabled machine" in body_text,
    }


def dashboard_smoke_index_no_landscape(
    page,
    base_url: str,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
) -> dict:
    out = dashboard_smoke_index(page, base_url, timeout_ms=timeout_ms)
    if out["has_landscape_volpca"]:
        raise RuntimeError(
            "expected inactive landscape card without analyze_landscape outputs"
        )
    body_text = page.inner_text("body")
    if "Volume sketched landscape explorer" not in body_text:
        raise RuntimeError("missing inactive landscape card title")
    if "analyze_landscape" not in body_text:
        raise RuntimeError("missing analyze_landscape hint on inactive landscape card")
    return out


def dashboard_smoke_latent_3d_camera_on_covariate_change(
    page,
    base_url: str,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
) -> dict:
    """Discrete colour-covariate reload must not reset the latent-3d orbit."""
    base = base_url.rstrip("/")
    page.goto(f"{base}/latent-3d", wait_until="domcontentloaded", timeout=timeout_ms)
    _dashboard_smoke_wait_plot_ready(page, "latent3d", timeout_ms=timeout_ms)
    _dashboard_smoke_wait_latent3d_overlay_hidden(page, timeout_ms=timeout_ms)
    cam_before = _dashboard_smoke_wait_scene_camera(
        page, "latent3d", timeout_ms=timeout_ms
    )

    page.wait_for_function(
        """() => {
          var el = document.getElementById('sc');
          return el && el.options && el.options.length > 1;
        }""",
        timeout=timeout_ms,
    )
    _dashboard_smoke_set_select_value(page, "sc", "labels")
    page.wait_for_function(
        """() => {
          var switches = document.getElementById('latent3d-color-discrete-switches');
          return switches && switches.querySelectorAll('button, label, input').length > 0;
        }""",
        timeout=timeout_ms,
    )
    _dashboard_smoke_wait_latent3d_overlay_hidden(page, timeout_ms=timeout_ms)
    cam_after = _dashboard_smoke_wait_scene_camera(
        page, "latent3d", timeout_ms=timeout_ms
    )
    stable = _dashboard_smoke_cameras_match(cam_before, cam_after)
    if not stable:
        raise RuntimeError(
            f"latent3d camera changed after covariate switch: {cam_before!r} -> {cam_after!r}"
        )
    return {
        "camera_stable": stable,
        "eye_before": cam_before["eye"],
        "eye_after": cam_after["eye"],
    }


def dashboard_smoke_latent_3d(
    page,
    base_url: str,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
) -> dict:
    base = base_url.rstrip("/")
    page.goto(f"{base}/latent-3d", wait_until="domcontentloaded", timeout=timeout_ms)
    info = _dashboard_smoke_wait_plot_ready(page, "latent3d", timeout_ms=timeout_ms)
    plot0 = _dashboard_smoke_eval_plot(page, "latent3d")
    if plot0.get("trace0type") != "scatter3d":
        raise RuntimeError(f"expected scatter3d trace, got {plot0!r}")

    page.wait_for_function(
        """() => {
          var el = document.getElementById('sc');
          return el && el.options && el.options.length > 1;
        }""",
        timeout=timeout_ms,
    )
    _dashboard_smoke_set_select_value(page, "sc", "labels")
    page.wait_for_function(
        """() => {
          var switches = document.getElementById('latent3d-color-discrete-switches');
          return switches && switches.querySelectorAll('button, label, input').length > 0;
        }""",
        timeout=timeout_ms,
    )
    discrete_n = page.locator(
        "#latent3d-color-discrete-switches button, #latent3d-color-discrete-switches label"
    ).count()
    return {
        "scatter_points": info.get("n"),
        "discrete_legend_toggles": discrete_n,
        **plot0,
    }


def dashboard_smoke_pairplot(
    page,
    base_url: str,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
) -> dict:
    base = base_url.rstrip("/")
    page.goto(f"{base}/pairplot", wait_until="domcontentloaded", timeout=timeout_ms)
    src0 = _dashboard_smoke_wait_pairplot_image(page, timeout_ms=timeout_ms)
    page.click('input[name="upper_style"][value="hex"]')
    page.wait_for_function(
        """(prev) => {
          var img = document.getElementById('pairplot');
          return img && img.src && img.src !== prev;
        }""",
        arg=src0,
        timeout=timeout_ms,
    )
    src1 = page.evaluate("() => document.getElementById('pairplot').src")
    return {
        "initial_src_len": len(src0),
        "hex_src_len": len(src1),
        "src_changed": src1 != src0,
    }


def dashboard_smoke_trajectory(
    page,
    base_url: str,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
) -> dict | None:
    base = base_url.rstrip("/")
    page.goto(f"{base}/trajectory", wait_until="domcontentloaded", timeout=timeout_ms)
    if page.locator("#scatter").count() == 0:
        body = page.inner_text("body")
        if "CUDA" in body or "GPU" in body or "weights" in body.lower():
            return None
        raise RuntimeError(
            "trajectory page missing #scatter without known ineligible message"
        )
    info = _dashboard_smoke_wait_plot_ready(page, "scatter", timeout_ms=timeout_ms)
    page.wait_for_function(
        """() => {
          var overlay = document.getElementById('traj-glyph-overlay');
          var picker = document.getElementById('vslice-volume-picker-rows');
          var hasGlyph = !!(overlay && overlay.querySelector('.cryo-traj-glyph-path'));
          var pickerBtns = picker
            ? picker.querySelectorAll('.cryo-vslice-vol-btn').length
            : 0;
          return hasGlyph && pickerBtns >= 2;
        }""",
        timeout=timeout_ms,
    )
    before = _dashboard_smoke_trajectory_state(page)
    page.click("#btn-anchor-random")
    page.wait_for_function(
        """() => {
          var st = document.getElementById('traj-status');
          var txt = st ? st.textContent || '' : '';
          return /manual selection ready|anchor path|latent z ready|selection and volumes ready/i.test(txt);
        }""",
        timeout=timeout_ms,
    )
    after = _dashboard_smoke_trajectory_state(page)
    return {
        "scatter_points": info.get("n"),
        "glyph_markers_before_anchor": before.get("glyph_markers"),
        "glyph_markers_after_anchor": after.get("glyph_markers"),
        "active_picker_before_anchor": before.get("active_picker_buttons"),
        "active_picker_after_anchor": after.get("active_picker_buttons"),
    }


def _dashboard_smoke_trajectory_state(page) -> dict:
    return page.evaluate(
        """() => {
          var overlay = document.getElementById('traj-glyph-overlay');
          return {
            glyph_markers: overlay
              ? overlay.querySelectorAll('.cryo-traj-glyph-marker').length
              : 0,
            active_picker_buttons: document.querySelectorAll(
              '.cryo-vslice-vol-btn--active'
            ).length,
          };
        }"""
    )


def _dashboard_smoke_volume_viewer_ready(
    page,
    base_url: str,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
) -> dict | None:
    """Navigate to trajectory creator (integrated volume viewer) and wait for manual picker."""
    base = base_url.rstrip("/")
    page.goto(
        f"{base}/trajectory",
        wait_until="domcontentloaded",
        timeout=timeout_ms,
    )
    if page.locator("#vslice-canvas").count() == 0:
        body = page.inner_text("body")
        if "CUDA" in body or "GPU" in body or "weights" in body.lower():
            return None
        raise RuntimeError(
            "volume slice viewer missing #vslice-canvas without known ineligible message"
        )
    scatter_info = _dashboard_smoke_wait_plot_ready(
        page, "scatter", timeout_ms=timeout_ms
    )
    page.wait_for_function(
        """() => {
          var manual = document.getElementById('traj-mode-manual');
          var picker = document.getElementById('vslice-volume-picker-rows');
          if (!manual || !manual.checked) return false;
          if (!picker) return false;
          return picker.querySelectorAll('.cryo-vslice-vol-btn').length > 0;
        }""",
        timeout=timeout_ms,
    )
    picker_count = page.evaluate(
        """() => {
          var picker = document.getElementById('vslice-volume-picker-rows');
          return picker
            ? picker.querySelectorAll('.cryo-vslice-vol-btn').length
            : 0;
        }"""
    )
    return {
        "scatter_points": scatter_info.get("n"),
        "volume_picker_buttons": int(picker_count),
        "canvas_present": page.locator("#vslice-canvas").count() > 0,
    }


def dashboard_smoke_rerender_manual_volumes(
    page,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
) -> None:
    """Trajectory creator: click Render volumes and wait for analyze load to finish."""
    page.wait_for_function(
        """() => {
          var btn = document.getElementById('btn-rerender-volumes');
          return !!(btn && !btn.hidden);
        }""",
        timeout=timeout_ms,
    )
    page.wait_for_function(
        """() => {
          var btn = document.getElementById('btn-rerender-volumes');
          return !!(btn && !btn.disabled);
        }""",
        timeout=timeout_ms,
    )
    page.click("#btn-rerender-volumes")
    page.wait_for_function(
        """() => {
          var btn = document.getElementById('btn-rerender-volumes');
          if (btn && btn.textContent.indexOf('Rendering volumes') >= 0) return false;
          var overlay = document.getElementById('vslice-rendering-overlay');
          if (overlay && !overlay.hidden) return false;
          var preview = document.getElementById('vslice-chimerax-preview');
          if (preview && !preview.hidden && preview.src) return true;
          var cx = document.getElementById('traj-vol-backend-chimerax');
          if (cx && cx.checked && !cx.disabled) return true;
          return !!(btn && btn.disabled
            && (btn.title || '').indexOf('already match') >= 0);
        }""",
        timeout=timeout_ms,
    )


def dashboard_smoke_volume_viewer(
    page,
    base_url: str,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
) -> dict | None:
    """Trajectory creator: manual k-means/PC picker and slice canvas."""
    return _dashboard_smoke_volume_viewer_ready(page, base_url, timeout_ms=timeout_ms)


def dashboard_smoke_volume_viewer_picker_switch(
    page,
    base_url: str,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
) -> dict | None:
    """Trajectory creator: switch active volume via manual picker buttons."""
    ready = _dashboard_smoke_volume_viewer_ready(page, base_url, timeout_ms=timeout_ms)
    if ready is None:
        return None
    if ready["volume_picker_buttons"] < 2:
        return {**ready, "switched": False, "reason": "fewer_than_two_volumes"}

    pick_state = page.evaluate(
        """() => {
          var buttons = document.querySelectorAll('.cryo-vslice-vol-btn');
          if (buttons.length < 2) return { ok: false };
          buttons[0].click();
          return {
            ok: true,
            fromLabel: buttons[0].textContent.trim(),
            targetIdx: 1
          };
        }"""
    )
    if not pick_state.get("ok"):
        return {**ready, "switched": False, "reason": "fewer_than_two_volumes"}
    page.wait_for_function(
        """() => {
          var buttons = document.querySelectorAll('.cryo-vslice-vol-btn');
          return buttons.length >= 2
            && buttons[0].classList.contains('cryo-vslice-vol-btn--active');
        }""",
        timeout=timeout_ms,
    )
    page.evaluate(
        """() => {
          var buttons = document.querySelectorAll('.cryo-vslice-vol-btn');
          if (buttons.length < 2) return;
          if (buttons[0].classList.contains('cryo-vslice-vol-btn--active')) {
            buttons[0].click();
          }
          buttons[1].click();
        }"""
    )
    page.wait_for_function(
        """() => {
          var buttons = document.querySelectorAll('.cryo-vslice-vol-btn');
          if (buttons.length < 2) return false;
          return !buttons[0].classList.contains('cryo-vslice-vol-btn--active')
            && buttons[1].classList.contains('cryo-vslice-vol-btn--active');
        }""",
        timeout=timeout_ms,
    )
    second_label = page.evaluate(
        """() => {
          var buttons = document.querySelectorAll('.cryo-vslice-vol-btn');
          return buttons.length >= 2 ? buttons[1].textContent.trim() : '';
        }"""
    )
    return {
        **ready,
        "switched": True,
        "from_label": pick_state.get("fromLabel", ""),
        "to_label": second_label,
    }


def dashboard_smoke_particle_explorer_panels(
    page,
    base_url: str,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
) -> dict:
    """Verify volume-explorer panel chrome is present (needs mocked eligibility)."""
    base = base_url.rstrip("/")
    page.goto(f"{base}/explorer", wait_until="domcontentloaded", timeout=timeout_ms)
    _dashboard_smoke_wait_plot_ready(page, "scatter", timeout_ms=timeout_ms)
    panel_ids = [
        "cache-panel-toggle",
        "image-grid-menu-toggle",
        "color-selection-panel-toggle",
        "volumes-panel-toggle",
        "btn-expand-cache",
        "montage-cache-size-input",
        "image-cache-progress",
        "cryo-explorer-cache-panel-body",
        "color-selection-panel-body",
        "volumes-panel-body",
        "color-discrete-switches",
        "image-grid-panel-shell",
    ]
    for pid in panel_ids:
        page.wait_for_selector(f"#{pid}", state="attached", timeout=timeout_ms)
    checks = page.evaluate(
        """() => ({
          progressbar: !!document.querySelector('#image-cache-progress[role="progressbar"]'),
          legendPrimitives: !!document.querySelector('script[src*="cryo_cc_legend_primitives"]'),
          legendModule: !!document.querySelector('script[src*="color_covariate_legend"]'),
        })"""
    )
    if not checks.get("progressbar"):
        raise RuntimeError("image-cache-progress missing progressbar role")
    if not checks.get("legendPrimitives") or not checks.get("legendModule"):
        raise RuntimeError(f"explorer legend scripts missing: {checks!r}")
    return {"panels": len(panel_ids), "volumes_panel": True, **checks}


def dashboard_smoke_particle_explorer_no_volumes_panel(
    page,
    base_url: str,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
) -> dict:
    base = base_url.rstrip("/")
    page.goto(f"{base}/explorer", wait_until="domcontentloaded", timeout=timeout_ms)
    _dashboard_smoke_wait_plot_ready(page, "scatter", timeout_ms=timeout_ms)
    has_volumes = page.locator("#volumes-panel-toggle").count() > 0
    if has_volumes:
        raise RuntimeError("volumes panel present on ineligible explorer page")
    return {"volumes_panel": False}


def dashboard_smoke_command_builder(
    page,
    base_url: str,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
) -> dict:
    base = base_url.rstrip("/")
    page.goto(
        f"{base}/command-builder", wait_until="domcontentloaded", timeout=timeout_ms
    )
    # ``#cmd-out`` lives in a max-height dock and is often clipped (not "visible").
    page.wait_for_selector("#cmd-form", state="attached", timeout=timeout_ms)
    page.wait_for_selector("#cmd-out", state="attached", timeout=timeout_ms)
    out0 = page.evaluate(
        """() => {
          var el = document.getElementById('cmd-out');
          return el ? (el.textContent || '') : '';
        }"""
    )
    if "cryodrgn" not in out0:
        raise RuntimeError(f"cmd-out missing cryodrgn prefix: {out0[:80]!r}")
    _dashboard_smoke_set_select_value(page, "cmd-type", "train_vae")
    page.wait_for_function(
        """() => {
          var out = document.getElementById('cmd-out');
          return out && /train_vae/.test(out.textContent || '');
        }""",
        timeout=timeout_ms,
    )
    out1 = page.evaluate(
        """() => {
          var el = document.getElementById('cmd-out');
          return el ? (el.textContent || '') : '';
        }"""
    )
    return {
        "initial_has_abinit": "abinit" in out0,
        "switched_to_train_vae": "train_vae" in out1,
    }


def dashboard_smoke_particle_explorer_color_covariate(
    page,
    base_url: str,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
) -> dict:
    base = base_url.rstrip("/")
    page.goto(f"{base}/explorer", wait_until="domcontentloaded", timeout=timeout_ms)
    info = _dashboard_smoke_wait_plot_ready(page, "scatter", timeout_ms=timeout_ms)
    page.wait_for_function(
        """() => {
          var el = document.getElementById('sc');
          return el && el.options && el.options.length > 1;
        }""",
        timeout=timeout_ms,
    )
    _dashboard_smoke_set_select_value(page, "sc", "labels")
    page.wait_for_function(
        """() => {
          var t = document.getElementById('color-selection-panel-toggle');
          return t && t.getAttribute('aria-disabled') === 'false';
        }""",
        timeout=timeout_ms,
    )
    _dashboard_smoke_ensure_panel_open(page, "color-selection-panel-toggle")
    page.wait_for_function(
        """() => {
          var switches = document.getElementById('color-discrete-switches');
          return switches && switches.querySelectorAll('button, label, input').length > 0;
        }""",
        timeout=timeout_ms,
    )
    toggles = page.locator(
        "#color-discrete-switches button, #color-discrete-switches label"
    ).count()
    plot = _dashboard_smoke_eval_plot(page, "scatter")
    return {
        "scatter_points": info.get("n"),
        "discrete_toggles": toggles,
        "points_after_color": plot.get("points"),
    }


def dashboard_smoke_particle_explorer_cache_expand(
    page,
    base_url: str,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
    initial_size: int = 10,
    expanded_size: int = 20,
) -> dict:
    base = base_url.rstrip("/")
    page.goto(f"{base}/explorer", wait_until="domcontentloaded", timeout=timeout_ms)
    _dashboard_smoke_wait_plot_ready(page, "scatter", timeout_ms=timeout_ms)
    _dashboard_smoke_clear_explorer_cache_if_any(page, timeout_ms=timeout_ms)
    _dashboard_smoke_ensure_panel_open(page, "cache-panel-toggle")
    page.fill("#montage-cache-size-input", str(initial_size))
    expand = page.locator("#btn-expand-cache")
    expand.wait_for(state="visible", timeout=timeout_ms)
    if expand.is_disabled():
        raise RuntimeError(
            f"expand cache disabled before initial build (size={initial_size})"
        )
    expand.click()
    _dashboard_smoke_preload_overlay_hidden(page, timeout_ms=timeout_ms)
    page.wait_for_function(
        """(minN) => {
          var st = document.getElementById('preload-status');
          if (!st) return false;
          var m = (st.textContent || '').match(/([\\d,]+)\\s+images\\s+cached/i);
          return m && parseInt(m[1].replace(/,/g, ''), 10) >= minN;
        }""",
        arg=initial_size,
        timeout=timeout_ms,
    )
    have0 = page.evaluate(
        """() => {
          var st = document.getElementById('preload-status');
          var m = (st && st.textContent || '').match(/([\\d,]+)\\s+images\\s+cached/i);
          return m ? parseInt(m[1].replace(/,/g, ''), 10) : 0;
        }"""
    )
    page.fill("#montage-cache-size-input", str(expanded_size))
    expand.click()
    _dashboard_smoke_preload_overlay_hidden(page, timeout_ms=timeout_ms)
    page.wait_for_function(
        """(minN) => {
          var st = document.getElementById('preload-status');
          if (!st) return false;
          var m = (st.textContent || '').match(/([\\d,]+)\\s+images\\s+cached/i);
          return m && parseInt(m[1].replace(/,/g, ''), 10) >= minN;
        }""",
        arg=expanded_size,
        timeout=timeout_ms,
    )
    have1 = page.evaluate(
        """() => {
          var st = document.getElementById('preload-status');
          var m = (st && st.textContent || '').match(/([\\d,]+)\\s+images\\s+cached/i);
          return m ? parseInt(m[1].replace(/,/g, ''), 10) : 0;
        }"""
    )
    if have1 <= have0:
        raise RuntimeError(f"cache did not expand: {have0} -> {have1}")
    return {"initial_cached": have0, "expanded_cached": have1}


def dashboard_smoke_landscape_volpca_clear_selection(
    page,
    base_url: str,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
) -> dict:
    base = base_url.rstrip("/")
    page.goto(
        f"{base}/landscape-volpca", wait_until="domcontentloaded", timeout=timeout_ms
    )
    _dashboard_smoke_wait_plot_ready(page, "volsketch", timeout_ms=timeout_ms)
    page.click("#volsketch-random-sel")
    page.wait_for_function(
        """() => {
          var gd = document.getElementById('volsketch');
          if (!gd || !gd.data) return false;
          for (var i = 0; i < gd.data.length; i++) {
            if (gd.data[i] && gd.data[i].name === 'cdrgnVolSelectionOverlay') return true;
          }
          return false;
        }""",
        timeout=30_000,
    )
    page.click("#volsketch-clear-sel")
    page.wait_for_function(
        """() => {
          var summary = document.getElementById('volsketch-sel-summary');
          var txt = summary ? summary.textContent || '' : '';
          return /0\\/\\d+ sketched|No sketched volumes selected/i.test(txt);
        }""",
        timeout=timeout_ms,
    )
    return {"overlay_cleared": True}


def dashboard_smoke_landscape_volpca_axis_reload(
    page,
    base_url: str,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
) -> dict:
    base = base_url.rstrip("/")
    page.goto(
        f"{base}/landscape-volpca", wait_until="domcontentloaded", timeout=timeout_ms
    )
    info = _dashboard_smoke_wait_plot_ready(page, "volsketch", timeout_ms=timeout_ms)
    before = _dashboard_smoke_eval_plot(page, "volsketch")
    options = page.evaluate(
        """() => Array.from(document.querySelectorAll('#volsketch-pcx option'))
          .map(function(o) { return o.value; })
          .filter(Boolean)"""
    )
    if len(options) < 2:
        raise RuntimeError(f"need >=2 vol PCA axis options, got {options!r}")
    alt = options[1] if options[0] == page.input_value("#volsketch-pcx") else options[0]
    page.select_option("#volsketch-pcx", alt)
    page.wait_for_function(
        """(prev) => {
          var gd = document.getElementById('volsketch');
          var P = window.CryoPlotlyArrays;
          if (!gd || !P || !gd.data || !gd.data[0]) return false;
          var pts = P.length(gd.data[0].x);
          return pts > 0 && pts === prev;
        }""",
        arg=before.get("points"),
        timeout=timeout_ms,
    )
    after = _dashboard_smoke_eval_plot(page, "volsketch")
    return {"scatter_points": info.get("n"), "axis": alt, **after}


def dashboard_smoke_landscape_full_3d_clear_selection(
    page,
    base_url: str,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
) -> dict | None:
    base = base_url.rstrip("/")
    page.goto(
        f"{base}/landscape-full-3d", wait_until="domcontentloaded", timeout=timeout_ms
    )
    if page.locator("#latent3d").count() == 0:
        body = page.inner_text("body")
        if "analyze_landscape_full" in body or "need" in body.lower():
            return None
        raise RuntimeError("unexpected landscape-full-3d page without #latent3d")
    _dashboard_smoke_wait_plot_ready(page, "latent3d", timeout_ms=timeout_ms)
    page.click("#l3dva-random-sel")
    page.wait_for_function(
        """() => {
          var gd = document.getElementById('latent3d');
          var anns = gd && gd.layout && gd.layout.scene && gd.layout.scene.annotations;
          return anns && anns.length > 0;
        }""",
        timeout=timeout_ms,
    )
    page.click("#l3dva-clear-sel")
    page.wait_for_function(
        """() => {
          var gd = document.getElementById('latent3d');
          var anns = gd && gd.layout && gd.layout.scene && gd.layout.scene.annotations;
          return !anns || anns.length === 0;
        }""",
        timeout=timeout_ms,
    )
    return {"annotations_cleared": True}


def run_dashboard_plotly_smoke(
    base_url: str,
    *,
    timeout_ms: int = DASHBOARD_BROWSER_SMOKE_TIMEOUT_MS,
    cache_size: int | None = None,
) -> int:
    """Run all dashboard browser smoke steps (for ``scripts/dashboard_plotly_smoke.py``)."""
    import json
    import sys

    pytest.importorskip("playwright")
    from playwright.sync_api import Error as PlaywrightError
    from playwright.sync_api import TimeoutError as PlaywrightTimeout
    from playwright.sync_api import sync_playwright

    results: dict[str, object] = {}
    errors: list[str] = []
    skipped: list[str] = []

    with sync_playwright() as p:
        try:
            browser = p.chromium.launch(**_playwright_chromium_launch_kwargs())
        except PlaywrightError as exc:
            print(
                f"Playwright Chromium not installed ({exc}). "
                "Run: playwright install chromium",
                file=sys.stderr,
            )
            return 1
        context = browser.new_context(viewport={"width": 1400, "height": 900})
        page = context.new_page()

        steps = [
            (
                "index",
                lambda: dashboard_smoke_index(page, base_url, timeout_ms=timeout_ms),
            ),
            (
                "command_builder",
                lambda: dashboard_smoke_command_builder(
                    page, base_url, timeout_ms=timeout_ms
                ),
            ),
            (
                "latent_3d",
                lambda: dashboard_smoke_latent_3d(
                    page, base_url, timeout_ms=timeout_ms
                ),
            ),
            (
                "pairplot",
                lambda: dashboard_smoke_pairplot(page, base_url, timeout_ms=timeout_ms),
            ),
            (
                "trajectory",
                lambda: dashboard_smoke_trajectory(
                    page, base_url, timeout_ms=timeout_ms
                ),
            ),
            (
                "particle_explorer",
                lambda: dashboard_smoke_particle_explorer(
                    page, base_url, timeout_ms=timeout_ms, cache_size=cache_size
                ),
            ),
            (
                "landscape_volpca",
                lambda: dashboard_smoke_landscape_volpca(
                    page, base_url, timeout_ms=timeout_ms
                ),
            ),
        ]
        page.goto(
            f"{base_url.rstrip('/')}/landscape-full-3d",
            wait_until="domcontentloaded",
            timeout=60_000,
        )
        if page.locator("#latent3d").count():
            steps.append(
                (
                    "landscape_full_3d",
                    lambda: dashboard_smoke_landscape_full_3d(
                        page, base_url, timeout_ms=timeout_ms
                    ),
                )
            )
        else:
            skipped.append("landscape_full_3d (no analyze_landscape_full outputs)")

        skip_labels = {
            "trajectory": "trajectory (no CUDA GPU / weights)",
            "landscape_volpca": "landscape_volpca (no analyze_landscape outputs)",
            "landscape_full_3d": "landscape_full_3d (no analyze_landscape_full outputs)",
        }

        for name, fn in steps:
            try:
                result = fn()
                if result is None:
                    skipped.append(skip_labels.get(name, f"{name} (skipped)"))
                    continue
                results[name] = result
                print(f"PASS {name}: {json.dumps(results[name], sort_keys=True)}")
            except (RuntimeError, PlaywrightTimeout) as exc:
                errors.append(f"{name}: {exc}")
                print(f"FAIL {name}: {exc}", file=sys.stderr)
            page.wait_for_timeout(500)

        browser.close()

    for s in skipped:
        print(f"SKIP {s}")
    if errors:
        print("\nSmoke test FAILED:", "; ".join(errors), file=sys.stderr)
        return 1
    print("\nSmoke test PASSED (dashboard browser smoke)")
    return 0
