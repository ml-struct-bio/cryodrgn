"""Integration tests of ab initio reconstruction and downstream analyses.

The training done here has unrealistically low parameter values to allow the tests to
run in any environment in a reasonable amount of time with or without GPUs.

"""
import pytest
import os
import shutil
from typing import Optional, Any
from cryodrgn.utils import run_command

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "testing", "data")


class AbInitioDir:
    def __init__(
        self,
        zdim: int,
        dataset: str = "hand",
        epochs: int = 2,
        seed: Optional[int] = None,
    ) -> None:
        self.zdim = zdim
        self.dataset = dataset
        self.epochs = epochs
        self.seed = seed

        self.outdir = os.path.abspath(f"test-output_{dataset}")
        os.makedirs(self.outdir)

        if self.dataset == "toy":
            self.particles = os.path.join(DATA_DIR, "toy_projections.mrcs")
        elif self.dataset == "hand":
            self.particles = os.path.join(DATA_DIR, "hand.mrcs")
        else:
            raise ValueError(f"Unrecognized dataset label `{self.dataset}`!")

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


@pytest.mark.parametrize(
    "abinit_dir", [{"zdim": zdim} for zdim in [0, 4, 8]], indirect=True
)
def test_analysis_and_backproject(abinit_dir):
    abinit_dir.train()
    abinit_dir.train(load_epoch=0)
    abinit_dir.backproject()
    abinit_dir.view_config()
