"""The setup command used to help create output directories for reconstruction."""

import pytest
import os.path
import shutil
import argparse
from cryodrgn.commands import setup, train
import cryodrgn.utils
from typing import Optional, Any


class SetupRequest:
    """A set of parameters passed to the `cryodrgn setup `command."""

    def __init__(
        self,
        outdir: str,
        model: Optional[str],
        dataset: Optional[str],
        particles: Optional[str],
        ctf: Optional[str],
        poses: Optional[str],
        ind: Optional[str],
        datadir: Optional[str],
        reconstruction_type: Optional[str],
        z_dim: Optional[str],
        pose_estimation: Optional[str],
        tilt: Optional[bool],
        cfgs: Optional[list[str]],
        include_cfgs: Optional[dict[str, Any]],
    ):
        args = [outdir]
        if model is not None:
            args += ["--model", model]
        if dataset is not None:
            args += ["--dataset", dataset]
        if particles is not None:
            args += ["--particles", particles]
        if ctf is not None:
            args += ["--ctf", ctf]
        if poses is not None:
            args += ["--poses", poses]
        if ind is not None:
            args += ["--ind", ind]
        if datadir is not None:
            args += ["--datadir", datadir]
        if reconstruction_type is not None:
            args += ["--reconstruction-type", reconstruction_type]
        if z_dim is not None:
            args += ["--z-dim", z_dim]
        if pose_estimation is not None:
            args += ["--pose-estimation", pose_estimation]
        if tilt is not None and tilt is True:
            args += ["--tilt"]

        add_cfgs = [
            "pretrain=10",
            "hidden_dim=8",
            "checkpoint=1",
            "learning_rate=0.1",
            "t_ngrid=2",
            "pe_dim=4",
        ]
        if model == "cryodrgn":
            add_cfgs += ["ps_freq=2"]

        if include_cfgs is not None:
            os.makedirs(outdir, exist_ok=True)
            include_fl = os.path.join(outdir, "test-include.yml")
            cryodrgn.utils.save_yaml(include_cfgs, include_fl)
            args += ["--include", include_fl]

        if cfgs:
            add_cfgs += cfgs
        args += ["--cfgs"] + add_cfgs

        self.args = args
        self.outdir = outdir
        self.model = model
        self.particles = particles
        self.ctf = ctf
        self.poses = poses
        self.ind = ind
        self.datadir = datadir
        self.reconstruction_type = reconstruction_type
        self.z_dim = z_dim
        self.pose_estimation = pose_estimation
        self.tilt = tilt
        self.cfgs = cfgs if cfgs is not None else list()
        self.include_cfgs = include_cfgs if include_cfgs is not None else dict()


@pytest.mark.parametrize(
    "model, particles, ctf, poses, indices, datadir, dataset, pose_estimation",
    [
        ("cryodrgn-ai", "toy.mrcs", "CTF-Test", None, "5", None, None, None),
        ("cryodrgn-ai", "toy.mrcs", "CTF-Test", "toy-poses", None, None, None, None),
        pytest.param(
            "cryodrgn-ai",
            "toy.mrcs",
            "CTF-Test",
            None,
            None,
            None,
            None,
            "fixed",
            marks=pytest.mark.xfail(
                raises=(ValueError, FileNotFoundError),
                reason="fixed estimation but no poses given",
            ),
        ),
        ("cryodrgn", "toy.mrcs", None, None, "5", None, None, "abinit"),
        ("cryodrgn-ai", "toy.txt", "CTF-Test", "toy-poses", None, None, None, "fixed"),
        ("cryodrgn", "hand", None, "hand-poses", "5", None, None, "fixed"),
        ("cryodrgn", "toy.txt", "CTF-Test", "toy-poses", "5", None, None, "abinit"),
    ],
    indirect=["particles", "ctf", "poses", "indices", "datadir"],
)
class TestSetupThenRun:
    @pytest.fixture
    def setup_request(
        self,
        tmpdir_factory,
        model,
        dataset,
        particles,
        ctf,
        poses,
        indices,
        datadir,
        reconstruction_type,
        z_dim,
        pose_estimation,
        cfgs,
        include,
    ):

        dirname = os.path.join(
            "SetupThenRun",
            model if model is not None else "None",
            dataset if dataset is not None else "None",
            particles.label,
            poses.label,
            ctf.label,
            poses.label,
            indices.label,
            datadir.label,
            reconstruction_type if reconstruction_type is not None else "None",
            z_dim if z_dim is not None else "None",
            pose_estimation if pose_estimation is not None else "None",
            str(hash(tuple(cfgs)))[1:10] if cfgs is not None else "None",
            str(hash(tuple(include))) if include is not None else "None",
        )
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)

        return SetupRequest(
            odir,
            model,
            dataset,
            particles.path,
            ctf.path,
            poses.path,
            indices.path,
            datadir.path,
            reconstruction_type,
            z_dim,
            pose_estimation,
            False,
            cfgs,
            include,
        )

    @pytest.mark.parametrize(
        "reconstruction_type, z_dim, cfgs, include",
        [
            (None, None, None, None),
            (None, None, None, {"weight_decay": 0.05}),
            (None, None, ["z_dim=2"], None),
            (None, "0", ["window_r=0.80"], None),
            (None, "4", ["window_r=0.80", "z_dim=2"], None),
            (None, "4", ["window_r=0.80", "z_dim=2"], {"weight_decay": 0.05}),
            ("homo", None, ["window_r=0.80"], None),
            ("homo", None, ["window_r=0.80"], {"weight_decay": 0.05}),
            ("het", None, ["window_r=0.75", "z_dim=2"], None),
            pytest.param(
                "homo",
                "4",
                None,
                None,
                marks=pytest.mark.xfail(
                    raises=ValueError,
                    reason="cannot specify both reconstruction-type and z_dim",
                ),
            ),
            pytest.param(
                "het",
                "12",
                None,
                {"weight_decay": 0.05},
                marks=pytest.mark.xfail(
                    raises=ValueError,
                    reason="cannot specify both reconstruction-type and z_dim",
                ),
            ),
        ],
    )
    def test_setup_directory(self, setup_request):
        """Create a reconstruction directory using the setup command."""

        parser = argparse.ArgumentParser()
        setup.add_args(parser)
        setup.main(parser.parse_args(setup_request.args))
        assert os.path.isdir(setup_request.outdir), "Missing output directory!"
        assert os.path.isfile(
            os.path.join(setup_request.outdir, "config.yaml")
        ), "Missing config file!"
        configs = cryodrgn.utils.load_yaml(
            os.path.join(setup_request.outdir, "config.yaml")
        )

        assert configs["dataset_args"]["particles"] == setup_request.particles
        assert configs["dataset_args"]["ctf"] == setup_request.ctf
        assert configs["dataset_args"]["poses"] == setup_request.poses
        assert configs["dataset_args"]["ind"] == setup_request.ind
        assert configs["dataset_args"]["datadir"] == setup_request.datadir
        assert configs["model_args"]["model"] == setup_request.model
        if setup_request.z_dim is not None:
            use_zdim = int(setup_request.z_dim)
        elif setup_request.reconstruction_type is None:
            use_zdim = 0
        elif setup_request.reconstruction_type == "homo":
            use_zdim = 0
        elif setup_request.reconstruction_type == "het":
            use_zdim = 8
        else:
            raise ValueError(
                f"Unrecognized reconstruction type "
                f"`{setup_request.reconstruction_type}`!"
            )
        for cfg in setup_request.cfgs:
            if cfg.startswith("z_dim="):
                use_zdim = int(cfg.split("=")[1])

            elif cfg.startswith("window_r="):
                assert configs["dataset_args"]["window_r"] == float(
                    cfg.split("=")[1]
                ), "Incorrect window radius!"

        assert configs["model_args"]["z_dim"] == use_zdim
        assert configs["model_args"]["learning_rate"] == 0.1
        if setup_request.pose_estimation is None:
            if setup_request.poses is not None:
                assert configs["model_args"]["pose_estimation"] == "fixed"
            else:
                assert configs["model_args"]["pose_estimation"] == "abinit"
        else:
            assert (
                configs["model_args"]["pose_estimation"]
                == setup_request.pose_estimation
            )

        for par_k, par_v in setup_request.include_cfgs.items():
            assert configs["model_args"][par_k] == par_v

    @pytest.mark.parametrize(
        "reconstruction_type, z_dim, cfgs, include",
        [
            (None, None, None, None),
            (None, None, None, {"weight_decay": 0.05}),
            (None, "0", ["window_r=0.80"], None),
            (None, "4", ["window_r=0.80", "z_dim=2"], {"weight_decay": 0.05}),
            ("homo", None, ["window_r=0.80"], None),
        ],
    )
    def test_then_train(self, setup_request):
        """Run the reconstruction experiment using the directory created by setup."""

        args = [setup_request.outdir, "--num-epochs", "3", "--no-analysis"]
        parser = argparse.ArgumentParser()
        train.add_args(parser)
        train.main(parser.parse_args(args))

        out_files = os.listdir(setup_request.outdir)
        for epoch in range(1, 10):
            if epoch <= 3:
                assert (
                    f"weights.{epoch}.pkl" in out_files
                ), f"Missing output model weights for epoch {epoch}!"
                if (
                    setup_request.pose_estimation == "abinit"
                    or setup_request.poses is None
                ):
                    assert (
                        f"pose.{epoch}.pkl" in out_files
                    ), f"Missing output model poses for epoch {epoch}!"
                else:
                    assert (
                        f"pose.{epoch}.pkl" not in out_files
                    ), f"Extra output model poses for epoch {epoch}!"
            else:
                assert (
                    f"weights.{epoch}.pkl" not in out_files
                ), f"Extra output model weights for epoch {epoch}!"

    @pytest.mark.parametrize(
        "reconstruction_type, z_dim, cfgs, include",
        [
            (None, None, None, None),
            (None, "4", ["window_r=0.80", "z_dim=2"], {"weight_decay": 0.05}),
        ],
    )
    def test_then_rerun(self, setup_request):
        """Copy the experiment config file to an empty directory and run again."""

        rerun_dir = os.path.join(setup_request.outdir, "rerun")
        os.makedirs(rerun_dir)
        shutil.copyfile(
            os.path.join(setup_request.outdir, "config.yaml"),
            os.path.join(rerun_dir, "config.yaml"),
        )
        args = [rerun_dir, "--num-epochs", "1", "--no-analysis"]
        parser = argparse.ArgumentParser()
        train.add_args(parser)
        train.main(parser.parse_args(args))

        out_files = os.listdir(rerun_dir)
        for epoch in range(1, 5):
            if epoch <= 1:
                assert (
                    f"weights.{epoch}.pkl" in out_files
                ), f"Missing output model weights for epoch {epoch}!"
                if (
                    setup_request.pose_estimation == "abinit"
                    or setup_request.poses is None
                ):
                    assert (
                        f"pose.{epoch}.pkl" in out_files
                    ), f"Missing output model poses for epoch {epoch}!"
            else:
                assert (
                    f"weights.{epoch}.pkl" not in out_files
                ), f"Extra output model weights for epoch {epoch}!"
