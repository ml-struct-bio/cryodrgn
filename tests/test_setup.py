"""Volume reconstruction training with fixed poses, as well as downstream analyses."""

import pytest
import os.path
import argparse
from cryodrgn.commands import setup, train
import cryodrgn.utils


class SetupRequest:
    def __init__(
        self,
        outdir,
        model,
        dataset,
        particles,
        ctf,
        poses,
        ind,
        datadir,
        reconstruction_type,
        z_dim,
        pose_estimation,
        tilt,
        cfgs,
        defaults,
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

        if cfgs:
            add_cfgs += cfgs
        args += ["--cfgs"] + add_cfgs
        if defaults:
            args += ["--defaults"]

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
        self.cfgs = cfgs
        self.defaults = defaults


@pytest.mark.parametrize(
    "model, ctf",
    [("cryodrgn", "CTF-Test"), ("cryodrgn", None), ("cryodrgn-ai", "CTF-Test")],
    indirect=["ctf"],
)
@pytest.mark.parametrize("dataset", [None])
@pytest.mark.parametrize(
    "particles, datadir", [("toy.mrcs", None), ("toy.txt", None)], indirect=True
)
@pytest.mark.parametrize(
    "pose_estimation, poses, indices",
    [
        (None, None, "5"),
        (None, "toy-poses", None),
        pytest.param(
            "fixed",
            None,
            None,
            marks=pytest.mark.xfail(
                raises=(ValueError, FileNotFoundError),
                reason="fixed estimation but no poses given",
            ),
        ),
        ("abinit", None, "5"),
        ("fixed", "toy-poses", None),
        ("abinit", "toy-poses", "5"),
    ],
    indirect=["poses", "indices"],
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
        )
        odir = os.path.join(tmpdir_factory.getbasetemp(), dirname)
        os.makedirs(odir, exist_ok=True)

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
            None,
            None,
        )

    @pytest.mark.parametrize(
        "reconstruction_type, z_dim",
        [
            (None, None),
            (None, "0"),
            (None, "4"),
            ("homo", None),
            ("het", None),
            pytest.param(
                "homo",
                "4",
                marks=pytest.mark.xfail(
                    raises=ValueError,
                    reason="cannot specify both reconstruction-type and z_dim",
                ),
            ),
        ],
    )
    def test_use_setup(self, setup_request):
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
            assert configs["model_args"]["z_dim"] == int(setup_request.z_dim)
        elif setup_request.reconstruction_type is None:
            assert configs["model_args"]["z_dim"] == 0
        elif setup_request.reconstruction_type == "homo":
            assert configs["model_args"]["z_dim"] == 0
        elif setup_request.reconstruction_type == "het":
            assert configs["model_args"]["z_dim"] == 8
        else:
            assert False

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

    @pytest.mark.parametrize(
        "reconstruction_type, z_dim",
        [(None, "0"), (None, "4"), ("homo", None)],
    )
    def test_then_run(self, setup_request):
        args = [setup_request.outdir, "--num-epochs", "2", "--no-analysis"]
        parser = argparse.ArgumentParser()
        train.add_args(parser)
        train.main(parser.parse_args(args))

        out_files = os.listdir(setup_request.outdir)
        assert "weights.2.pkl" in out_files, "Missing output model weights!"
