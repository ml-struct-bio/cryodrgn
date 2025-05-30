"""The setup command used to help create output directories for reconstruction."""

import pytest
import os.path
import shutil
import argparse
from pathlib import Path
from typing import Optional, Any
from cryodrgn.commands import setup, train
import cryodrgn.utils
from cryodrgn.models.utils import get_model_configurations
from cryodrgn.trainers.sgd_trainer import SGDPoseSearchConfigurations


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
        zdim: Optional[str],
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
        if zdim is not None:
            args += ["--zdim", zdim]
        if pose_estimation is not None:
            args += ["--pose-estimation", pose_estimation]
        if tilt is not None and tilt is True:
            args += ["--tilt"]

        add_cfgs = [
            "pretrain=10",
            "hidden_dim=8",
            "learning_rate=0.1",
            "t_ngrid=2",
            "pe_dim=4",
        ]
        # TODO: find a way to not hard-code these
        if model == "autoenc":
            add_cfgs += ["ps_freq=2"]
        else:
            add_cfgs += ["n_imgs_pose_search=10"]

        if include_cfgs is not None:
            os.makedirs(outdir, exist_ok=True)
            include_fl = os.path.normpath(
                os.path.join(outdir, "..", "test-include.yml")
            )
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
        self.zdim = zdim
        self.pose_estimation = pose_estimation
        self.tilt = tilt
        self.cfgs = cfgs if cfgs is not None else list()
        self.include_cfgs = include_cfgs if include_cfgs is not None else dict()


@pytest.fixture
def setup_request(
    tmpdir_factory,
    model,
    dataset,
    particles,
    ctf,
    poses,
    indices,
    datadir,
    reconstruction_type,
    zdim,
    pose_estimation,
    cfgs,
    include,
) -> SetupRequest:
    dirname = Path(
        "SetupThenRun",
        "_".join(
            [
                model if model is not None else "None",
                dataset if dataset is not None else "None",
                particles.label,
                poses.label,
                ctf.label,
                poses.label,
                indices.label,
                datadir.label,
                reconstruction_type if reconstruction_type is not None else "None",
                zdim if zdim is not None else "None",
                pose_estimation if pose_estimation is not None else "None",
                str(hash(tuple(cfgs))) if cfgs is not None else "None",
                str(hash(tuple(include))) if include is not None else "None",
            ]
        ),
    )
    odir = Path(tmpdir_factory.getbasetemp(), dirname)

    return SetupRequest(
        str(odir),
        model,
        dataset,
        particles.path,
        ctf.path,
        poses.path,
        indices.path,
        datadir.path,
        reconstruction_type,
        zdim,
        pose_estimation,
        False,
        cfgs,
        include,
    )


def test_empty_setup(tmpdir_factory):
    parser = argparse.ArgumentParser()
    setup.add_args(parser)
    outdir = tmpdir_factory.mktemp("empty_setup").strpath
    setup.main(parser.parse_args([outdir]))
    assert os.path.isfile(os.path.join(outdir, "config.yaml"))

    cfgs = cryodrgn.utils.load_yaml(os.path.join(outdir, "config.yaml"))
    configs = get_model_configurations(cfgs)
    assert configs == SGDPoseSearchConfigurations(zdim=8, seed=configs.seed)


@pytest.mark.xfail(
    reason="Output directory already exists and is not empty!",
    raises=ValueError,
)
def test_nonempty_setup(tmpdir_factory):
    parser = argparse.ArgumentParser()
    setup.add_args(parser)
    outdir = tmpdir_factory.mktemp("nonempty_setup").strpath
    setup.main(parser.parse_args([outdir]))
    assert os.path.isfile(os.path.join(outdir, "config.yaml"))
    setup.main(parser.parse_args([outdir]))


def setup_directory(setup_request):
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

    # Default model if not given should always be autodec
    if setup_request.model is not None:
        assert configs["model_args"]["model"] == setup_request.model
    else:
        assert configs["model_args"]["model"] == "autodec"

    if setup_request.particles is None:
        assert configs["dataset_args"]["particles"] is None
    else:
        assert configs["dataset_args"]["particles"] == os.path.abspath(
            setup_request.particles
        )
    if setup_request.ctf is None:
        assert configs["dataset_args"]["ctf"] is None
    else:
        assert configs["dataset_args"]["ctf"] == os.path.abspath(setup_request.ctf)
    if setup_request.poses is None:
        assert configs["dataset_args"]["poses"] is None
    else:
        assert configs["dataset_args"]["poses"] == os.path.abspath(setup_request.poses)

    for cfg in setup_request.cfgs:
        if cfg.startswith("ind="):
            assert configs["dataset_args"]["ind"] == int(cfg.split("=")[1])
            break
    else:
        if setup_request.ind is None:
            assert configs["dataset_args"]["ind"] is None
        elif not setup_request.ind.isdigit():
            assert configs["dataset_args"]["ind"] == os.path.abspath(setup_request.ind)

    if setup_request.zdim is not None:
        use_zdim = int(setup_request.zdim)
    elif setup_request.reconstruction_type is None:
        use_zdim = 8
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
        if cfg.startswith("zdim="):
            use_zdim = int(cfg.split("=")[1])

        elif cfg.startswith("window_r="):
            assert configs["dataset_args"]["window_r"] == float(
                cfg.split("=")[1]
            ), "Incorrect window radius!"

    assert configs["model_args"]["zdim"] == use_zdim
    assert configs["model_args"]["learning_rate"] == 0.1
    if setup_request.pose_estimation is None:
        if setup_request.poses is not None:
            assert configs["model_args"]["pose_estimation"] == "fixed"
        else:
            assert configs["model_args"]["pose_estimation"] == "abinit"
    else:
        assert configs["model_args"]["pose_estimation"] == setup_request.pose_estimation

    for par_k, par_v in setup_request.include_cfgs.items():
        assert configs["model_args"][par_k] == par_v


@pytest.mark.parametrize(
    "model, particles, ctf, poses, indices, datadir, dataset, pose_estimation",
    [
        ("autoenc", "toy.mrcs", "CTF-Test", None, "5", None, None, None),
        ("autodec", "toy.mrcs", "CTF-Test", None, "5", None, None, None),
        ("autodec", "toy.mrcs", "CTF-Test", "toy-poses", None, None, None, None),
        pytest.param(
            "autodec",
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
        ("autoenc", "toy.mrcs", None, None, "5", None, None, "abinit"),
        ("autodec", "toy.txt", "CTF-Test", "toy-poses", None, None, None, "fixed"),
        ("autoenc", "hand", None, "hand-poses", "5", None, None, "fixed"),
        ("autoenc", "toy.txt", "CTF-Test", "toy-poses", "5", None, None, "abinit"),
        (None, "toy.txt", "CTF-Test", None, "5", None, None, "abinit"),
    ],
    indirect=["particles", "ctf", "poses", "indices", "datadir"],
)
@pytest.mark.parametrize(
    "reconstruction_type, zdim, cfgs, include",
    [
        (None, None, None, None),
        ("homo", None, None, None),
        (None, "2", None, None),
        ("homo", None, None, {"weight_decay": 0.05}),
        (None, None, ["zdim=2"], None),
        (None, "0", ["window_r=0.80"], None),
        (None, "4", ["window_r=0.80", "zdim=2"], None),
        (None, "4", ["window_r=0.80", "zdim=2", "ind=4"], {"weight_decay": 0.05}),
        ("homo", None, ["window_r=0.80"], None),
        ("homo", None, ["window_r=0.80"], {"weight_decay": 0.05}),
        ("het", None, ["window_r=0.75", "zdim=2"], None),
        pytest.param(
            "homo",
            "4",
            None,
            None,
            marks=pytest.mark.xfail(
                raises=ValueError,
                reason="cannot specify both reconstruction-type and zdim",
            ),
        ),
        pytest.param(
            "het",
            "12",
            None,
            {"weight_decay": 0.05},
            marks=pytest.mark.xfail(
                raises=ValueError,
                reason="cannot specify both reconstruction-type and zdim",
            ),
        ),
    ],
)
def test_setup_standalone(setup_request):
    setup_directory(setup_request)
    shutil.rmtree(setup_request.outdir)


@pytest.mark.parametrize(
    "model, particles, ctf, poses, indices, datadir, dataset, pose_estimation",
    [
        (
            "autoenc",
            "data/toy_projections.mrcs",
            "CTF-Test",
            None,
            "just-5",
            None,
            None,
            None,
        ),
        ("autodec", "toy.mrcs", "data/test_ctf.pkl", None, "5", None, None, None),
        ("autodec", "toy.mrcs", "CTF-Test", "toy-poses", None, None, None, None),
        ("autoenc", "toy.mrcs", None, None, "5", None, None, "abinit"),
        ("autodec", "toy.txt", "CTF-Test", "toy-poses", None, None, None, "fixed"),
        ("autoenc", "hand", None, "hand-poses", "just-5", None, None, "fixed"),
        (
            "autoenc",
            "data/toy_projections.txt",
            "CTF-Test",
            "data/toy_rot_trans.pkl",
            "5",
            None,
            None,
            "abinit",
        ),
        (None, "toy.txt", "CTF-Test", None, "just-5", None, None, "abinit"),
    ],
    indirect=["particles", "ctf", "poses", "indices", "datadir"],
    ids=[
        "v3,relative_particles,absolute_ctf,indices_file,no_poses,implicit_abinit",
        "v4,absolute_particles,relative_ctf,indices_int,no_poses,implicit_abinit",
        "v4,absolute_particles,absolute_ctf,no_indices,with_poses,implicit_abinit",
        "v3,absolute_particles,no_ctf,indices_int,no_poses,explicit_abinit",
        "v4,abstxt_particles,absolute_ctf,no_indices,with_poses,explicit_fixed",
        "v3,hand_particles,no_ctf,indices_file,with_poses,explicit_fixed",
        "v4,abstxt_particles,absolute_ctf,indices_int,relative_poses,explicit_abinit",
        "implicit_v4,absolute_particles,absolute_ctf,indices_file,with_poses,explicit_abinit",
    ],
)
@pytest.mark.parametrize(
    "reconstruction_type, zdim, cfgs, include",
    [
        ("homo", None, None, {"weight_decay": 0.05}),
        (None, "4", ["window_r=0.80", "zdim=2", "ind=4"], {"weight_decay": 0.05}),
        ("homo", None, ["window_r=0.80"], {"weight_decay": 0.05}),
    ],
    ids=["homo_with_wd", "het_through_cfgs", "homo_with_include_cfgs"],
)
class TestSetupThenRun:

    num_epochs = 5

    def check_outputs(
        self,
        setup_request,
        new_outdir=None,
        num_epochs=None,
        check_epochs=None,
        checkpoint=1,
        start_epoch=1,
        load=False,
    ):
        if new_outdir is None:
            out_files = set(os.listdir(setup_request.outdir))
        else:
            out_files = set(os.listdir(new_outdir))
        assert "config.yaml" in out_files

        if num_epochs is None:
            num_epochs = self.num_epochs
        if check_epochs is None:
            check_epochs = self.num_epochs * 2

        abinit = setup_request.poses is None
        abinit |= setup_request.pose_estimation == "abinit"
        four_ind = False
        if setup_request.ind is not None:
            four_ind |= "4" in os.path.basename(setup_request.ind)
        if setup_request.cfgs is not None:
            four_ind |= any(cfg.startswith("ind=4") for cfg in setup_request.cfgs)

        for epoch in range(0, check_epochs + 1):
            if epoch == num_epochs or (
                start_epoch <= epoch < num_epochs
                and (
                    (
                        setup_request.model in {None, "autodec"}
                        and (
                            epoch % checkpoint == 0
                            or abinit
                            and (
                                four_ind
                                and epoch == 3
                                or "5" in os.path.basename(setup_request.ind)
                                and epoch <= 2
                                or epoch <= 2
                            )
                        )
                    )
                    or (
                        setup_request.model == "autoenc"
                        and (epoch % checkpoint == 0 or (epoch - 1) % 2 == 0 and abinit)
                    )
                )
                or (
                    epoch == 0
                    and not load
                    and (setup_request.model in {None, "autodec"} or abinit)
                )
            ):
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

    def test_setup_directory(self, setup_request):
        """Create a reconstruction directory using the setup command."""
        setup_directory(setup_request)

    def test_then_train(self, setup_request):
        """Run the reconstruction experiment using the directory created by setup."""

        args = [
            setup_request.outdir,
            "--num-epochs",
            str(self.num_epochs),
            "--no-analysis",
            "--seed",
            "733",
        ]
        parser = argparse.ArgumentParser()
        train.add_args(parser)
        train.main(parser.parse_args(args))
        self.check_outputs(setup_request)

    @pytest.mark.parametrize("epoch, ksample", [(1, 3), (None, 4)])
    def test_then_analyze(self, setup_request, epoch, ksample):
        """Produce standard analyses for a particular epoch."""

        args = [setup_request.outdir]
        if epoch is not None:
            args += [str(epoch)]
        else:
            epoch = 6
        args += ["--ksample", str(ksample)]

        parser = argparse.ArgumentParser()
        cryodrgn.commands.analyze.add_args(parser)
        cryodrgn.commands.analyze.main(parser.parse_args(args))

        if setup_request.reconstruction_type == "het":
            anlzdir = os.path.join(setup_request.outdir, f"analyze.{epoch}")
            assert os.path.exists(anlzdir)
            assert os.path.exists(os.path.join(anlzdir, "z_pca.png"))
            assert os.path.exists(os.path.join(anlzdir, "pc2_10"))
            assert os.path.exists(os.path.join(anlzdir, f"kmeans{ksample}"))

    @pytest.mark.parametrize("checkpoint", [1, 2, 3])
    def test_then_rerun_with_checkpointing(self, setup_request, checkpoint):
        """Copy the experiment config file to an empty directory and run again."""

        rerun_dir = os.path.join(setup_request.outdir, "rerun")
        os.makedirs(rerun_dir)
        shutil.copyfile(
            os.path.join(setup_request.outdir, "config.yaml"),
            os.path.join(rerun_dir, "config.yaml"),
        )
        args = [
            rerun_dir,
            "--num-epochs",
            str(self.num_epochs),
            "--checkpoint",
            str(checkpoint),
            "--no-analysis",
            "--seed",
            "8990",
        ]
        parser = argparse.ArgumentParser()
        train.add_args(parser)
        train.main(parser.parse_args(args))
        self.check_outputs(setup_request, rerun_dir, checkpoint=checkpoint)
        shutil.rmtree(rerun_dir)

    @pytest.mark.parametrize("load_epoch", [None, 2, 4])
    def test_then_reload(self, setup_request, load_epoch):
        """Copy the experiment config file to an empty directory and run again."""
        args = [
            setup_request.outdir,
            "--num-epochs",
            str(self.num_epochs + 1),
            "--no-analysis",
            "--seed",
            "8990",
            "--load",
        ]
        if load_epoch is not None:
            args += [os.path.join(setup_request.outdir, f"weights.{load_epoch}.pkl")]

        parser = argparse.ArgumentParser()
        train.add_args(parser)
        train.main(parser.parse_args(args))

        out_files = os.listdir(setup_request.outdir)
        assert (
            f"weights.{self.num_epochs + 1}.pkl" in out_files
        ), f"Missing output model weights for epoch {self.num_epochs + 1}!"
        os.unlink(
            os.path.join(setup_request.outdir, f"weights.{self.num_epochs + 1}.pkl")
        )

        if setup_request.pose_estimation == "abinit" or setup_request.poses is None:
            assert (
                f"pose.{self.num_epochs + 1}.pkl" in out_files
            ), f"Missing output model poses for epoch {self.num_epochs + 1}!"
            os.unlink(
                os.path.join(setup_request.outdir, f"pose.{self.num_epochs + 1}.pkl")
            )

    def test_then_rerun_analyze(self, setup_request):
        """Produce analyses for restarted experiment after changing working dir."""
        args = [setup_request.outdir, "--ksample", "2"]
        parser = argparse.ArgumentParser()
        cryodrgn.commands.analyze.add_args(parser)

        curdir = os.getcwd()
        os.chdir(os.path.join(setup_request.outdir, ".."))
        cryodrgn.commands.analyze.main(parser.parse_args(args))
        os.chdir(curdir)

        if setup_request.reconstruction_type == "het":
            anlzdir = os.path.join(
                setup_request.outdir, f"analyze.{self.num_epochs + 1}"
            )
            assert os.path.exists(anlzdir)
            assert os.path.exists(os.path.join(anlzdir, "z_pca.png"))
            assert os.path.exists(os.path.join(anlzdir, "pc2_10"))
            assert os.path.exists(os.path.join(anlzdir, "kmeans2"))

    def test_new_outdir(self, setup_request, tmpdir_factory):
        """Copy configurations to a new output directory and run again."""
        new_outdir = tmpdir_factory.mktemp("new_outdir").strpath
        args = [
            new_outdir,
            "--from-outdir",
            setup_request.outdir,
            "--no-analysis",
            "--num-epochs",
            str(self.num_epochs - 1),
        ]
        parser = argparse.ArgumentParser()
        train.add_args(parser)
        train.main(parser.parse_args(args))

        self.check_outputs(setup_request, new_outdir, num_epochs=self.num_epochs - 1)
        shutil.rmtree(new_outdir)

    def test_load_and_new_outdir(self, setup_request, tmpdir_factory):
        """Copy configurations and run again, this time loading from the first run."""
        new_outdir = tmpdir_factory.mktemp("new_outdir_and_load").strpath
        args = [
            new_outdir,
            "--from-outdir",
            setup_request.outdir,
            "--no-analysis",
            "--load",
            "latest",
            "--num-epochs",
            str(self.num_epochs + 2),
        ]
        parser = argparse.ArgumentParser()
        train.add_args(parser)
        train.main(parser.parse_args(args))

        self.check_outputs(
            setup_request,
            new_outdir,
            num_epochs=self.num_epochs + 2,
            start_epoch=self.num_epochs + 1,
            load=True,
        )

        shutil.rmtree(new_outdir)
