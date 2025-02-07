"""Utilities for creating experiment output folders and configuration files."""

import argparse
from cryodrgn.trainers.amortinf_trainer import AmortizedInferenceTrainer
from cryodrgn.trainers.hps_trainer import HierarchicalPoseSearchTrainer

TRAINER_CLASSES = {
    "amort": AmortizedInferenceTrainer,
    "hps": HierarchicalPoseSearchTrainer,
}


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("config_file", help="experiment config file (.yaml)")

    parser.add_argument(
        "--model",
        "-m",
        default="amort",
        choices=["amort", "hps"],
        help="which generation of cryoDRGN learning models to apply",
    )

    parser.add_argument("--dataset", help="which dataset to run the experiment on")
    parser.add_argument(
        "--particles", help="path to the picked particles (.mrcs/.star /.txt)"
    )
    parser.add_argument("--ctf", help="path to the CTF parameters (.pkl)")
    parser.add_argument("--poses", help="path to the poses (.pkl)")
    parser.add_argument("--ind", help="path to filtering indices (.pkl)")

    parser.add_argument(
        "--reconstruction-type",
        choices=["het", "homo", None],
        default=None,
        help="homogeneous or heterogeneous reconstruction with z-dim=8?",
    )
    parser.add_argument(
        "--z_dim",
        type=int,
        default=0,
        help="homogeneous (default) or heterogeneous reconstruction with z-dim=8?",
    )
    parser.add_argument(
        "--pose-estimation",
        default="abinit",
        choices=["abinit", "refine", "fixed"],
        help="`abinit` for no initialization (default), `refine` to refine "
        "ground truth poses by gradient descent or `fixed` to use ground "
        "truth poses",
    )
    parser.add_argument(
        "--tilt",
        action="store_true",
        help="specify if using a tilt series",
    )

    parser.add_argument(
        "--cfgs",
        "-c",
        nargs="+",
        help="additional configuration parameters to pass to the model "
        "in the form of 'CFG_KEY1=CFG_VAL1' 'CFG_KEY2=CFG_VAL2' ... ",
    )


def main(args: argparse.Namespace) -> None:
    trainer_cls = TRAINER_CLASSES[args.model]

    if args.reconstruction_type == "het":
        z_dim = 8
    elif args.reconstruction_type == "homo":
        z_dim = 0
    elif args.reconstruction_type is None:
        z_dim = int(args.z_dim)
    else:
        raise ValueError(
            f"Unrecognized reconstruction type `{args.reconstruction_type}`!"
        )

    cfgs = {
        "model": args.model,
        "dataset": args.dataset,
        "particles": args.particles,
        "ctf": args.ctf,
        "poses": args.poses,
        "ind": args.ind,
        "z_dim": z_dim,
        "pose_estimation": args.pose_estimation,
    }
    cfgs.update(trainer_cls.config_cls.parse_cfg_keys(args.cfgs))
    configs = trainer_cls.config_cls(**cfgs)
    configs.write(args.config_file)
