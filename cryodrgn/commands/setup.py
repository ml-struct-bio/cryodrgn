"""Utilities for creating experiment output folders and configuration files."""

import argparse
import pandas as pd
from cryodrgn.utils import save_yaml
from cryodrgn.models.utils import get_model_configurations


def add_args(parser: argparse.ArgumentParser) -> None:
    """The command-line arguments for use with `cryodrgn setup`."""

    parser.add_argument(
        "config_file",
        nargs="?",
        help="experiment config file (.yaml); if not given, "
        "will print configurations to screen",
    )

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
        "--datadir",
        help="Path prefix to particle stack if loading relative "
        "paths from a .star or .cs file",
    )

    parser.add_argument(
        "--reconstruction-type",
        choices=["het", "homo"],
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
        default=list(),
        help="additional configuration parameters to pass to the model "
        "in the form of 'CFG_KEY1=CFG_VAL1' 'CFG_KEY2=CFG_VAL2' ... ",
    )
    parser.add_argument(
        "--defaults",
        action="store_true",
        help="get the default values of all configuration parameters "
        "used in this engine",
    )


def main(args: argparse.Namespace) -> None:
    """Running the `cryodrgn setup` command (see `add_args` above for arguments)."""

    if args.reconstruction_type is None:
        z_dim = int(args.z_dim)
    elif args.reconstruction_type == "het":
        z_dim = 8
    elif args.reconstruction_type == "homo":
        z_dim = 0
    else:
        raise ValueError(
            f"Unrecognized reconstruction type `{args.reconstruction_type}`!"
        )
    if args.pose_estimation is None:
        pose_estimation = "abinit" if args.poses is None else "fixed"
    else:
        pose_estimation = args.pose_estimation

    cfgs = {
        "model": args.model,
        "particles": args.particles,
        "z_dim": z_dim,
        "pose_estimation": pose_estimation,
    }
    if args.datadir:
        cfgs["datadir"] = args.datadir
    if args.ctf:
        cfgs["ctf"] = args.ctf
    if args.poses:
        cfgs["poses"] = args.poses
    if args.ind:
        cfgs["ind"] = args.ind
    if args.dataset:
        cfgs["dataset"] = args.dataset
    if args.tilt:
        cfgs["tilt"] = True

    configs = get_model_configurations(cfgs, add_cfgs=args.cfgs)
    if args.defaults:
        if args.config_file:
            configs.write(args.config_file)
        else:
            print(configs)
    else:
        cfgs.update(configs.parse_cfg_keys(args.cfgs))
        if args.config_file:
            save_yaml(cfgs, args.config_file)
        else:
            print(pd.Series(cfgs).to_string())
