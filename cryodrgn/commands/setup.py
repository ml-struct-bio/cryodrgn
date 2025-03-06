"""Utilities for creating experiment output folders and configuration files."""

import os
import argparse
import cryodrgn.utils
from cryodrgn.models.utils import get_model_configurations


def add_args(parser: argparse.ArgumentParser) -> None:
    """The command-line arguments for use with `cryodrgn setup`."""

    parser.add_argument(
        "outdir",
        type=os.path.abspath,
        help="Path to the directory that will be created for this experiment",
    )

    parser.add_argument(
        "--model",
        "-m",
        default="cryodrgn-ai",
        choices=["cryodrgn-ai", "cryodrgn"],
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
        help="homogeneous or heterogeneous (default) reconstruction with z-dim=8?",
    )
    parser.add_argument(
        "--z-dim",
        type=int,
        help="homogeneous or heterogeneous (default) reconstruction when z-dim>0?",
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
        "--include",
        "-i",
        type=os.path.abspath,
        help="Path to a .yaml file containing additional configuration parameters.",
    )


def main(args: argparse.Namespace) -> None:
    """Running the `cryodrgn setup` command (see `add_args` above for arguments)."""

    if args.reconstruction_type is None:
        z_dim = int(args.z_dim) if args.z_dim is not None else 8
    else:
        if args.z_dim is not None:
            raise ValueError("Cannot specify both --reconstruction-type and --z-dim!")

        if args.reconstruction_type == "het":
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

    # handling different ways of specifying the input data, starting with a
    # file containing the data files
    paths_file = os.environ.get("CRYODRGN_DATASETS")
    data_paths = cryodrgn.utils.load_yaml(paths_file) if paths_file else None

    # TODO: merge with similar logic in reconstruction engine?
    if args.dataset is not None:
        if os.path.exists(args.dataset):
            these_paths = cryodrgn.utils.load_yaml(args.dataset)

            # resolve paths relative to the dataset file if they look relative
            for k in list(these_paths):
                if these_paths[k] and not os.path.isabs(these_paths[k]):
                    these_paths[k] = os.path.abspath(
                        os.path.join(args.dataset, these_paths[k])
                    )

        elif data_paths and args.dataset not in data_paths:
            raise ValueError(
                f"Given dataset {args.dataset} is not a "
                "label in the list of known datasets!"
            )

        elif data_paths is None:
            raise ValueError(
                "To specify datasets using a label, first specify"
                "a .yaml catalogue of datasets using the "
                "environment variable $CRYODRGN_DATASETS!"
            )

        # you can also give the dataset as a label in the global dataset list
        else:
            these_paths = data_paths[args.dataset]
    else:
        these_paths = {"particles": args.particles}

    cfgs = {"model": args.model, "z_dim": z_dim, "pose_estimation": pose_estimation}
    if args.tilt:
        cfgs["tilt"] = args.tilt

    for data_lbl in ["particles", "poses", "ctf", "datadir", "ind"]:
        if getattr(args, data_lbl) is not None:
            cfgs[data_lbl] = getattr(args, data_lbl)
        elif these_paths.get(data_lbl) is not None:
            cfgs[data_lbl] = these_paths[data_lbl]

    if args.include:
        cfgs.update(cryodrgn.utils.load_yaml(args.include))

    configs = get_model_configurations(cfgs, add_cfgs=args.cfgs)
    os.makedirs(args.outdir, exist_ok=True)
    configs.write(os.path.join(args.outdir, "config.yaml"))
