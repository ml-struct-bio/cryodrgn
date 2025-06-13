"""Setup a cryoDRGN experiment

This command creates a new cryoDRGN experiment directory containing a `config.yaml`
file that stores the parameters used to run the experiment.
This folder can then be used to run the `cryodrgn train` command to
run the reconstruction model.

Example usage
-------------
# Setup cryoDRGN with fixed poses
$ cryodrgn setup outdir/ --particles particles.128.mrcs --ctf ctf.pkl --poses pose.pkl

# Setup cryoDRGN-AI with ab initio pose estimation
$ cryodrgn setup outdir/ --particles particles.128.mrcs --ctf ctf.pkl \
                         --model autodec --pose-estimation abinit

# Setup cryoDRGN-ET for heterogeneous subtomogram averaging
# Additional arguments to be added to the config.yaml file can be specified from
# the command line with --cfgs/-c
$ cryodrgn setup outdir/ --particles particles.128.star --datadir /path/to/subtilts \
                         --ctf ctf.pkl --poses pose.pkl --tilt \
                         --cfgs 't_extent=40' 'dose_per_tilt=3.5'

See also
--------
commands/train      `cryodrgn train` command used to train model after this command

"""
import os
import argparse
import logging
import cryodrgn.utils
from cryodrgn.models.utils import get_model_configurations, get_model_trainer_class

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    """The command-line arguments for use with `cryodrgn setup`."""

    parser.add_argument(
        "outdir",
        type=os.path.abspath,
        help="Path to the output directory to create for this cryoDRGN job",
    )

    inputs_group = parser.add_argument_group("Input Dataset")
    inputs_group.add_argument(
        "--particles",
        type=os.path.abspath,
        help="Path to the particle stack (.mrcs/.txt/.star/.cs)",
    )
    inputs_group.add_argument(
        "--ctf", type=os.path.abspath, help="Path to the CTF parameters (.pkl)"
    )
    inputs_group.add_argument(
        "--poses", type=os.path.abspath, help="Path to the poses (.pkl)"
    )
    inputs_group.add_argument("--ind", help="Path to filtering indices (.pkl)")
    inputs_group.add_argument(
        "--datadir",
        type=os.path.abspath,
        help="Path prefix to particle stack if loading relative "
        "paths from a .star or .cs file",
    )
    inputs_group.add_argument(
        "--dataset",
        help="Alternatively, specify a dataset name from an environmental "
        " variable `$CRYODRGN_DATASETS` that define dataset paths",
    )

    training_group = parser.add_argument_group("Training Settings")
    training_group.add_argument(
        "--model",
        "-m",
        choices=["autodec", "autoenc"],
        default="autoenc",
        help="Version of the cryoDRGN model to use: 'autoenc' for a standard cryoDRGN "
        "autoencoder, 'autodec' for a cryoDRGN-AI autodecoder. (default: %(default)s)",
    )
    training_group.add_argument(
        "--reconstruction-type",
        choices=["het", "homo"],
        default="het",
        help="Homogeneous or heterogeneous reconstruction (default: %(default)s)",
    )
    training_group.add_argument(
        "--zdim",
        type=int,
        default=8,
        help="Latent variable dimension if using heterogeneous reconstruction "
        "(default: %(default)s)",
    )
    training_group.add_argument(
        "--pose-estimation",
        choices=["abinit", "refine", "fixed"],
        default="fixed",
        help="Pose estimation procedure: `abinit` for ab initio reconstruction, "
        "`refine` to refine input poses by gradient descent or `fixed` to use fixed "
        "input poses (default: %(default)s)",
    )
    training_group.add_argument(
        "--tilt",
        action="store_true",
        help="Flag to specify if using a tilt series for subtomogram averaging",
    )
    training_group.add_argument(
        "--cfgs",
        "-c",
        nargs="+",
        default=list(),
        help="Additional configuration parameters to pass to the model "
        "in the form of 'CFG_KEY1=CFG_VAL1' 'CFG_KEY2=CFG_VAL2' ... ",
    )
    training_group.add_argument(
        "--include",
        "-i",
        type=os.path.abspath,
        help="Path to a .yaml file containing additional configuration parameters.",
    )


def main(args: argparse.Namespace) -> None:
    """Running the `cryodrgn setup` command (see `add_args` above for arguments)."""

    if os.path.exists(args.outdir) and os.listdir(args.outdir):
        raise ValueError(
            f"Output directory `{args.outdir}` already exists and is not empty!"
        )

    if args.reconstruction_type is None:
        args.reconstruction_type = "het" if args.zdim > 0 else "homo"
    elif args.reconstruction_type == "het":
        if not args.zdim > 0:
            raise ValueError(
                f"chosen zdim must be > 0 for heterogeneous reconstruction! "
                f"(got {args.zdim})"
            )
    elif args.reconstruction_type == "homo":
        if args.zdim != 0:
            logger.warning(
                f"zdim is set to {args.zdim}, but since using homogeneous "
                "reconstruction, this is ignored!"
            )
        args.zdim = 0
    else:
        raise ValueError(
            f"Unrecognized reconstruction type `{args.reconstruction_type}`!"
        )

    if args.pose_estimation in ("fixed", "refine"):
        if args.poses is None:
            raise ValueError(
                f"Must provide input poses for "
                f"--pose-estimation {args.pose_estimation} "
            )

    # handling different ways of specifying the input data, starting with a
    # file containing the data files
    if args.dataset is not None:  # TODO: We should document this somewhere
        paths_file = os.environ.get("CRYODRGN_DATASETS")
        data_paths = cryodrgn.utils.load_yaml(paths_file) if paths_file else None
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

    cfgs = {
        "model": args.model,
        "zdim": args.zdim,
        "pose_estimation": args.pose_estimation,
    }
    if args.tilt:
        cfgs["tilt"] = args.tilt  # True or False

    for data_label in ["particles", "poses", "ctf", "datadir", "ind"]:
        if getattr(args, data_label) is not None:  # From argparse
            cfgs[data_label] = getattr(args, data_label)
        elif these_paths.get(data_label) is not None:  # From CRYODRGN_DATASETS
            cfgs[data_label] = these_paths[data_label]

        # Check that input path is valid
        if data_label in cfgs and cfgs[data_label] is not None:
            if data_label != "ind" or not cfgs[data_label].isnumeric():
                assert os.path.exists(
                    cfgs[data_label]
                ), f"{cfgs[data_label]} does not exist"

    if args.include:
        include_cfgs = cryodrgn.utils.load_yaml(args.include)
        include_cfgs["model"] = args.model
        trainer_cls = get_model_trainer_class(include_cfgs)
        cfgs.update(trainer_cls.config_cls.parse_config(include_cfgs))

    configs = get_model_configurations(cfgs, add_cfgs=args.cfgs)
    os.makedirs(args.outdir, exist_ok=True)
    configs.write(os.path.join(args.outdir, "config.yaml"))
    print(f"Wrote {os.path.join(args.outdir, 'config.yaml')}")
