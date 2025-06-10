"""Train a volume reconstruction model described by a config file in an output folder.

The `train` command is used to run a reconstruction experiment using the configurations
previously saved using a command such as `cryodrgn setup`.

This command can be used to train any model included in cryoDRGN, including:
    - homogeneous and heterogeneous reconstruction with given poses (cryoDRGN v1)
    - homo and het reconstruction with ab initio poses (cryoDRGN v2)
    - reconstruction using tilt series particle stacks (cryoDRGN v3)
    - ab-initio pose estimation using cryoDRGN-AI (cryoDRGN v4)

Example usage
--------------
$ cryodrgn train new-test

# Submit reconstruction task to a compute cluster with Slurm installed
$ sbatch -t 3:00:00 --mem=16G -o new-test.out --wrap='cryodrgn train new-test'

See also
--------
commands/setup      `cryodrgn setup` command used as precursor to this command

"""
import os
import shutil
import argparse
from typing import Optional, Any
import cryodrgn.config
import cryodrgn.utils
from cryodrgn.models.utils import get_model_trainer_class
from cryodrgn.commands.analyze import ModelAnalyzer


def add_args(parser: argparse.ArgumentParser) -> None:
    """The command-line arguments for use with the `cryodrgn train` command."""

    outputs_group = parser.add_argument_group("Managing Model Outputs and Logging")
    outputs_group.add_argument(
        "outdir",
        type=os.path.abspath,
        help="Path to the directory that has been created for this experiment",
    )
    outputs_group.add_argument(
        "--checkpoint",
        type=int,
        help="Checkpointing interval in N_EPOCHS (default: %(default)s)",
    )
    outputs_group.add_argument(
        "--log-interval",
        type=int,
        help="Logging interval in N_IMGS",
    )
    outputs_group.add_argument(
        "-v", "--verbose", action="store_true", help="Increase verbosity"
    )

    loading_group = parser.add_argument_group("Managing Dataset Loading and Parsing")
    loading_group.add_argument(
        "--lazy",
        action="store_true",
        help="Lazy loading if full dataset is too large to fit in memory",
    )
    loading_group.add_argument(
        "--shuffler-size",
        type=int,
        help="If non-zero, will use a data shuffler for faster lazy data loading.",
    )
    loading_group.add_argument(
        "--max-threads",
        type=int,
        help="Maximum number of CPU cores for data loading",
    )

    training_group = parser.add_argument_group("Reconstruction Training Settings")
    training_group.add_argument(
        "-n",
        "--num-epochs",
        type=int,
        help="Number of training epochs (default: %(default)s)",
    )
    training_group.add_argument(
        "--load",
        "-l",
        nargs="?",
        const="latest",
        help="Initialize training from a checkpoint; "
        "give with no arguments or with `latest` to load latest available epoch",
    )
    training_group.add_argument(
        "-b",
        "--batch-size",
        type=int,
        help="Minibatch size (default: %(default)s)",
    )
    training_group.add_argument(
        "--no-amp",
        action="store_false",
        dest="amp",
        help="Do not use mixed-precision training for accelerating training",
    )
    training_group.add_argument(
        "--multigpu",
        action="store_true",
        help="Parallelize training across all detected GPUs",
    )
    training_group.add_argument(
        "--seed", type=int, help="Random seed for use with e.g. numpy"
    )
    training_group.add_argument(
        "--cfgs",
        "-c",
        nargs="+",
        help="additional configuration parameters to pass to the model "
        "in the form of 'CFG_KEY1=CFG_VAL1' 'CFG_KEY2=CFG_VAL2' ... ",
    )
    training_group.add_argument(
        "--include",
        "-i",
        type=os.path.abspath,
        help="Path to a .yaml file containing additional configuration parameters.",
    )

    parser.add_argument(
        "--from-outdir",
        "-f",
        type=os.path.abspath,
        help="Copy configs from this directory instead of using any configs "
        "found in the primary output directory.",
    )
    parser.add_argument(
        "--no-analysis",
        action="store_false",
        dest="do_analysis",
        help="just do the training stage",
    )


def main(
    args: argparse.Namespace, additional_configs: Optional[dict[str, Any]] = None
) -> None:
    """Running the `cryodrgn train` command (see `add_args` above for arguments).

    An additional `configs` dictionary of configuration values can also be passed, which
    will be appended to (and override) the values in the given config file and any
    values specified through the `--cfgs` command-line argument.

    """
    if args.from_outdir:
        from_cfgs = os.path.join(args.from_outdir, "config.yaml")
        if not os.path.exists(from_cfgs):
            raise ValueError(
                f"No CryoDRGN configuration file `config.yaml` "
                f"found in `{args.from_outdir=}`!"
            )
        os.makedirs(args.outdir, exist_ok=True)
        shutil.copy(from_cfgs, os.path.join(args.outdir, "config.yaml"))

    configs = cryodrgn.utils.load_yaml(os.path.join(args.outdir, "config.yaml"))
    trainer_cls = get_model_trainer_class(configs)
    configs = trainer_cls.config_cls.parse_config(configs)

    if additional_configs is not None:
        configs.update(trainer_cls.config_cls.parse_config(additional_configs))
    if args.include:
        configs.update(
            trainer_cls.config_cls.parse_config(cryodrgn.utils.load_yaml(args.include))
        )
    if args.cfgs:
        configs.update(trainer_cls.config_cls.parse_cfg_keys(args.cfgs))

    if args.seed is not None:
        configs["seed"] = args.seed
    train_args = configs["train_args"] if "train_args" in configs else configs
    model_args = configs["model_args"] if "model_args" in configs else configs

    if args.load is not None:
        if args.from_outdir is not None:
            train_args["load"] = cryodrgn.utils.find_latest_output(args.from_outdir)
            if model_args["pose_estimation"] == "abinit":
                train_args["load_poses"] = cryodrgn.utils.find_latest_output(
                    args.from_outdir, outlbl="pose"
                )
        else:
            train_args["load"] = args.load

    if args.num_epochs is not None:
        train_args["num_epochs"] = args.num_epochs
    if args.checkpoint is not None:
        train_args["checkpoint"] = args.checkpoint
    if args.log_interval is not None:
        train_args["log_interval"] = args.log_interval
    if args.lazy:
        train_args["lazy"] = True
    if args.multigpu:
        train_args["multigpu"] = True
    if args.shuffler_size is not None:
        train_args["shuffler_size"] = args.shuffler_size
    if args.max_threads is not None:
        train_args["max_threads"] = args.max_threads
    if args.batch_size is not None:
        train_args["batch_size"] = args.batch_size
    if not args.amp:
        train_args["amp"] = False

    cryodrgn.utils._verbose = False
    trainer = trainer_cls(configs, outdir=args.outdir)
    trainer.train()

    if args.do_analysis:
        ModelAnalyzer(trainer.outdir).analyze()
