"""Create, train, and analyze a reconstruction experiment.

This command can be used to train any model included in cryoDRGN, including:
    - homogeneous and heterogeneous reconstruction with given poses (cryoDRGN v1)
    - homo and het reconstruction with ab initio poses (cryoDRGN v2)
    - reconstruction using tilt series particle stacks (cryoDRGN v3)
    - pose estimation using amortized inference (drgnai; cryoDRGN v4)

Example usages
--------------

$ cryodrgn train new-test --model=cryodrgn
$ sbatch -t 3:00:00 --wrap='cryodrgn train new-test' --mem=16G -o new-test.out

"""
import os
import argparse
from typing import Optional, Any
import cryodrgn.config
import cryodrgn.utils
from cryodrgn.models.utils import get_model_trainer
from cryodrgn.commands.analyze import ModelAnalyzer


def add_args(parser: argparse.ArgumentParser) -> None:
    """The command-line arguments for use with the `cryodrgn train` command."""

    parser.add_argument(
        "outdir",
        type=os.path.abspath,
        help="Path to the directory that has been created for this experiment",
    )

    parser.add_argument(
        "-n",
        "--num-epochs",
        type=int,
        help="Number of training epochs (default: %(default)s)",
    )
    parser.add_argument(
        "--checkpoint",
        type=int,
        help="Checkpointing interval in N_EPOCHS (default: %(default)s)",
    )
    parser.add_argument(
        "--log-interval",
        type=int,
        help="Logging interval in N_IMGS",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Increase verbosity"
    )
    parser.add_argument("--seed", type=int, help="Random seed for use with e.g. numpy")

    parser.add_argument(
        "--lazy",
        action="store_true",
        help="Lazy loading if full dataset is too large to fit in memory",
    )
    parser.add_argument(
        "--shuffler-size",
        type=int,
        help="If non-zero, will use a data shuffler for faster lazy data loading.",
    )
    parser.add_argument(
        "--max-threads",
        type=int,
        help="Maximum number of CPU cores for data loading",
    )

    parser.add_argument(
        "-b",
        "--batch-size",
        type=int,
        help="Minibatch size (default: %(default)s)",
    )

    parser.add_argument(
        "--no-amp",
        action="store_false",
        dest="amp",
        help="Do not use mixed-precision training for accelerating training",
    )
    parser.add_argument(
        "--multigpu",
        action="store_true",
        help="Parallelize training across all detected GPUs",
    )

    parser.add_argument(
        "--cfgs",
        "-c",
        nargs="+",
        help="additional configuration parameters to pass to the model "
        "in the form of 'CFG_KEY1=CFG_VAL1' 'CFG_KEY2=CFG_VAL2' ... ",
    )
    parser.add_argument(
        "--include",
        "-i",
        type=os.path.abspath,
        help="Path to a .yaml file containing additional configuration parameters.",
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
    configs = cryodrgn.utils.load_yaml(os.path.join(args.outdir, "config.yaml"))

    if additional_configs is not None:
        configs.update(additional_configs)
    if args.include:
        configs.update(cryodrgn.utils.load_yaml(args.include))
    if args.seed is not None:
        configs["seed"] = args.seed

    train_args = configs["train_args"] if "train_args" in configs else configs
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

    trainer = get_model_trainer(configs, outdir=args.outdir, add_cfgs=args.cfgs)
    cryodrgn.utils._verbose = args.verbose
    trainer.train()

    if args.do_analysis:
        ModelAnalyzer(trainer.outdir).analyze()
