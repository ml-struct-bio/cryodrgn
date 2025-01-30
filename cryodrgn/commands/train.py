"""Create, train, and analyze a reconstruction experiment.

This command can be used to train any model included in cryoDRGN, including:
    - homogeneous and heterogeneous reconstruction with given poses (cryoDRGN v1)
    - homo and het reconstruction with ab initio poses (cryoDRGN v2)
    - reconstruction using tilt series particle stacks (cryoDRGN v3)
    - pose estimation using amortized inference (drgnai; cryoDRGN v4)

Example usages
--------------

$ cryodrgn train new-test --model=hps
$ sbatch -t 3:00:00 --wrap='cryodrgn train new-test' --mem=16G -o new-test.out

"""
import argparse
from typing import Optional, Any
import cryodrgn.utils
from cryodrgn.commands.analyze import ModelAnalyzer
from cryodrgn.commands.setup import SetupHelper
from cryodrgn.trainers.amortinf_trainer import AmortizedInferenceTrainer
from cryodrgn.trainers.hps_trainer import HierarchicalPoseSearchTrainer

TRAINER_CLASSES = {
    "amort": AmortizedInferenceTrainer,
    "hps": HierarchicalPoseSearchTrainer,
}


def add_args(parser: argparse.ArgumentParser) -> None:
    """The command-line arguments for use with the `cryodrgn train` command."""

    parser.add_argument("config_file", help="experiment config file")
    parser.add_argument(
        "--outdir", "-o", required=True, help="where to store experiment output"
    )

    parser.add_argument(
        "--model",
        choices=set(TRAINER_CLASSES),
        help="which model to use for reconstruction",
    )
    parser.add_argument(
        "--no-analysis",
        action="store_false",
        dest="do_analysis",
        help="just do the training stage",
    )
    parser.add_argument(
        "--cfgs",
        "-c",
        nargs="+",
        help="additional configuration parameters to pass to the model "
        "in the form of 'CFG_KEY1=CFG_VAL1' 'CFG_KEY2=CFG_VAL2' ... ",
    )


def main(args: argparse.Namespace, configs: Optional[dict[str, Any]] = None) -> None:
    """Running the `cryodrgn train` command (see `add_args` above for arguments).

    An additional `configs` dictionary of configuration values can also be passed, which
    will be appended to (and override) the values in the given config file and any
    values specified through the `--cfgs` command-line argument.

    """
    if configs is None:
        configs = dict()

    file_configs = SetupHelper(args.config_file, update_existing=False).create_configs(
        model=args.model,
    )

    trainer_cls = TRAINER_CLASSES[configs["model"]]
    if args.cfgs:
        configs = {**configs, **trainer_cls.config_cls.parse_cfg_keys(args.cfgs)}

    configs = {**file_configs, **configs, "outdir": args.outdir}
    trainer = trainer_cls(configs)
    cryodrgn.utils._verbose = False
    trainer.train()

    if args.do_analysis:
        ModelAnalyzer(trainer.outdir).analyze()
