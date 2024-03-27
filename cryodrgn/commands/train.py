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
    parser.add_argument("outdir", help="experiment output location")

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


def main(args: argparse.Namespace, configs: Optional[dict[str, Any]] = None) -> None:
    if configs is None:
        configs = SetupHelper(args.outdir, update_existing=False).create_configs(
            model=args.model
        )

    cryodrgn.utils._verbose = False
    trainer = TRAINER_CLASSES[configs["model"]](configs)
    trainer.train()

    if args.do_analysis:
        ModelAnalyzer(trainer.outdir).analyze()
