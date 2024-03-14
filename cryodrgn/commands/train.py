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
from cryodrgn.trainers.amortinf_trainer import AmortizedInferenceTrainer
from cryodrgn.trainers.hps_trainer import HierarchicalPoseSearchTrainer
from cryodrgn.commands.analyze import ModelAnalyzer
from cryodrgn.commands.setup import SetupHelper


def add_args(parser):
    parser.add_argument("outdir", help="experiment output location")

    parser.add_argument(
        "--model",
        choices={"hps", "amort"},
        help="which model to use for reconstruction"
        )
    parser.add_argument(
        "--no-analysis",
        action="store_false",
        dest="do_analysis",
        help="just do the training stage",
    )

    return parser


def main(args: argparse.Namespace, configs: Optional[dict[str, Any]] = None):
    if configs is None:
        configs = SetupHelper(args.outdir, update_existing=False).create_configs(
            model=args.model
        )

    cryodrgn.utils._verbose = False

    if configs["model"] == "amort":
        trainer = AmortizedInferenceTrainer(configs)
    elif configs["model"] == "hps":
        trainer = HierarchicalPoseSearchTrainer(configs)

    else:
        raise ValueError(f"Unrecognized model: {configs['model']}")

    trainer.train()
    if args.do_analysis:
        ModelAnalyzer(configs).analyze()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    main(add_args(parser).parse_args())
