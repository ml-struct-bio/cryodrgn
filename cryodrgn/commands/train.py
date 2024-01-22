import argparse
from typing import Optional, Any
import cryodrgn.utils
from cryodrgn.models.amortinf_trainer import AmortizedInferenceTrainer
from cryodrgn.models.hps_trainer import HierarchicalSearchTrainer
from cryodrgn.commands.analyze import ModelAnalyzer
from cryodrgn.commands.setup import SetupHelper


def add_args(parser):
    parser.add_argument("outdir", help="experiment output location")

    parser.add_argument(
        "--no-analysis",
        action="store_true",
        help="just do the training stage",
    )


def main(args: argparse.Namespace, configs: Optional[dict[str, Any]] = None):
    if configs is None:
        configs = SetupHelper(args.outdir, update_existing=False).create_configs()

    cryodrgn.utils._verbose = False

    if configs["model"] == "ai":
        trainer = AmortizedInferenceTrainer(configs)
    elif configs["model"] == "hps":
        trainer = HierarchicalSearchTrainer(configs)

    else:
        raise ValueError(f"Unrecognized model: {configs['model']}")

    trainer.train()
    if not args.no_analysis:
        ModelAnalyzer(configs).analyze()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
