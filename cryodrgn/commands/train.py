"""Train a cryoDRGN model to perform reconstruction on a particle stack.

Example usages
--------------

"""
import argparse
from cryodrgn.models.config import load_configs
from cryodrgn.commands.train_ai import main as run_train_ai
from cryodrgn.commands.train_vae import main as run_train_vae
from cryodrgn.commands.train_nn import main as run_train_nn


def add_args(parser):
    parser.add_argument("outdir", help="experiment output location")

    parser.add_argument(
        "--no-analysis",
        action="store_true",
        help="just do the training stage",
        dest="no_anlz",
    )


def main(args):
    cfg = load_configs(args.outdir)

    if "model" in cfg:
        if cfg["model"] == "AI":
            run_train_ai(args, cfg)
        elif cfg["model"] == "VAE":
            run_train_vae(args, cfg)
        elif cfg["model"] == "NN":
            run_train_nn(args, cfg)

        else:
            raise ValueError(f'Unrecognized model `{cfg["model"]}`!')

    else:
        run_train_ai(args, cfg)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
