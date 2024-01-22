import argparse
from typing import Optional, Any
from cryodrgn.commands.train import main as training_main
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

    if "z_dim" in configs and configs["z_dim"] == 0:
        default_zdim = 8
        configs["z_dim"] = default_zdim

        print(
            "WARNING: given configurations specify homogeneous reconstruction, "
            "updating them to specify heterogeneous reconstruction "
            f"with z_dim={default_zdim}!"
        )

    training_main(args, configs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
