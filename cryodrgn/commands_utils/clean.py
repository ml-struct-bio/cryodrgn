"""Remove extraneous files from experiment output directories

This utility removes output files from cryoDRGN output directories
except for those belonging to every `n` epoch, and those whose epochs have analysis
folders associated with them. This includes files containing weights, poses,
reconstructions, and latent variables.

Running this script without any positional arguments will just clean the current
directory. You can also provide any number of directories instead, in which case the
script will scan each of them. Finally, any directory can also be given
as a glob wildcard expression, and the script will scan any directory matching it. In
this case the script will also force the user to confirm cleaning of any directory
that looks like cryoDRGN output unless additional arguments are provided (see below).

Example usage
-------------

@ Scan current directory
$ cryodrgn_utils clean

# See how many files would be removed with every 3 epochs but don't remove them
$ cryodrgn_utils clean -n 3 -d

# Scan a single directory
$ cryodrgn_utils clean outdir

# Scan multiple directories
$ cryodrgn_utils clean outdir outdir2

# Scan every directory found in the current folder
# Note use of double-quotes to prevent bash from evaluating the glob expression itself
$ cryodrgn_utils clean "*"

# Recursively scan current folder and its children for cryoDRGN outputs
# containing the substring "old"
$ cryodrgn_utils clean "**/*old*"

"""
import os
import argparse
from pathlib import Path
import yaml
import pickle
from itertools import product
from math import log10

NORMAL_EXCEPTIONS = (
    EOFError,
    BufferError,
    pickle.UnpicklingError,
    ImportError,
    IndexError,
    AttributeError,
)


def add_args(parser):
    parser.add_argument(
        "outglobs",
        nargs="*",
        help="Path patterns to scan for experiment outputs."
        "Must be relative to current directory.",
    )

    parser.add_argument(
        "--every-n-epochs",
        "-n",
        type=int,
        default=5,
        help="Only save output from every certain number of epochs.",
    )
    parser.add_argument(
        "--force",
        "-f",
        action="count",
        default=0,
        help="Clean automatically without prompting the user.",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="count",
        default=0,
        help="Print more messages about ignored directories, etc.",
    )
    parser.add_argument(
        "--dry-run",
        "-d",
        action="store_true",
        help="Only scan directories and identify their status, don't update.",
    )


def clean_dir(d: Path, args: argparse.Namespace) -> None:
    analysis_epochs = {
        int(out_dir.name.split(".")[1])
        for out_dir in d.iterdir()
        if (
            out_dir.is_dir()
            and out_dir.name[:8] == "analyze."
            and out_dir.name.split(".")[1].isnumeric()
        )
    }

    rmv_count = 0
    for out_lbl in ["weights", "z", "pose", "reconstruct"]:
        for out_fl in d.glob(f"{out_lbl}.*"):
            epoch = out_fl.name.split(".")[1]

            if out_fl.is_file() and (
                epoch.isnumeric()
                and (int(epoch) % args.every_n_epochs) != 0
                and int(epoch) not in analysis_epochs
            ):
                rmv_count += 1
                if not args.dry_run:
                    os.remove(out_fl)

    if args.dry_run:
        print(f"\tWould remove {rmv_count} files!")
    else:
        print(f"\tRemoved {rmv_count} files!")


def _prompt_dir(
    d: Path,
    cfg: dict,
    maxlen: int,
    args: argparse.Namespace,
) -> None:
    msg = "is a cryoDRGN directory"

    if not cfg:
        msg = msg.replace("is a", "is not a")

    if (cfg and (args.dry_run or args.force)) or (not cfg and args.verbose > 0):
        print("".join(["".join(["`", str(d), "`"]).ljust(maxlen + 7, "."), msg]))

    if cfg:
        if args.force or args.dry_run:
            clean_dir(d, args)

        else:
            prompt_msg = ', enter 1) "(s)kip" or 2) any other key to clean:'
            prompt = input("".join(["`", str(d), " ", msg, prompt_msg]))

            if prompt not in {"s", "skip"}:
                clean_dir(d, args)


def check_open_config(d: Path) -> dict:
    """Safely gets configs from a potential cryoDRGN directory."""
    dir_files = {p.name for p in d.iterdir() if p.is_file()}
    config = None
    cfg = None

    for cfg_lbl, cfg_ext in product(["config", "configs"], [".yaml", ".yml"]):
        if "".join([cfg_lbl, cfg_ext]) in dir_files:
            config = "".join([cfg_lbl, cfg_ext])
            break

    if config:
        corrupted = False

        try:
            with open(Path(d, config), "r") as f:
                cfg = yaml.safe_load(f)

        except NORMAL_EXCEPTIONS:
            cfg = None
            corrupted = True

        if not corrupted:
            cryodrgn_cmds = {"abinit_homo", "abinit_het", "train_nn", "train_vae"}
            if (
                "cmd" not in cfg
                or "cryodrgn" not in cfg["cmd"][0]
                or cfg["cmd"][1] not in cryodrgn_cmds
                or "dataset_args" not in cfg
            ):
                cfg = dict()

    return cfg


def main(args):
    if not args.outglobs:
        clean_dir(Path(os.getcwd()), args)

    else:
        for outglob in args.outglobs:
            scan_dirs = sorted(set(p for p in Path().glob(outglob) if p.is_dir()))

            if len(scan_dirs) == 1:
                clean_dir(Path(tuple(scan_dirs)[0]), args)

            else:
                maxlen = len(str(max(scan_dirs, key=lambda d: len(str(d)))))
                print(f"Found {len(scan_dirs)} directories matching {outglob}")

                if len(scan_dirs) > 10**args.force:
                    if args.force:
                        force_str = "f" * (int(log10(len(scan_dirs))) + 1)
                        print(
                            f"Reverting to interactivity — use -{force_str} to "
                            f"disable interactivity in this case!"
                        )
                        args.force = False

                while scan_dirs:
                    cur_dir = scan_dirs.pop(0)
                    cfg = check_open_config(cur_dir)
                    _prompt_dir(cur_dir, cfg, maxlen, args)

                    # don't scan subdirectories of already identified cryoDRGN folders
                    if cfg:
                        scan_dirs = [p for p in scan_dirs if cur_dir not in p.parents]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
