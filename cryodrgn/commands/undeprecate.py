"""Update existing experiment output directories to current format used by cryoDRGN.

Example usages
--------------

Scan a single directory:
$ cryodrgn undeprecate outdir

Scan every directory found in the current folder.
Note use of double-quotes to prevent bash from evaluating the glob expression itself.
$ cryodrgn undeprecate "*"

Scan every directory recursively found in the current folder
AND all of its children subdirectories.
$ cryodrgn undeprecate "**"

Recursively scan current folder and its children for cryoDRGN outputs
containing the substring "mike".
$ cryodrgn undeprecate "**/*mike*"

Recursively scan current folder and its children for cryoDRGN outputs containing the
substring "mike" anywhere under folders containing the substring "paper".
$ cryodrgn undeprecate "**/*paper*/**/*mike*"

"""
import os
import argparse
from pathlib import Path
import yaml
from cryodrgn.commands_utils.configs import check_open_config
from cryodrgn.commands.clean import clean_dir


def add_args(parser):
    parser.add_argument("outglob", help="Path pattern to scan for experiment outputs.")

    parser.add_argument(
        "--force",
        "-f",
        action="store_true",
        help="Update automatically without prompting the user.",
    )
    parser.add_argument(
        "--clean",
        "-c",
        nargs="?",
        const=5,
        help="Remove extraneous files from outputs while updating.",
    )
    parser.add_argument(
        "--dry-run",
        "-d",
        action="store_true",
        help="Only scan directories and identify their status, don't update.",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="count",
        default=0,
        help="Print more messages about ignored directories, etc.",
    )


def update_dir(d: Path, version: str, cfg: dict) -> None:
    if version == "3":
        new_cfg = {
            **cfg["dataset_args"],
            **cfg["model_args"],
            **cfg["lattice_args"],
            "model": "cryoDRGN_v3",
        }

        new_cfg["pose"] = new_cfg["poses"]
        del new_cfg["poses"]

        origwd = os.getcwd()
        os.chdir(d)

        for f in Path(".").iterdir():
            if f.name == "run.log":
                f.rename("training.log")

            elif f.name[:2] == "z." and f.name[-4:] == ".pkl":
                f.rename(f"conf.{f.name[2:]}")

        os.makedirs("out", exist_ok=True)
        for f in Path(".").iterdir():
            if not f.is_dir() or f.name != "out":
                f.rename("out" / f)

        with open(Path("out", "train-configs.yaml"), "w") as f:
            yaml.dump(new_cfg, f)

        os.chdir(origwd)

    elif version == "2":
        pass
    else:
        raise ValueError(f"Unrecognized version {version}")


def _prompt_dir(
    d: Path,
    cfg: dict,
    version_code: str,
    maxlen,
    args: argparse.Namespace,
) -> None:
    msg = f"is a cryoDRGNv{version_code} directory"

    if version_code == "4":
        msg = msg.replace("is a", "is already a")
    elif version_code == "N":
        msg = msg.replace("is a", "is not a")
        msg = msg.replace("cryoDRGNvN", "cryoDRGN")

    elif version_code == "2C":
        msg = msg.replace("is a", "might be a")
        msg = "".join([msg, ", but config.pkl is corrupted"])
    elif version_code == "3C":
        msg = msg.replace("is a", "might be a")
        msg = "".join([msg, ", but configs.yaml is corrupted"])

    if (
        args.dry_run
        or "C" in version_code
        or (args.verbose > 0 and version_code == "4")
        or (args.verbose > 1 and version_code == "N")
    ):
        print("\t".join(["".join(["`", str(d), "`"]).ljust(maxlen, "."), msg]))

    elif version_code in {"2", "3"}:
        prompt_msg = ', enter 1) "(s)kip" or 2) any other key to update:\n'
        prompt = input("".join(["`", str(d), " ", msg, prompt_msg]))

        if prompt not in {"s", "skip"}:
            update_dir(d, version_code, cfg)

            if args.clean:
                clean_dir(d, args.clean)


def main(args):
    scan_dirs = sorted(p for p in Path().glob(args.outglob) if p.is_dir())
    maxlen = len(str(max(scan_dirs, key=lambda d: len(str(d)))))

    while scan_dirs:
        cur_dir = scan_dirs.pop(0)
        dir_ls = set(Path.iterdir(cur_dir))

        if (
            Path(cur_dir, "configs.yaml") in dir_ls
            and Path(cur_dir, "out") in dir_ls
            and Path(cur_dir, "out").is_dir()
        ):
            cfg, version_code = check_open_config(cur_dir, "configs.yaml", "4")

        elif Path(cur_dir, "config.yaml") in dir_ls:
            cfg, version_code = check_open_config(cur_dir, "config.yaml", "3")
        elif Path(cur_dir, "config.yml") in dir_ls:
            cfg, version_code = check_open_config(cur_dir, "config.yml", "3")

        elif Path(cur_dir, "config.pkl") in dir_ls:
            cfg, version_code = check_open_config(cur_dir, "config.pkl", "2")

        else:
            cfg = None
            version_code = "N"

        _prompt_dir(cur_dir, cfg, version_code, maxlen, args)

        # don't scan subdirectories of already identified cryoDRGN folders
        if "N" not in version_code and "C" not in version_code:
            scan_dirs = [p for p in scan_dirs if cur_dir not in p.parents]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
