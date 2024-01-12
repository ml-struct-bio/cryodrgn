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
import argparse
from pathlib import Path
import pickle
import yaml

NORMAL_EXCEPTIONS = (
    EOFError,
    BufferError,
    pickle.UnpicklingError,
    ImportError,
    IndexError,
    AttributeError,
)


def add_args(parser):
    parser.add_argument("outglob", help="Path pattern to scan for experiment outputs.")

    parser.add_argument(
        "--force",
        "-f",
        action="store_true",
        help="Update automatically without prompting the user.",
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


def update_dir(d: Path, version: str) -> None:
    if version == "3":
        pass
    elif version == "2":
        pass
    else:
        raise ValueError(f"Unrecognized version {version}")


def main(args):
    scan_dirs = {p for p in Path().glob(args.outglob) if p.is_dir()}
    maxlen = len(str(max(scan_dirs, key=lambda d: len(str(d)))))

    def prompt_dir(d: Path, version_code: str) -> None:
        msg = "".join(f"is a cryoDRGNv{version_code} directory")

        if version_code == "4":
            msg = msg.replace("is a", "is already a")
        elif version_code == "N":
            msg = msg.replace("is a", "is not a")

        elif version_code == "2C":
            msg = msg.replace("is a", "might be a")
            msg = "".join([msg, ", but config.pkl is corrupted"])
        elif version_code == "3C":
            msg = msg.replace("is a", "might be a")
            msg = "".join([msg, ", but configs.yaml is corrupted"])

        if args.dry_run:
            print("\t".join(["".join(["`", str(d), "`"]).ljust(maxlen, "."), msg]))

        else:
            if version_code == "4":
                print(" ".join(["".join(["`", str(d), "`"]), msg]))

            else:
                prompt = input(
                    "".join(
                        [
                            "`",
                            str(d),
                            " ",
                            msg,
                            ', enter 1) "(s)kip" or 2) any other key to update:\n',
                        ]
                    )
                )

                if prompt not in {"s", "skip"}:
                    update_dir(d, version_code)

    while scan_dirs:
        cur_dir = scan_dirs.pop()
        is_cryodrgn = True
        dir_ls = set(Path.iterdir(cur_dir))

        if Path(cur_dir, "configs.yaml") in dir_ls:
            try:
                with open(Path(cur_dir, "configs.yaml"), "r") as f:
                    cfg = yaml.safe_load(f)
            except NORMAL_EXCEPTIONS:
                is_cryodrgn = False
                cfg = dict()

                if args.verbose > 1:
                    prompt_dir(cur_dir, "3C")

            if "particles" in cfg or "dataset" in cfg:
                if Path(cur_dir, "out") in dir_ls and Path(cur_dir, "out").is_dir():
                    if args.verbose > 0:
                        prompt_dir(cur_dir, "4")

                else:
                    prompt_dir(cur_dir, "3")

        elif Path(cur_dir, "config.pkl") in dir_ls:
            try:
                with open(Path(cur_dir, "config.pkl"), "rb") as f:
                    cfg = pickle.load(f)
            except NORMAL_EXCEPTIONS:
                is_cryodrgn = False
                cfg = dict()

                if args.verbose > 1:
                    prompt_dir(cur_dir, "2C")

            if "particles" in cfg or "dataset" in cfg:
                prompt_dir(cur_dir, "2")

        else:
            is_cryodrgn = False
            if args.verbose > 1:
                prompt_dir(cur_dir, "N")

        # don't scan subdirectories of already identified cryoDRGN folders
        if is_cryodrgn:
            scan_dirs -= {p for p in scan_dirs if cur_dir in p.parents}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
