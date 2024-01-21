"""Remove extraneous files from experiment output directories

Example usages
--------------

Scan a single directory:
$ cryodrgn clean outdir

Scan every directory found in the current folder.
Note use of double-quotes to prevent bash from evaluating the glob expression itself.
$ cryodrgn clean "*"

Scan every directory recursively found in the current folder
AND all of its children subdirectories.
$ cryodrgn clean "**"

Recursively scan current folder and its children for cryoDRGN outputs
containing the substring "mike".
$ cryodrgn clean "**/*mike*"

Recursively scan current folder and its children for cryoDRGN outputs containing the
substring "mike" anywhere under folders containing the substring "paper".
$ cryodrgn clean "**/*paper*/**/*mike*"

"""
import os
import argparse
from pathlib import Path
from cryodrgn.commands_utils.configs import check_open_config


def add_args(parser):
    parser.add_argument("outglob", help="Path pattern to scan for experiment outputs.")

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
        action="store_true",
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


def clean_dir(d: Path, every_n_epochs: int = 5) -> None:
    analysis_epochs = {
        int(out_dir.name.split(".")[1])
        for out_dir in Path(d / "out").iterdir()
        if (
            out_dir.is_dir()
            and out_dir.name[:8] == "analyze."
            and out_dir.name.split(".")[1].isnumeric()
        )
    }

    for out_lbl in ["pose", "conf", "reconstruct", "weights"]:
        out_fls = tuple(Path(d / "out").glob(f"{out_lbl}.*"))

        if out_fls:
            all_epochs = [
                int(out_fl.name.split(".")[1])
                for out_fl in out_fls
                if out_fl.name.split(".")[1].isnumeric()
            ]
            max_epoch = max([0] + all_epochs)

            for out_fl in out_fls:
                epoch = out_fl.name.split(".")[1]

                if out_fl.is_file() and (
                    not epoch.isnumeric()
                    or (
                        0 < int(epoch) < max_epoch
                        and (int(epoch) % every_n_epochs) != 0
                        and int(epoch) not in analysis_epochs
                    )
                ):
                    os.remove(out_fl)


def _prompt_dir(
    d: Path,
    cfg: dict,
    version_code: str,
    maxlen: int,
    args: argparse.Namespace,
) -> None:
    msg = "is a cryoDRGN directory"

    if version_code == "N":
        msg = msg.replace("is a", "is not a")

    if args.dry_run or (args.verbose > 1 and version_code == "N"):
        print("\t".join(["".join(["`", str(d), "`"]).ljust(maxlen, "."), msg]))

    elif version_code == "4":
        prompt_msg = ', enter 1) "(s)kip" or 2) any other key to clean:\n'
        prompt = input("".join(["`", str(d), " ", msg, prompt_msg]))

        if prompt not in {"s", "skip"}:
            clean_dir(d, args.every_n_epochs)


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
