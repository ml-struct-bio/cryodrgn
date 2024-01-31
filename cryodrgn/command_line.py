"""Creating commands installed with cryoDRGN for use from command line.

This module searches through the `commands` and `commands_utils` folders
for anything that matches the format of a cryoDRGN command module
and creates a `cryodrgn <x>` command line interface for each of the
former and a `cryodrgn_utils <x>` for each of the latter.

See the `[project.scripts]` entry in the `pyproject.toml` file for how this module
is used to create the commands during installation.

"""
import argparse
import os
from importlib import import_module
import cryodrgn


def _get_commands(cmd_dir: str) -> None:
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--version", action="version", version="cryoDRGN " + cryodrgn.__version__
    )

    subparsers = parser.add_subparsers(title="Choose a command")
    subparsers.required = True
    module_files = os.listdir(cmd_dir)
    dir_lbl = os.path.basename(cmd_dir)

    for module_file in module_files:
        if module_file != "__init__.py" and module_file[-3:] == ".py":
            module_name = ".".join(["cryodrgn", dir_lbl, module_file[:-3]])
            module = import_module(module_name)

            if hasattr(module, "add_args"):
                this_parser = subparsers.add_parser(
                    module_file[:-3], description=module.__doc__
                )
                module.add_args(this_parser)
                this_parser.set_defaults(func=module.main)

    args = parser.parse_args()
    args.func(args)


def main_commands():
    """Commands installed with cryoDRGN."""
    _get_commands(os.path.join(os.path.dirname(__file__), "commands"))


def util_commands():
    """Utility commands installed with cryoDRGN."""
    _get_commands(os.path.join(os.path.dirname(__file__), "commands_utils"))
