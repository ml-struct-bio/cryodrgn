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


def _get_commands(cmd_dir: str, doc_str: str = "") -> None:
    """Start up a parser using the modules in a cryoDRGN directory as subparsers.

    Arguments
    ---------
        cmd_dir: path to folder containing cryoDRGN command modules
        doc_str: documentation for this list of commands as a whole

    """
    parser = argparse.ArgumentParser(description=doc_str)
    parser.add_argument(
        "--version", action="version", version="cryoDRGN " + cryodrgn.__version__
    )

    subparsers = parser.add_subparsers(title="Choose a command")
    subparsers.required = True
    module_files = os.listdir(cmd_dir)
    dir_lbl = os.path.basename(cmd_dir)

    # look for Python modules that have the "add_args" method defined, which is what we
    # use to mark a module in these directories as added to the command namespace
    for module_file in module_files:
        if module_file != "__init__.py" and module_file[-3:] == ".py":
            module_name = ".".join(["cryodrgn", dir_lbl, module_file[:-3]])
            module = import_module(module_name)

            if hasattr(module, "add_args"):
                parsed_doc = module.__doc__.split("\n") if module.__doc__ else list()
                descr_txt = parsed_doc[0] if parsed_doc else ""
                epilog_txt = "" if len(parsed_doc) <= 1 else "\n".join(parsed_doc[1:])

                # we add documentation text parsed from the module's docstring
                this_parser = subparsers.add_parser(
                    module_file[:-3],
                    description=descr_txt,
                    epilog=epilog_txt,
                    formatter_class=argparse.RawTextHelpFormatter,
                )
                module.add_args(this_parser)
                this_parser.set_defaults(func=module.main)

    args = parser.parse_args()
    args.func(args)


def main_commands():
    """Commands installed with cryoDRGN."""
    _get_commands(
        cmd_dir=os.path.join(os.path.dirname(__file__), "commands"),
        doc_str="Commands installed with cryoDRGN",
    )


def util_commands():
    """Utility commands installed with cryoDRGN."""
    _get_commands(
        cmd_dir=os.path.join(os.path.dirname(__file__), "commands_utils"),
        doc_str="Utility commands installed with cryoDRGN",
    )
