"""Creating commands installed with cryoDRGN for use from command line using modules.

This module searches through the `commands` and `commands_utils`
folders for anything that matches the format of a cryoDRGN command module, and creates
a `cryodrgn <x>` command line interface for each such found in the former
and a `cryodrgn_utils <x>` for each found in the latter. This format is kept simple:
anything that is a .py file and has an `add_args` method defined is considered
a command module!

See the `[project.scripts]` entry in the `pyproject.toml` file for how this module
is used to create the commands during installation.

"""
import argparse
import os
from importlib import import_module
import re
import cryodrgn


def _get_commands(cmd_dir: str, doc_str: str = "") -> None:
    """Start up a command line interface using the modules in a directory as subparsers.

    Arguments
    ---------
        cmd_dir: path to folder containing cryoDRGN command modules
        doc_str: short documentation string describing this list of commands as a whole

    """
    parser = argparse.ArgumentParser(description=doc_str)
    parser.add_argument(
        "--version", action="version", version="cryoDRGN " + cryodrgn.__version__
    )

    subparsers = parser.add_subparsers(title="Choose a command")
    subparsers.required = True
    module_files = os.listdir(cmd_dir)
    dir_lbl = os.path.basename(cmd_dir)

    # look for Python modules that have the `add_args` method defined, which is what we
    # use to mark a module in these directories as added to the command namespace
    for module_file in module_files:
        if module_file != "__init__.py" and module_file[-3:] == ".py":
            module_name = ".".join(["cryodrgn", dir_lbl, module_file[:-3]])
            module = import_module(module_name)

            # if this module has the `add_args` method, parse its docstring to get the
            # help message for this command
            if hasattr(module, "add_args"):
                parsed_doc = module.__doc__.split("\n") if module.__doc__ else list()
                descr_txt = parsed_doc[0] if parsed_doc else ""
                epilog_txt = "" if len(parsed_doc) <= 1 else "\n".join(parsed_doc[1:])

                # we have to manually re-add the backslashes used to break up lines
                # for multi-line examples as these get parsed into spaces by .__doc__
                epilog_txt = re.sub(" ([ ]+)", " \\\n\\1", epilog_txt)

                # the docstring header becomes the help message "description", while
                # the rest of the docstring becomes the "epilog"
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
    """Primary commands installed with cryoDRGN as `cryodrgn <cmd_module_name>."""
    _get_commands(
        cmd_dir=os.path.join(os.path.dirname(__file__), "commands"),
        doc_str="Commands installed with cryoDRGN",
    )


def util_commands():
    """Utility commands installed with cryoDRGN as `cryodrgn_utils <cmd_module_name>."""
    _get_commands(
        cmd_dir=os.path.join(os.path.dirname(__file__), "commands_utils"),
        doc_str="Utility commands installed with cryoDRGN",
    )
