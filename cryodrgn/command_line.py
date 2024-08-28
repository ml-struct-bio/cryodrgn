"""Creating the commands installed with cryoDRGN using the package's modules.

Here we add modules under the `cryodrgn.commands` and `cryodrgn.commands_utils` folders
to the namespace of commands that are installed as part of the cryoDRGN package.
Each module in the former folder thus corresponds to a `cryodrgn <module_name>` command,
while those in the latter folder correspond to a `cryodgn_utils <module_name>` command.

See the `[project.scripts]` entry in the `pyproject.toml` file for how this module
is used to create the commands during installation. We list the modules to use
explicitly for each folder in case the namespace is inadvertantly polluted, and also
since automated scanning for command modules is computationally non-trivial.

"""
import argparse
import os
from importlib import import_module
import re
import cryodrgn


def _get_commands(cmd_dir: str, cmds: list[str], doc_str: str = "") -> None:
    """Start up a command line interface using given modules as subcommands.

    Arguments
    ---------
    cmd_dir:    path to folder containing cryoDRGN command modules
    cmds:       list of commands in the above directory we want to use in the package
    doc_str:    short documentation string describing this list of commands as a whole

    """
    parser = argparse.ArgumentParser(description=doc_str)
    parser.add_argument(
        "--version", action="version", version="cryoDRGN " + cryodrgn.__version__
    )

    subparsers = parser.add_subparsers(title="Choose a command")
    subparsers.required = True
    dir_lbl = os.path.basename(cmd_dir)

    # look for Python modules that have the `add_args` method defined, which is what we
    # use to mark a module in these directories as added to the command namespace
    for cmd in cmds:
        module_name = ".".join(["cryodrgn", dir_lbl, cmd])
        module = import_module(module_name)

        if not hasattr(module, "add_args"):
            raise RuntimeError(
                f"Module `{cmd}` under `{cmd_dir}` does not have the required "
                f"`add_args()` function defined; see other modules under the "
                f"same directory for examples!"
            )

        # Parse the module-level documentation appearing at the top of the file
        parsed_doc = module.__doc__.split("\n") if module.__doc__ else list()
        descr_txt = parsed_doc[0] if parsed_doc else ""
        epilog_txt = "" if len(parsed_doc) <= 1 else "\n".join(parsed_doc[1:])

        # We have to manually re-add the backslashes used to break up lines
        # for multi-line examples as these get parsed into spaces by .__doc__
        # NOTE: This means command docstrings shouldn't otherwise have
        # consecutive spaces!
        epilog_txt = re.sub(" ([ ]+)", " \\\n\\1", epilog_txt)

        # the docstring header becomes the help message "description", while
        # the rest of the docstring becomes the "epilog"
        this_parser = subparsers.add_parser(
            cmd,
            description=descr_txt,
            epilog=epilog_txt,
            formatter_class=argparse.RawTextHelpFormatter,
        )
        module.add_args(this_parser)
        this_parser.set_defaults(func=module.main)

    args = parser.parse_args()
    args.func(args)


def main_commands() -> None:
    """Primary commands installed with cryoDRGN as `cryodrgn <cmd_module_name>."""
    _get_commands(
        cmd_dir=os.path.join(os.path.dirname(__file__), "commands"),
        cmds=[
            "abinit_het",
            "abinit_homo",
            "analyze",
            "analyze_landscape",
            "analyze_landscape_full",
            "backproject_voxel",
            "direct_traversal",
            "downsample",
            "eval_images",
            "eval_vol",
            "filter",
            "graph_traversal",
            "parse_ctf_csparc",
            "parse_ctf_star",
            "parse_pose_csparc",
            "parse_pose_star",
            "pc_traversal",
            "train_nn",
            "train_vae",
            "view_config",
        ],
        doc_str="Commands installed with cryoDRGN",
    )


def util_commands() -> None:
    """Utility commands installed with cryoDRGN as `cryodrgn_utils <cmd_module_name>."""
    _get_commands(
        cmd_dir=os.path.join(os.path.dirname(__file__), "commands_utils"),
        cmds=[
            "add_psize",
            "clean",
            "concat_pkls",
            "filter_mrcs",
            "filter_pkl",
            "filter_star",
            "flip_hand",
            "fsc",
            "gen_mask",
            "invert_contrast",
            "phase_flip",
            "plot_classes",
            "plot_fsc",
            "select_clusters",
            "select_random",
            "translate_mrcs",
            "view_cs_header",
            "view_header",
            "view_mrcs",
            "write_cs",
            "write_star",
        ],
        doc_str="Utility commands installed with cryoDRGN",
    )
