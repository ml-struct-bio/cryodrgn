"""CryoDRGN utilities"""


def main():
    import argparse
    import os

    parser = argparse.ArgumentParser(description=__doc__)
    import cryodrgn

    parser.add_argument(
        "--version", action="version", version="cryoDRGN " + cryodrgn.__version__
    )

    import cryodrgn.commands_utils.add_psize
    import cryodrgn.commands_utils.concat_pkls
    import cryodrgn.commands_utils.filter_mrcs
    import cryodrgn.commands_utils.filter_pkl
    import cryodrgn.commands_utils.filter_star
    import cryodrgn.commands_utils.flip_hand
    import cryodrgn.commands_utils.invert_contrast
    import cryodrgn.commands_utils.phase_flip
    import cryodrgn.commands_utils.select_clusters
    import cryodrgn.commands_utils.select_random
    import cryodrgn.commands_utils.translate_mrcs
    import cryodrgn.commands_utils.view_cs_header
    import cryodrgn.commands_utils.view_header
    import cryodrgn.commands_utils.view_mrcs
    import cryodrgn.commands_utils.write_cs
    import cryodrgn.commands_utils.write_star

    modules = [
        cryodrgn.commands_utils.add_psize,
        cryodrgn.commands_utils.concat_pkls,
        cryodrgn.commands_utils.filter_mrcs,
        cryodrgn.commands_utils.filter_pkl,
        cryodrgn.commands_utils.filter_star,
        cryodrgn.commands_utils.flip_hand,
        cryodrgn.commands_utils.invert_contrast,
        cryodrgn.commands_utils.phase_flip,
        cryodrgn.commands_utils.select_clusters,
        cryodrgn.commands_utils.select_random,
        cryodrgn.commands_utils.translate_mrcs,
        cryodrgn.commands_utils.view_cs_header,
        cryodrgn.commands_utils.view_header,
        cryodrgn.commands_utils.view_mrcs,
        cryodrgn.commands_utils.write_star,
        cryodrgn.commands_utils.write_cs,
    ]

    subparsers = parser.add_subparsers(title="Choose a command")
    subparsers.required = True

    def get_str_name(module):
        return os.path.splitext(os.path.basename(module.__file__))[0]

    for module in modules:
        this_parser = subparsers.add_parser(
            get_str_name(module), description=module.__doc__
        )
        module.add_args(this_parser)
        this_parser.set_defaults(func=module.main)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
