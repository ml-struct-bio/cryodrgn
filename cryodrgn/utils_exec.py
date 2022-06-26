'''CryoDRGN utilities'''

def main():
    import argparse, os
    parser = argparse.ArgumentParser(description=__doc__)
    import cryodrgn
    parser.add_argument('--version', action='version', version='cryoDRGN '+cryodrgn.__version__)

    import cryodrgn.commands_utils.filter_star

    modules = [
        cryodrgn.commands_utils.filter_star,
        ]

    subparsers = parser.add_subparsers(title='Choose a command')
    subparsers.required = 'True'

    def get_str_name(module):
        return os.path.splitext(os.path.basename(module.__file__))[0]

    for module in modules:
        this_parser = subparsers.add_parser(get_str_name(module), description=module.__doc__)
        module.add_args(this_parser)
        this_parser.set_defaults(func=module.main)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()

