
def main():
    import argparse
    parser = argparse.ArgumentParser()
    import cryodrgn
    parser.add_argument('--version', action='version', version='cryoDRGN '+cryodrgn.__version__)

    import cryodrgn.commands.resize_images
    import cryodrgn.commands.parse_pose
    import cryodrgn.commands.parse_ctf
    import cryodrgn.commands.backproject_nn
    import cryodrgn.commands.backproject_voxel
    import cryodrgn.commands.train_vae
    import cryodrgn.commands.eval_vol
    import cryodrgn.commands.eval_images

    modules = [cryodrgn.commands.resize_images,
        cryodrgn.commands.parse_pose,
        cryodrgn.commands.parse_ctf,
        cryodrgn.commands.backproject_nn,
        cryodrgn.commands.backproject_voxel,
        cryodrgn.commands.train_vae,
        cryodrgn.commands.eval_vol,
        cryodrgn.commands.eval_images]

    subparsers = parser.add_subparsers(title='commands', metavar='<command>')
    subparsers.required = 'True'

    def get_str_name(module):
        return os.path.basename(module.__file__).splitext()[0]

    for module in modules:
        this_parser = subparsers.add_parser(get_str_name(module))
        module.add_args(this_parser)
        this_parser.set_defaults(func=module.main)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()

