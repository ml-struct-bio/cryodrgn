"""Evaluate the decoder of a heterogeneous model at given z-latent-space co-ordinates.

Example usage
-------------
# This model used the default of zdim=8
$ cryodrgn eval_vol 004_vae128/weights.pkl -c 004_vae128/config.yaml \
                                           -o zero-vol.mrc -z 0 0 0 0 0 0 0 0

# We can instead specify a z-latent-space path instead of a single location
# Here the model was trained using zdim=4
$ cryodrgn eval_vol 004_vae128/weights.pkl -c 004_vae128/config.yaml -o zero-vol.mrc \
                                           --z-start 0 -1 0 0 --z-end 1 1 1 1

"""
import argparse
import os
import pprint
from datetime import datetime as dt
import logging
import numpy as np
import torch
import cryodrgn.config
from cryodrgn.models.utils import load_model
from cryodrgn.source import write_mrc

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("weights", help="Model weights")

    parser.add_argument(
        "-c",
        "--config",
        metavar="YAML",
        required=True,
        help="CryoDRGN config.yaml file",
    )
    parser.add_argument(
        "-o", type=os.path.abspath, required=True, help="Output .mrc or directory"
    )
    parser.add_argument("--device", type=int, help="Optionally specify CUDA device")
    parser.add_argument(
        "--prefix",
        default="vol_",
        help="Prefix when writing out multiple .mrc files (default: %(default)s)",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Increase verbosity"
    )

    group = parser.add_argument_group("Specify z values")
    group.add_argument("-z", type=np.float32, nargs="*", help="Specify one z-value")
    group.add_argument(
        "--z-start", type=np.float32, nargs="*", help="Specify a starting z-value"
    )
    group.add_argument(
        "--z-end", type=np.float32, nargs="*", help="Specify an ending z-value"
    )
    group.add_argument(
        "-n", type=int, default=10, help="Number of structures between [z_start, z_end]"
    )
    group.add_argument("--zfile", help="Text file with z-values to evaluate")

    group = parser.add_argument_group("Volume arguments")
    group.add_argument(
        "--Apix",
        type=float,
        default=1,
        help="Pixel size to add to .mrc header (default: %(default)s A/pix)",
    )
    group.add_argument(
        "--flip", action="store_true", help="Flip handedness of output volume"
    )
    group.add_argument(
        "--invert", action="store_true", help="Invert contrast of output volume"
    )
    group.add_argument(
        "-d",
        "--downsample",
        type=int,
        help="Downsample volumes to this box size (pixels)",
    )
    group.add_argument(
        "--vol-start-index",
        type=int,
        default=0,
        help="Default value of start index for volume generation (default: %(default)s)",
    )

    group = parser.add_argument_group(
        "Overwrite architecture hyperparameters in config.yaml"
    )
    group.add_argument("--norm", nargs=2, type=float)
    group.add_argument("-D", type=int, help="Box size")
    group.add_argument(
        "--enc-layers", dest="qlayers", type=int, help="Number of hidden layers"
    )
    group.add_argument(
        "--enc-dim", dest="qdim", type=int, help="Number of nodes in hidden layers"
    )
    group.add_argument("--zdim", type=int, help="Dimension of latent variable")
    group.add_argument(
        "--encode-mode",
        choices=("conv", "resid", "mlp", "tilt"),
        help="Type of encoder network",
    )
    group.add_argument(
        "--dec-layers", dest="players", type=int, help="Number of hidden layers"
    )
    group.add_argument(
        "--dec-dim", dest="pdim", type=int, help="Number of nodes in hidden layers"
    )
    group.add_argument(
        "--enc-mask", type=int, help="Circular mask radius for image encoder"
    )
    group.add_argument(
        "--pe-type",
        choices=(
            "geom_ft",
            "geom_full",
            "geom_lowf",
            "geom_nohighf",
            "linear_lowf",
            "none",
        ),
        help="Type of positional encoding",
    )
    group.add_argument(
        "--feat-sigma", type=float, help="Scale for random Gaussian features"
    )
    group.add_argument(
        "--pe-dim",
        type=int,
        help="Num sinusoid features in positional encoding (default: D/2)",
    )
    group.add_argument("--domain", choices=("hartley", "fourier"))
    group.add_argument("--l-extent", type=float, help="Coordinate lattice size")
    group.add_argument(
        "--activation",
        choices=("relu", "leaky_relu"),
        default="relu",
        help="Activation (default: %(default)s)",
    )


def check_inputs(args: argparse.Namespace) -> None:
    if args.z_start:
        assert args.z_end, "Must provide --z-end with argument --z-start"
    assert (
        sum((bool(args.z), bool(args.z_start), bool(args.zfile))) == 1
    ), "Must specify either -z OR --z-start/--z-end OR --zfile"


def main(args: argparse.Namespace) -> None:
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    check_inputs(args)
    t1 = dt.now()

    # set the device
    if args.device is not None:
        device = torch.device(f"cuda:{args.device}")
    else:
        use_cuda = torch.cuda.is_available()
        device = torch.device("cuda" if use_cuda else "cpu")
        logger.info("Use cuda {}".format(use_cuda))
        if not use_cuda:
            logger.warning("WARNING: No GPUs detected")

    logger.info(args)
    cfg = cryodrgn.config.overwrite_config(args.config, args)
    logger.info("Loaded configuration:")
    pprint.pprint(cfg)

    D = cfg["lattice_args"]["D"]  # image size + 1
    zdim = cfg["model_args"]["zdim"]
    norm = [float(x) for x in cfg["dataset_args"]["norm"]]

    if args.downsample:
        if args.downsample % 2 != 0:
            raise ValueError(f"Boxsize {args.downsample} is not even!")
        if args.downsample >= D:
            raise ValueError(
                f"New boxsize {args.downsample=} must be "
                f"smaller than original box size {D=}!"
            )

    model, lattice = load_model(cfg, args.weights, device=device)
    model.eval()

    # parse user inputs for location(s) in the latent space
    if args.zfile:
        z = np.loadtxt(args.zfile).reshape(-1, zdim)

    elif args.z_start:
        z_start = np.array(args.z_start)
        z_end = np.array(args.z_end)
        z = np.repeat(np.arange(args.n, dtype=np.float32), zdim).reshape((args.n, zdim))
        z *= (z_end - z_start) / (args.n - 1)  # type: ignore
        z += z_start

    else:
        z = np.array(args.z)

    if args.downsample:
        coords = lattice.get_downsample_coords(args.downsample + 1)
        D = args.downsample + 1
        extent = lattice.extent * (args.downsample / (D - 1))
    else:
        coords = lattice.coords
        D = lattice.D
        extent = lattice.extent

    def transform_volume(vol):
        if args.flip:
            vol = vol.flip([0])
        if args.invert:
            vol *= -1

        return vol

    # multiple latent space co-ordinates
    if len(z.shape) > 1:
        if not os.path.exists(args.o):
            os.makedirs(args.o)

        logger.info(f"Generating {len(z)} volumes")
        for i, zz in enumerate(z, start=args.vol_start_index):
            logger.info(zz)
            volume = transform_volume(model.eval_volume(coords, D, extent, norm, zz))

            out_mrc = os.path.join(args.o, "{}{:03d}.mrc".format(args.prefix, i))
            write_mrc(
                out_mrc, np.array(volume.cpu()).astype(np.float32), Apix=args.Apix
            )

    # single location in latent space
    else:
        z = np.array(args.z)
        logger.info(z)
        volume = transform_volume(model.eval_volume(coords, D, extent, norm, z))
        write_mrc(args.o, np.array(volume).astype(np.float32), Apix=args.Apix)

    td = dt.now() - t1
    logger.info(f"Finished in {td}")
