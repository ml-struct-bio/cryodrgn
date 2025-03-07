"""Evaluate the decoder of a heterogeneous model at given z-latent-space co-ordinates.

Example usage
-------------
# This model used the default of z_dim=8
$ cryodrgn eval_vol 004_vae128/weights.pkl -c 004_vae128/config.yaml \
                                           -o zero-vol.mrc -z 0 0 0 0 0 0 0 0

# We can instead specify a z-latent-space path instead of a single location
# Here the model was trained using z_dim=4
$ cryodrgn eval_vol 004_vae128/weights.pkl -c 004_vae128/config.yaml -o zero-vol.mrc \
                                           --z-start 0 -1 0 0 --z-end 1 1 1 1

"""
import argparse
import os
import pprint
from datetime import datetime as dt
import logging
from typing import Optional
import numpy as np
import torch
import cryodrgn.config
from cryodrgn.lattice import Lattice
from cryodrgn.source import write_mrc
from cryodrgn.models.utils import get_model

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    """The command-line arguments for use with the command `cryodrgn eval_vol`."""

    parser.add_argument("weights", help="Model weights")

    parser.add_argument(
        "--output",
        "-o",
        type=os.path.abspath,
        required=True,
        help="Output .mrc or directory",
    )
    parser.add_argument(
        "--config",
        "-c",
        metavar="YAML",
        required=True,
        help="CryoDRGN config.yaml file",
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
    group.add_argument(
        "--z-val", "-z", type=np.float32, nargs="*", help="Specify one z-value"
    )
    group.add_argument(
        "--z-start", type=np.float32, nargs="*", help="Specify a starting z-value"
    )
    group.add_argument(
        "--z-end", type=np.float32, nargs="*", help="Specify an ending z-value"
    )
    group.add_argument(
        "--volume-count",
        "-n",
        type=int,
        default=10,
        help="Number of structures between [z_start, z_end]",
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
    group.add_argument("--z_dim", type=int, help="Dimension of latent variable")
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


class VolumeEvaluator:
    """An engine for generating volumes from a given model."""

    def __init__(
        self,
        weights,
        cfg_data,
        device=None,
        verbose=False,
        apix=None,
        flip=False,
        invert=False,
        downsample=None,
        **architecture_args,
    ):

        # set the device
        if isinstance(device, torch.device):
            self.device = device
        elif device is not None:
            self.device = torch.device(f"cuda:{device}")

        else:
            use_cuda = torch.cuda.is_available()
            self.device = torch.device("cuda" if use_cuda else "cpu")
            logger.info(f"Use cuda {use_cuda}")
            if not use_cuda:
                logger.warning("WARNING: No GPUs detected")

        cfg_data = cryodrgn.config.overwrite_config(
            cfg_data, argparse.Namespace(**architecture_args)
        )
        logger.info("Loaded configuration:")
        cfg_data["load"] = weights
        pprint.pprint(cfg_data)

        if apix is not None:
            self.apix = apix
        elif "apix" in cfg_data["dataset_args"]:
            self.apix = cfg_data["dataset_args"]["apix"]
        else:
            self.apix = 1.0

        model = get_model(cfg_data, outdir=os.path.dirname(weights))
        if downsample:
            if downsample % 2 != 0:
                raise ValueError("Boxsize must be even")
            if downsample > model.lattice.D - 1:
                raise ValueError(
                    "Downsampling size must be smaller than original box size"
                )

            self.coords = model.lattice.get_downsample_coords(downsample + 1)
            self.D = downsample + 1
            self.extent = model.lattice.extent * (downsample / (model.lattice.D - 1))
            self.lattice = Lattice(
                self.D, extent=model.lattice.extent, device=self.device
            )
        else:
            self.lattice = model.lattice
            self.coords = model.lattice.coords
            self.D = model.lattice.D
            self.extent = model.lattice.extent

        self.verbose = verbose
        self.apix = apix
        self.flip = flip
        self.invert = invert
        self.model = model
        self.model.eval()
        self.norm = cfg_data["dataset_args"]["norm"]

    def transform_volume(self, vol):
        if self.flip:
            vol = vol.flip([0])
        if self.invert:
            vol *= -1

        return vol

    def evaluate_volume(self, z):
        return self.transform_volume(
            self.model.eval_volume(
                lattice=self.lattice,
                coords=self.coords,
                resolution=self.D,
                extent=self.extent,
                norm=self.norm,
                zval=z,
                radius=None,
            )
        )

    def produce_volumes(
        self,
        z_values: np.array,
        outpath: str,
        prefix: str = "vol_",
        suffix: Optional[str] = None,
    ) -> None:

        # multiple latent space co-ordinates
        if len(z_values.shape) > 1:
            os.makedirs(outpath, exist_ok=True)

            logger.info(f"Generating {len(z_values)} volumes")
            for i, z_val in enumerate(z_values):
                logger.info(z_val)
                volume = self.evaluate_volume(z_val)
                suffix_str = "" if suffix is None else suffix
                write_mrc(
                    os.path.join(
                        outpath,
                        f"{prefix}{(i+1):03d}{suffix_str}.mrc".format(
                            prefix, i, suffix_str
                        ),
                    ),
                    np.array(volume.cpu()).astype(np.float32),
                    Apix=self.apix,
                )

        # single location in latent space
        else:
            logger.info(z_values)
            volume = self.evaluate_volume(z_values)
            write_mrc(outpath, np.array(volume).astype(np.float32), Apix=self.apix)


def main(args: argparse.Namespace) -> None:
    """Running the command `cryodrgn eval_vol` (see `add_args` above for arguments)."""

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.info(args)
    t0 = dt.now()

    cfg = cryodrgn.config.load(args.config)
    evaluator = VolumeEvaluator(
        args.weights,
        cfg,
        args.device,
        args.verbose,
        args.Apix,
        args.flip,
        args.invert,
        args.downsample,
        **{k: v for k, v in vars(args).items() if k in cfg},
    )

    z_bounds = (args.z_start, args.z_end) if args.z_start is not None else None
    if (args.z_val is None) + (z_bounds is None) + (args.zfile is None) != 2:
        raise ValueError(
            "Must specify either a single z value (-z) "
            "OR z bounds (--z-start AND --z-end) OR z file (--zfile)"
        )

    if z_bounds is not None:
        if len(z_bounds) != 2:
            raise ValueError(
                "`z_bounds` must be given as a list of length two (z-start, z-end)"
            )
        z_start, z_end = z_bounds

        if len(z_start) != len(z_end):
            raise ValueError(
                "`z_bounds` must be given as a list of z-start "
                "and z-end of equal length!"
            )

    # Parse location(s) specified by the user in the latent space in various formats
    if args.zfile:
        z_vals = np.loadtxt(args.zfile).reshape(-1, evaluator.model.z_dim)
    elif args.z_start:
        z_start, z_end = np.array(args.z_start), np.array(args.z_end)
        z_dim = cfg["model_args"]["z_dim"]
        z_vals = np.repeat(
            np.arange(args.volume_count, dtype=np.float32), z_dim
        ).reshape((args.volume_count, z_dim))
        z_vals *= (z_end - z_start) / (args.volume_count - 1)  # type: ignore
        z_vals += z_start
    else:
        z_vals = np.array(args.z_val)

    # Evaluate the volumes at these locations and save them to file
    if len(z_vals):
        evaluator.produce_volumes(z_vals, args.output, args.prefix)
    elif args.zfile:
        logger.warning(f"Given z-values file `{args.zfile}`is empty!")
    elif args.z_start:
        logger.warning("Given z-values range produces empty list of values!")
    else:
        logger.warning("Given z-values list is empty!")

    logger.info(f"Finished in {dt.now() - t0}")
