"""Generate volumes from a trained decoder for latent space co-ordinates z."""

import argparse
import os
import pprint
from datetime import datetime as dt
import logging
import numpy as np
import torch
import cryodrgn.config
from cryodrgn.mrc import MRCFile
from cryodrgn.models.utils import load_model

logger = logging.getLogger(__name__)


def add_args(parser):
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
    group.add_argument(
        "--vol-start-index",
        type=int,
        default=0,
        help="Default value of start index "
        "for volume generation (default: %(default)s)",
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
        if device is not None:
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
        pprint.pprint(cfg_data)

        orig_d = cfg_data["lattice_args"]["D"]  # image size + 1
        self.zdim = cfg_data["model_args"]["zdim"]
        self.norm = [float(x) for x in cfg_data["dataset_args"]["norm"]]

        if downsample:
            if downsample % 2 != 0:
                raise ValueError("Boxsize must be even")
            if downsample > orig_d - 1:
                raise ValueError(
                    "Downsampling size must be smaller than original box size"
                )

        self.model, self.lattice = load_model(cfg_data, weights, device=self.device)
        self.model.eval()

        if downsample:
            self.coords = self.lattice.get_downsample_coords(downsample + 1)
            self.D = downsample + 1
            self.extent = self.lattice.extent * (downsample / (orig_d - 1))
        else:
            self.coords = self.lattice.coords
            self.D = self.lattice.D
            self.extent = self.lattice.extent

        self.verbose = verbose
        self.apix = apix
        self.flip = flip
        self.invert = invert

    def transform_volume(self, vol):
        if self.flip:
            vol = vol.flip([0])
        if self.invert:
            vol *= -1

        return vol

    def evaluate_volume(self, z):
        return self.transform_volume(
            self.model.eval_volume(self.coords, self.D, self.extent, self.norm, z)
        )

    def produce_volumes(
        self,
        z_values: np.array,
        outpath: str,
        prefix: str = "vol_",
        vol_start_index: int = 0,
    ) -> None:

        if vol_start_index > (len(z_values) - 1):
            raise ValueError(
                f"Cannot use vol-start-index={vol_start_index} with only "
                f"{len(z_values)} latent space co-ordinates!"
            )

        # multiple latent space co-ordinates
        if len(z_values.shape) > 1:
            os.makedirs(outpath, exist_ok=True)

            logger.info(f"Generating {len(z_values)} volumes")
            for i, zz in enumerate(z_values, start=vol_start_index):
                logger.info(zz)
                volume = self.evaluate_volume(zz)

                MRCFile.write(
                    os.path.join(outpath, "{}{:03d}.mrc".format(prefix, i)),
                    np.array(volume.cpu()).astype(np.float32),
                    Apix=self.apix,
                )

        # single location in latent space
        else:
            logger.info(z_values)
            volume = self.evaluate_volume(z_values)
            MRCFile.write(outpath, np.array(volume).astype(np.float32), Apix=self.apix)


def main(args):
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.info(args)
    t0 = dt.now()

    cfg = cryodrgn.config.load(args.config)
    arch_args = {a.dest: getattr(args, a.dest, None) for a in args if a in cfg}

    evaluator = VolumeEvaluator(
        args.weights,
        cfg,
        args.device,
        args.verbose,
        args.apix,
        args.flip,
        args.invert,
        args.downsample,
        **arch_args,
    )

    z_bounds = (args.z_start, args.z_end) if args.z_start is not None else None
    if (args.z_val is None) + (z_bounds is None) + (args.z_file is None) != 2:
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

    # parse user inputs for location(s) in the latent space
    if args.zfile:
        z_vals = np.loadtxt(args.zfile).reshape(-1, evaluator.zdim)

    elif args.z_start:
        z_start = np.array(args.z_start)
        z_end = np.array(args.z_end)

        z_vals = np.repeat(
            np.arange(args.volume_count, dtype=np.float32), evaluator.zdim
        ).reshape((args.volume_count, evaluator.zdim))
        z_vals *= (z_end - z_start) / (args.volume_count - 1)  # type: ignore
        z_vals += z_start

    else:
        z_vals = np.array(args.z_val)

    evaluator.produce_volumes(z_vals, args.output, args.prefix, args.vol_start_index)
    logger.info(f"Finished in {dt.now() - t0}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
