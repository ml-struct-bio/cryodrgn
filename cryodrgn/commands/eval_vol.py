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
from cryodrgn import config
from cryodrgn.lattice import Lattice
from cryodrgn.models import HetOnlyVAE, load_decoder
from cryodrgn.models_ai import HyperVolume, eval_volume_method as eval_volume_method_ai
from cryodrgn.source import write_mrc
from cryodrgn import utils

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
        "--low-pass",
        type=float,
        help="Low-pass filter resolution in Angstroms (need to specify --Apix)",
    )
    group.add_argument(
        "--crop",
        type=int,
        help="crop volume to this box size after downsampling or low-pass filtering (pixels)",
    )
    group.add_argument(
        "--vol-start-index",
        type=int,
        default=1,
        help="Default value of start index for volume generation (default: %(default)s)",
    )


def check_inputs(args: argparse.Namespace) -> None:
    if args.z_start:
        assert args.z_end, "Must provide --z-end with argument --z-start"
    assert (
        sum((bool(args.z), bool(args.z_start), bool(args.zfile))) == 1
    ), "Must specify either -z OR --z-start/--z-end OR --zfile"


def postprocess_vol(vol, args):
    if args.flip:
        vol = vol.flip([0])
    if args.invert:
        vol *= -1
    if args.low_pass:
        vol = utils.low_pass_filter(vol, args.Apix, args.low_pass)
    if args.crop:
        vol = utils.crop_real_space(vol, args.crop)
    return vol


def reset_origin(oldD, cropD, Apix):
    """Reset origin for cropped volume from (0,0,0) to align with uncropped volume"""
    org = {}
    a = int(oldD / 2 - cropD / 2)
    org["xorg"] = a * Apix
    org["yorg"] = a * Apix
    org["zorg"] = a * Apix
    return org


def main(args: argparse.Namespace) -> None:
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    check_inputs(args)
    t1 = dt.now()

    # Find whether there is a GPU device to compute on and set the device
    if args.device is not None:
        device = torch.device(f"cuda:{args.device}")
    else:
        use_cuda = torch.cuda.is_available()
        device = torch.device("cuda" if use_cuda else "cpu")
        logger.info("Use cuda {}".format(use_cuda))
        if not use_cuda:
            logger.warning("WARNING: No GPUs detected")

    logger.info(args)
    cfg = config.load(args.config)
    logger.info("Loaded configuration:")
    pprint.pprint(cfg)

    D = cfg["lattice_args"]["D"]  # image size + 1
    zdim = cfg["model_args"]["zdim"]
    dataset_args = cfg.get("dataset_args") or {}
    if "norm" in dataset_args:
        norm = [float(x) for x in dataset_args["norm"]]
    elif "data_norm_mean" in cfg and "data_norm_std" in cfg:
        norm = [float(cfg["data_norm_mean"]), float(cfg["data_norm_std"])]
    else:
        raise KeyError(
            "config must include dataset_args['norm'] or "
            "data_norm_mean and data_norm_std (abinit-style configs)"
        )

    if args.downsample:
        if args.downsample % 2 != 0:
            raise ValueError(f"Boxsize {args.downsample} is not even!")
        if args.downsample >= D:
            raise ValueError(
                f"New boxsize {args.downsample} must be "
                f"smaller than original box size {D}!"
            )

    use_abinit_hypervolume = "data_norm_mean" in cfg
    if use_abinit_hypervolume and args.downsample:
        raise ValueError(
            "cryodrgn eval_vol does not support --downsample for abinit "
            "(DrgnAI / HyperVolume) checkpoints"
        )

    # load model
    decoder = None
    hypervolume = None
    zdim_hv = None
    radius_mask = None

    if use_abinit_hypervolume:
        checkpoint = torch.load(args.weights, map_location=device, weights_only=False)
        hypervolume_params = checkpoint["hypervolume_params"]
        hypervolume = HyperVolume(**hypervolume_params)
        hypervolume.load_state_dict(checkpoint["hypervolume_state_dict"])
        hypervolume.eval()
        hypervolume.to(device)
        lattice = Lattice(
            checkpoint["hypervolume_params"]["resolution"],
            extent=0.5,
            device=device,
        )
        zdim_hv = int(checkpoint["hypervolume_params"]["z_dim"])
        if zdim_hv != zdim:
            logger.warning(
                "model_args zdim (%s) != hypervolume z_dim (%s); using config zdim "
                "for z-file layout",
                zdim,
                zdim_hv,
            )
        radius_mask = checkpoint.get("output_mask_radius")
    elif "players" in cfg["model_args"]:  # could be improved
        model, lattice = HetOnlyVAE.load(cfg, args.weights, device=device)
        decoder = model.decoder
        decoder.eval()
    else:
        decoder, lattice = load_decoder(cfg, args.weights, device=device)
        decoder.eval()

    def eval_volume_at_z(zz: np.ndarray):
        if use_abinit_hypervolume:
            assert hypervolume is not None and zdim_hv is not None
            return eval_volume_method_ai(
                hypervolume,
                lattice,
                z_dim=zdim_hv,
                norm=(norm[0], norm[1]),
                zval=zz,
                radius=radius_mask,
            )
        assert decoder is not None
        if args.downsample:
            extent = lattice.extent * (args.downsample / (D - 1))
            return decoder.eval_volume(
                lattice.get_downsample_coords(args.downsample + 1),
                args.downsample + 1,
                extent,
                norm,
                zz,
            )
        return decoder.eval_volume(lattice.coords, lattice.D, lattice.extent, norm, zz)

    # Multiple z
    if args.z_start or args.zfile:
        # Get z values
        if args.z_start:
            args.z_start = np.array(args.z_start)
            args.z_end = np.array(args.z_end)
            z = np.repeat(np.arange(args.n, dtype=np.float32), zdim).reshape(
                (args.n, zdim)
            )
            z *= (args.z_end - args.z_start) / (args.n - 1)  # type: ignore
            z += args.z_start
        else:
            z = np.loadtxt(args.zfile).reshape(-1, zdim)

        os.makedirs(args.o, exist_ok=True)
        logger.info(f"Generating {len(z)} volumes")
        for i, zz in enumerate(z, start=args.vol_start_index):
            logger.info(zz)
            vol = eval_volume_at_z(zz)
            out_mrc = "{}/{}{:03d}.mrc".format(args.o, args.prefix, i)
            org = reset_origin(vol.shape[0], args.crop, args.Apix) if args.crop else {}
            vol = postprocess_vol(vol, args)
            write_mrc(
                out_mrc, np.array(vol.cpu()).astype(np.float32), Apix=args.Apix, **org
            )

    # Single z
    else:
        z = np.array(args.z)
        logger.info(z)
        vol = eval_volume_at_z(z)
        org = reset_origin(vol.shape[0], args.crop, args.Apix) if args.crop else {}
        vol = postprocess_vol(vol, args)
        write_mrc(args.o, np.array(vol.cpu()).astype(np.float32), Apix=args.Apix, **org)

    td = dt.now() - t1
    logger.info(f"Finished in {td}")
