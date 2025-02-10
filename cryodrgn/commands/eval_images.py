"""Evaluate cryoDRGN model latent variables and loss for a stack of images.

Example usage
-------------

$ cryodrgn eval_images hand.mrcs het_weights.pkl --config config.pkl
                        -o output/out_eval_images_losses.pkl
                       --out-z output/out_eval_images_z.pkl
                       --poses hand_rot.pkl --log-interval 1 --verbose

"""
import argparse
import os
import pickle
import pprint
from datetime import datetime as dt
import logging
import numpy as np
import torch
from cryodrgn import config, ctf, dataset
from cryodrgn.models.utils import load_model
from cryodrgn.pose import PoseTracker
from cryodrgn.trainers import AmortizedInferenceTrainer, HierarchicalPoseSearchTrainer

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "particles",
        type=os.path.abspath,
        help="Input particles (.mrcs, .star, .cs, or .txt)",
    )
    parser.add_argument("weights", help="Model weights")
    parser.add_argument("-c", "--config", metavar="PKL", help="CryoDRGN configuration")
    parser.add_argument(
        "-o",
        metavar="PKL",
        type=os.path.abspath,
        required=True,
        help="Output pickle for losses",
    )
    parser.add_argument(
        "--out-z",
        metavar="PKL",
        type=os.path.abspath,
        required=True,
        help="Output pickle for z",
    )
    parser.add_argument(
        "--poses", type=os.path.abspath, required=True, help="Image poses (.pkl)"
    )
    parser.add_argument(
        "--ctf",
        metavar="pkl",
        type=os.path.abspath,
        help="CTF parameters (.pkl) if particle stack is not phase flipped",
    )
    parser.add_argument(
        "--log-interval",
        type=int,
        default=1000,
        help="Logging interval in N_IMGS (default: %(default)s)",
    )
    parser.add_argument(
        "-b",
        "--batch-size",
        type=int,
        default=64,
        help="Minibatch size (default: %(default)s)",
    )
    parser.add_argument("--beta", type=float, help="KLD weight (default: 1/zdim)")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Increaes verbosity"
    )

    group = parser.add_argument_group("Dataset loading")
    group.add_argument(
        "--uninvert-data",
        dest="invert_data",
        action="store_false",
        help="Do not invert data sign",
    )
    group.add_argument(
        "--no-window",
        dest="window",
        action="store_false",
        help="Turn off real space windowing of dataset",
    )
    group.add_argument(
        "--window-r",
        type=float,
        default=0.85,
        help="Windowing radius (default: %(default)s)",
    )
    group.add_argument(
        "--ind", type=os.path.abspath, help="Filter particle stack by these indices"
    )
    group.add_argument(
        "--lazy",
        action="store_true",
        help="Lazy loading if full dataset is too large to fit in memory",
    )
    group.add_argument(
        "--datadir",
        type=os.path.abspath,
        help="Path prefix to particle stack if loading relative paths from a .star or .cs file",
    )
    group.add_argument(
        "--max-threads",
        type=int,
        default=16,
        help="Maximum number of CPU cores for data loading (default: %(default)s)",
    )

    group = parser.add_argument_group("Tilt series paramters")
    group.add_argument(
        "--ntilts",
        type=int,
        help="Number of tilts to encode (default: %(default)s)",
    )
    group.add_argument(
        "--random-tilts",
        action="store_true",
        help="Randomize ordering of tilts series to encoder",
    )
    group.add_argument(
        "--t-emb-dim",
        type=int,
        default=128,
        help="Intermediate embedding dimension (default: %(default)s)",
    )
    group.add_argument(
        "--tlayers",
        type=int,
        default=3,
        help="Number of hidden layers (default: %(default)s)",
    )
    group.add_argument(
        "--tdim",
        type=int,
        default=1024,
        help="Number of nodes in hidden layers (default: %(default)s)",
    )
    group.add_argument(
        "--dose-per-tilt",
        type=float,
        help="Expected dose per tilt (electrons/A^2 per tilt) (default: %(default)s)",
    )
    group.add_argument(
        "--angle-per-tilt",
        type=float,
        default=3,
        help="Tilt angle increment per tilt in degrees (default: %(default)s)",
    )

    group = parser.add_argument_group(
        "Overwrite architecture hyperparameters in config.yaml"
    )
    group.add_argument("--z_dim", type=int, help="Dimension of latent variable")
    group.add_argument(
        "--norm", type=float, nargs=2, help="Data normalization as shift, 1/scale"
    )
    group.add_argument(
        "--enc-layers", dest="qlayers", type=int, help="Number of hidden layers"
    )
    group.add_argument(
        "--enc-dim", dest="qdim", type=int, help="Number of nodes in hidden layers"
    )
    group.add_argument(
        "--encode-mode",
        choices=("conv", "resid", "mlp", "tilt"),
        help="Type of encoder network",
    )
    group.add_argument(
        "--enc-mask", type=int, help="Circular mask of image for encoder"
    )
    group.add_argument(
        "--use-real",
        action="store_true",
        help="Use real space image for encoder (for convolutional encoder)",
    )
    group.add_argument(
        "--dec-layers", dest="players", type=int, help="Number of hidden layers"
    )
    group.add_argument(
        "--dec-dim", dest="pdim", type=int, help="Number of nodes in hidden layers"
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
        "--feat-sigma",
        type=float,
        default=0.5,
        help="Scale for random Gaussian features",
    )
    group.add_argument(
        "--pe-dim",
        type=int,
        help="Num sinusoid features in positional encoding (default: D/2)",
    )
    group.add_argument(
        "--domain", choices=("hartley", "fourier"), help="Decoder representation domain"
    )
    group.add_argument(
        "--activation",
        choices=("relu", "leaky_relu"),
        default="relu",
        help="Activation (default: %(default)s)",
    )


def main(args: argparse.Namespace) -> None:
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    t1 = dt.now()

    # make output directories
    if not os.path.exists(os.path.dirname(args.o)):
        os.makedirs(os.path.dirname(args.o))
    if not os.path.exists(os.path.dirname(args.out_z)):
        os.makedirs(os.path.dirname(args.out_z))

    # set the device
    use_cuda = torch.cuda.is_available()
    device = torch.device("cuda" if use_cuda else "cpu")
    logger.info("Use cuda {}".format(use_cuda))
    if not use_cuda:
        logger.warning("WARNING: No GPUs detected")

    logger.info(args)
    cfg = config.overwrite_config(args.config, args)
    logger.info("Loaded configuration:")
    pprint.pprint(cfg)

    z_dim = (
        cfg["model_args"]["z_dim"]
        if "z_dim" in cfg["model_args"]
        else cfg["model_args"]["zdim"]
    )
    beta = 1.0 / z_dim if args.beta is None else args.beta

    # load the particles
    if args.ind is not None:
        logger.info("Filtering image dataset with {}".format(args.ind))
        ind = pickle.load(open(args.ind, "rb"))
    else:
        ind = None

    if args.encode_mode != "tilt":
        args.use_real = args.encode_mode == "conv"  # Must be False
        data = dataset.ImageDataset(
            mrcfile=args.particles,
            lazy=args.lazy,
            norm=args.norm,
            invert_data=args.invert_data,
            ind=ind,
            keepreal=args.use_real,
            window=args.window,
            datadir=args.datadir,
            window_r=args.window_r,
            max_threads=args.max_threads,
            device=device,
        )
    else:
        assert args.encode_mode == "tilt"
        data = dataset.TiltSeriesData(  # FIXME: maybe combine with above?
            args.particles,
            args.ntilts,
            args.random_tilts,
            norm=args.norm,
            invert_data=args.invert_data,
            ind=ind,
            keepreal=args.use_real,
            window=args.window,
            datadir=args.datadir,
            max_threads=args.max_threads,
            window_r=args.window_r,
            device=device,
            dose_per_tilt=args.dose_per_tilt,
            angle_per_tilt=args.angle_per_tilt,
        )
    if args.encode_mode == "conv":
        assert data.D - 1 == 64, "Image size must be 64x64 for convolutional encoder"

    # instantiate model
    model, lattice, _ = load_model(cfg, args.weights, device=device)
    model.eval()
    z_mu_all = []
    z_logvar_all = []
    gen_loss_accum = 0
    kld_accum = 0
    loss_accum = 0
    batch_it = 0
    data_generator = dataset.make_dataloader(
        data, batch_size=args.batch_size, shuffle=False
    )

    # load ctf
    if args.ctf is not None:
        if args.use_real:
            raise NotImplementedError(
                "Not implemented with real-space encoder!"
                "Use phase-flipped images instead."
            )
        logger.info("Loading ctf params from {}".format(args.ctf))
        ctf_params = ctf.load_ctf_for_training(data.D - 1, args.ctf)
        if args.ind is not None:
            ctf_params = ctf_params[ind]
        ctf_params = torch.tensor(ctf_params, device=device)
    else:
        ctf_params = None

    if "model" not in cfg["model_args"]:
        cfg["model_args"]["model"] = "hps"
        trainer = HierarchicalPoseSearchTrainer.load_from_config(
            cfg, load_data=False, model=model
        )
    elif cfg["model_args"]["model"] == "amort":
        trainer = AmortizedInferenceTrainer.load_from_config(
            cfg, load_data=False, model=model
        )
    elif cfg["model_args"]["model"] == "hps":
        trainer = HierarchicalPoseSearchTrainer.load_from_config(
            cfg, load_data=False, model=model
        )
    else:
        raise ValueError(
            f"Unrecognized cryoDRGN model `{cfg['model_args']['model']}` specified "
            f"in config file `{args.cfg}`!"
        )

    trainer.current_epoch = -1
    trainer.use_point_estimates = True
    trainer.ctf_params = ctf_params
    trainer.pose_tracker = PoseTracker.load(
        args.poses, data.N, data.D, None, ind, device=device
    )
    if trainer.configs.pose_estimation == "abinit":
        trainer.predicted_rots = np.array(trainer.pose_tracker.rots)
        if trainer.pose_tracker.trans is not None:
            trainer.predicted_trans = np.array(trainer.pose_tracker.trans)
        else:
            trainer.predicted_trans = None

    for minibatch in data_generator:
        ind = minibatch["indices"]
        batch_size = len(ind)
        batch_it += batch_size

        losses, tilt_ind, ind, rot, trans, z_mu, z_logvar = trainer.train_batch(
            batch=minibatch, lattice=lattice
        )
        z_mu = z_mu.detach().cpu().numpy()
        if z_logvar is not None:
            z_logvar = z_logvar.detach().cpu().numpy()
        z_mu_all.append(z_mu)
        if z_logvar is not None:
            z_logvar_all.append(z_logvar)

        loss = losses["total"].item()
        gen_loss = losses["gen"]
        kld_loss = losses["kld"].item() if "kld" in losses else 0.0

        # logging
        gen_loss_accum += gen_loss * batch_size
        kld_accum += kld_loss * batch_size
        loss_accum += loss * batch_size
        if batch_it % args.log_interval == 0:
            logger.info(
                "# [{}/{} images] gen loss={:.4f}, kld={:.4f}, beta={:.4f}, loss={:.4f}".format(
                    batch_it, data.N, gen_loss, kld_loss, beta, loss
                )
            )

    logger.info(
        "# =====> Average gen loss = {:.6}, KLD = {:.6f}, total loss = {:.6f}".format(
            gen_loss_accum / data.N, kld_accum / data.N, loss_accum / data.N
        )
    )

    z_mu_all = np.vstack(z_mu_all)
    if z_logvar_all:
        z_logvar_all = np.vstack(z_logvar_all)

    with open(args.out_z, "wb") as f:
        pickle.dump(z_mu_all, f)
        if z_logvar_all is not None:
            pickle.dump(z_logvar_all, f)

    with open(args.o, "wb") as f:
        pickle.dump(
            {
                "loss": loss_accum / data.N,
                "recon": gen_loss_accum / data.N,
                "kld": kld_accum / data.N,
            },
            f,
        )

    logger.info("Finished in {}".format(dt.now() - t1))
