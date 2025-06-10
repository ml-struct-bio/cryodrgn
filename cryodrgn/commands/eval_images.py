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
import cryodrgn.config
from cryodrgn.models.utils import get_model_trainer

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
    group.add_argument("--zdim", type=int, help="Dimension of latent variable")
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

    logger.info(args)
    cfg = cryodrgn.config.overwrite_config(args.config, args)
    logger.info("Loaded configuration:")
    pprint.pprint(cfg)

    model_args = cfg["model_args"] if "model_args" in cfg else cfg
    if "model" not in model_args:
        model_args["model"] = "autoenc"
    cfg["load"] = args.weights
    cfg["shuffle"] = False
    cfg["load_poses"] = args.poses
    cfg["model_args"]["pose_estimation"] = "fixed"
    cfg["dataset_args"]["particles"] = args.particles
    cfg["dataset_args"]["poses"] = args.poses
    if args.ctf is not None:
        cfg["dataset_args"]["ctf"] = args.ctf
    if model_args["model"] != "autoenc" and "encode_mode" in model_args:
        del model_args["encode_mode"]

    # instantiate model
    trainer = get_model_trainer(cfg, outdir=os.path.dirname(args.weights))
    trainer.reconstruction_model.eval()
    z_mu_all = []
    z_logvar_all = []
    gen_loss_accum = 0
    kld_accum = 0
    loss_accum = 0
    batch_it = 0
    trainer.current_epoch = -1
    trainer.use_point_estimates = True

    for minibatch in trainer.data_iterator:
        ind = minibatch["indices"]
        batch_size = len(ind)
        batch_it += batch_size

        losses, tilt_ind, ind, rot, trans, z_mu, z_logvar = trainer.train_batch(
            batch=minibatch
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
                f"# [{batch_it}/{trainer.data.N} images] gen loss={gen_loss:.4f}, "
                f"kld={kld_loss:.4f}, loss={loss:.4f}"
            )
    logger.info(
        f"# =====> Average gen loss = {gen_loss_accum / trainer.data.N:.6}, "
        f"KLD = {kld_accum / trainer.data.N:.6f}, "
        f"total loss = {loss_accum / trainer.data.N:.6f}"
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
                "loss": loss_accum / trainer.data.N,
                "recon": gen_loss_accum / trainer.data.N,
                "kld": kld_accum / trainer.data.N,
            },
            f,
        )

    logger.info("Finished in {}".format(dt.now() - t1))
