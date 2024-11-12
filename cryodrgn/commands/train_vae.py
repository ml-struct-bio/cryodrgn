"""Train a variational auto-encoder for heterogeneous reconstruction with known poses.

This is a wrapper for the special case of `commands.train` when using hierarchical pose
search reconstruction with heterogeneous latent conformations (z_dim > 0).

Example usage
-------------
$ cryodrgn train_vae projections.mrcs -o outs/002_trainvae --lr 0.0001 --zdim 10 \
                                      --poses angles.pkl --ctf test_ctf.pkl -n 50

# Restart after already running the same command with some epochs completed
$ cryodrgn train_vae projections.mrcs -o outs/002_trainvae --lr 0.0001 --zdim 10 \
                                      --poses angles.pkl --ctf test_ctf.pkl \
                                      --load latest -n 100

# cryoDRGN-ET tilt series reconstruction
$ cryodrgn train_vae particles_from_M.star --datadir particleseries -o your-outdir \
                                           --ctf ctf.pkl --poses pose.pkl \
                                           --encode-mode tilt --dose-per-tilt 2.93 \
                                           --zdim 8 --num-epochs 50 --beta .025

"""
import argparse
import os
import pickle
import sys
import contextlib
import logging
from typing import Optional
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

try:
    import apex.amp as amp  # type: ignore  # PYR01
except ImportError:
    pass

import cryodrgn
from cryodrgn import ctf, dataset
from cryodrgn.lattice import Lattice
from cryodrgn.models.variational_autoencoder import HetOnlyVAE, unparallelize
import cryodrgn.config

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "particles",
        type=os.path.abspath,
        help="Input particles (.mrcs, .star, .cs, or .txt)",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=os.path.abspath,
        required=True,
        help="Output directory to save model",
    )
    parser.add_argument(
        "--zdim", type=int, required=True, help="Dimension of latent variable"
    )
    parser.add_argument(
        "--poses", type=os.path.abspath, required=True, help="Image poses (.pkl)"
    )
    parser.add_argument(
        "--ctf", metavar="pkl", type=os.path.abspath, help="CTF parameters (.pkl)"
    )
    parser.add_argument(
        "--load", metavar="WEIGHTS.PKL", help="Initialize training from a checkpoint"
    )
    parser.add_argument(
        "--checkpoint",
        type=int,
        default=1,
        help="Checkpointing interval in N_EPOCHS (default: %(default)s)",
    )
    parser.add_argument(
        "--log-interval",
        type=int,
        default=1000,
        help="Logging interval in N_IMGS (default: %(default)s)",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Increase verbosity"
    )
    parser.add_argument(
        "--seed", type=int, default=np.random.randint(0, 100000), help="Random seed"
    )

    group = parser.add_argument_group("Dataset loading")
    group.add_argument(
        "--ind",
        type=os.path.abspath,
        metavar="PKL",
        help="Filter particles by these indices",
    )
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
        "--datadir",
        type=os.path.abspath,
        help="Path prefix to particle stack if loading relative paths from a .star or .cs file",
    )
    group.add_argument(
        "--lazy",
        action="store_true",
        help="Lazy loading if full dataset is too large to fit in memory",
    )
    group.add_argument(
        "--shuffler-size",
        type=int,
        default=0,
        help="If non-zero, will use a data shuffler for faster lazy data loading.",
    )
    group.add_argument(
        "--num-workers",
        type=int,
        default=0,
        help="Number of subprocesses to use as DataLoader workers. If 0, then use the main process for data loading. (default: %(default)s)",
    )
    group.add_argument(
        "--max-threads",
        type=int,
        default=16,
        help="Maximum number of CPU cores for data loading (default: %(default)s)",
    )

    group = parser.add_argument_group("Tilt series parameters")
    group.add_argument(
        "--ntilts",
        type=int,
        default=10,
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
        default=64,
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
        "-d",
        "--dose-per-tilt",
        type=float,
        help="Expected dose per tilt (electrons/A^2 per tilt) (default: %(default)s)",
    )
    group.add_argument(
        "-a",
        "--angle-per-tilt",
        type=float,
        default=3,
        help="Tilt angle increment per tilt in degrees (default: %(default)s)",
    )

    group = parser.add_argument_group("Training parameters")
    group.add_argument(
        "-n",
        "--num-epochs",
        type=int,
        default=20,
        help="Number of training epochs (default: %(default)s)",
    )
    group.add_argument(
        "-b",
        "--batch-size",
        type=int,
        default=8,
        help="Minibatch size (default: %(default)s)",
    )
    group.add_argument(
        "--wd",
        type=float,
        default=0,
        help="Weight decay in Adam optimizer (default: %(default)s)",
    )
    group.add_argument(
        "--lr",
        type=float,
        default=1e-4,
        help="Learning rate in Adam optimizer (default: %(default)s)",
    )
    group.add_argument(
        "--beta",
        default=None,
        help="Choice of beta schedule or a constant for KLD weight (default: 1/zdim)",
    )
    group.add_argument(
        "--beta-control",
        type=float,
        help="KL-Controlled VAE gamma. Beta is KL target",
    )
    group.add_argument(
        "--norm",
        type=float,
        nargs=2,
        default=None,
        help="Data normalization as shift, 1/scale (default: mean, std of dataset)",
    )
    group.add_argument(
        "--no-amp",
        action="store_false",
        dest="amp",
        help="Do not use mixed-precision training for accelerating training",
    )
    group.add_argument(
        "--multigpu",
        action="store_true",
        help="Parallelize training across all detected GPUs",
    )

    group = parser.add_argument_group("Pose SGD")
    group.add_argument(
        "--do-pose-sgd", action="store_true", help="Refine poses with gradient descent"
    )
    group.add_argument(
        "--pretrain",
        type=int,
        default=1,
        help="Number of epochs with fixed poses before pose SGD (default: %(default)s)",
    )
    group.add_argument(
        "--emb-type",
        choices=("s2s2", "quat"),
        default="quat",
        help="SO(3) embedding type for pose SGD (default: %(default)s)",
    )
    group.add_argument(
        "--pose-lr",
        type=float,
        default=3e-4,
        help="Learning rate for pose optimizer (default: %(default)s)",
    )

    group = parser.add_argument_group("Encoder Network")
    group.add_argument(
        "--enc-layers",
        dest="qlayers",
        type=int,
        default=3,
        help="Number of hidden layers (default: %(default)s)",
    )
    group.add_argument(
        "--enc-dim",
        dest="qdim",
        type=int,
        default=1024,
        help="Number of nodes in hidden layers (default: %(default)s)",
    )
    group.add_argument(
        "--encode-mode",
        default="resid",
        choices=("conv", "resid", "mlp", "tilt"),
        help="Type of encoder network (default: %(default)s)",
    )
    group.add_argument(
        "--enc-mask",
        type=int,
        help="Circular mask of image for encoder (default: D/2; -1 for no mask)",
    )
    group.add_argument(
        "--use-real",
        action="store_true",
        help="Use real space image for encoder (for convolutional encoder)",
    )

    group = parser.add_argument_group("Decoder Network")
    group.add_argument(
        "--dec-layers",
        dest="players",
        type=int,
        default=3,
        help="Number of hidden layers (default: %(default)s)",
    )
    group.add_argument(
        "--dec-dim",
        dest="pdim",
        type=int,
        default=1024,
        help="Number of nodes in hidden layers (default: %(default)s)",
    )
    group.add_argument(
        "--pe-type",
        choices=(
            "geom_ft",
            "geom_full",
            "geom_lowf",
            "geom_nohighf",
            "linear_lowf",
            "gaussian",
            "none",
        ),
        default="gaussian",
        help="Type of positional encoding (default: %(default)s)",
    )
    group.add_argument(
        "--feat-sigma",
        type=float,
        default=0.5,
        help="Scale for random Gaussian features (default: %(default)s)",
    )
    group.add_argument(
        "--pe-dim",
        type=int,
        help="Num frequencies in positional encoding (default: image D/2)",
    )
    group.add_argument(
        "--domain",
        choices=("hartley", "fourier"),
        default="fourier",
        help="Volume decoder representation (default: %(default)s)",
    )
    group.add_argument(
        "--activation",
        choices=("relu", "leaky_relu"),
        default="relu",
        help="Activation (default: %(default)s)",
    )


def train_batch(
    model: nn.Module,
    lattice: Lattice,
    y,
    ntilts: Optional[int],
    rot,
    trans,
    optim,
    beta,
    beta_control=None,
    ctf_params=None,
    yr=None,
    use_amp: bool = False,
    scaler=None,
    dose_filters=None,
):
    optim.zero_grad()
    model.train()
    if trans is not None:
        y = preprocess_input(y, lattice, trans)

    # Cast operations to mixed precision if using torch.cuda.amp.GradScaler()
    if scaler is not None:
        amp_mode = torch.cuda.amp.autocast_mode.autocast()
    else:
        amp_mode = contextlib.nullcontext()

    with amp_mode:
        z_mu, z_logvar, z, y_recon, mask = run_batch(
            model, lattice, y, rot, ntilts, ctf_params, yr
        )
        loss, gen_loss, kld = loss_function(
            z_mu,
            z_logvar,
            y,
            ntilts,
            y_recon,
            mask,
            beta,
            beta_control,
            dose_filters,
        )

    if use_amp:
        if scaler is not None:  # torch mixed precision
            scaler.scale(loss).backward()
            scaler.step(optim)
            scaler.update()
        else:  # apex.amp mixed precision
            with amp.scale_loss(loss, optim) as scaled_loss:
                scaled_loss.backward()
            optim.step()
    else:
        loss.backward()
        optim.step()

    return loss.item(), gen_loss.item(), kld.item()


def loss_function(
    z_mu,
    z_logvar,
    y,
    ntilts: Optional[int],
    y_recon,
    mask,
    beta: float,
    beta_control=None,
    dose_filters=None,
):
    # reconstruction error
    B = y.size(0)
    y = y.view(B, -1)[:, mask]
    if dose_filters is not None:
        y_recon = torch.mul(y_recon, dose_filters[:, mask])
    gen_loss = F.mse_loss(y_recon, y)

    # latent loss
    kld = torch.mean(
        -0.5 * torch.sum(1 + z_logvar - z_mu.pow(2) - z_logvar.exp(), dim=1), dim=0
    )
    if torch.isnan(kld):
        logger.info(z_mu[0])
        logger.info(z_logvar[0])
        raise RuntimeError("KLD is nan")

    # total loss
    if beta_control is None:
        loss = gen_loss + beta * kld / mask.sum().float()
    else:
        loss = gen_loss + beta_control * (beta - kld) ** 2 / mask.sum().float()

    return loss, gen_loss, kld


def eval_z(
    model,
    lattice,
    data,
    batch_size,
    device,
    trans=None,
    use_tilt: bool = False,
    ctf_params=None,
    use_real=False,
    shuffler_size=0,
):
    logger.info("Evaluating z")
    assert not model.training
    z_mu_all = []
    z_logvar_all = []
    data_generator = dataset.make_dataloader(
        data, batch_size=batch_size, shuffler_size=shuffler_size, shuffle=False
    )
    for i, minibatch in enumerate(data_generator):
        ind = minibatch[-1]
        y = minibatch[0].to(device)
        D = lattice.D
        if use_tilt:
            y = y.view(-1, D, D)
            ind = minibatch[1].to(device).view(-1)
        B = len(ind)

        c = None
        if ctf_params is not None:
            freqs = lattice.freqs2d.unsqueeze(0).expand(
                B, *lattice.freqs2d.shape
            ) / ctf_params[ind, 0].view(B, 1, 1)
            c = ctf.compute_ctf(freqs, *torch.split(ctf_params[ind, 1:], 1, 1)).view(
                B, D, D
            )
        if trans is not None:
            y = lattice.translate_ht(y.view(B, -1), trans[ind].unsqueeze(1)).view(
                B, D, D
            )

        if use_real:
            input_ = (torch.from_numpy(data.particles_real[ind]).to(device),)
        else:
            input_ = (y,)
        if c is not None:
            assert not use_real, "Not implemented"
            input_ = (x * c.sign() for x in input_)  # phase flip by the ctf
        _model = unparallelize(model)
        assert isinstance(_model, HetOnlyVAE)
        z_mu, z_logvar = _model.encode(*input_)
        z_mu_all.append(z_mu.detach().cpu().numpy())
        z_logvar_all.append(z_logvar.detach().cpu().numpy())
    z_mu_all = np.vstack(z_mu_all)
    z_logvar_all = np.vstack(z_logvar_all)
    return z_mu_all, z_logvar_all


def save_checkpoint(model, optim, epoch, z_mu, z_logvar, out_weights, out_z):
    """Save model weights, latent encoding z, and decoder volumes"""
    # save model weights
    torch.save(
        {
            "epoch": epoch,
            "model_state_dict": unparallelize(model).state_dict(),
            "optimizer_state_dict": optim.state_dict(),
        },
        out_weights,
    )
    # save z
    with open(out_z, "wb") as f:
        pickle.dump(z_mu, f)
        pickle.dump(z_logvar, f)


def save_config(args, dataset, lattice, model, out_config):
    dataset_args = dict(
        particles=args.particles,
        norm=dataset.norm,
        invert_data=args.invert_data,
        ind=args.ind,
        keepreal=args.use_real,
        window=args.window,
        window_r=args.window_r,
        datadir=args.datadir,
        ctf=args.ctf,
        poses=args.poses,
        do_pose_sgd=args.do_pose_sgd,
    )
    if args.encode_mode == "tilt":
        dataset_args["ntilts"] = args.ntilts

    lattice_args = dict(D=lattice.D, extent=lattice.extent, ignore_DC=lattice.ignore_DC)
    model_args = dict(
        qlayers=args.qlayers,
        qdim=args.qdim,
        players=args.players,
        pdim=args.pdim,
        zdim=args.zdim,
        encode_mode=args.encode_mode,
        enc_mask=args.enc_mask,
        pe_type=args.pe_type,
        feat_sigma=args.feat_sigma,
        pe_dim=args.pe_dim,
        domain=args.domain,
        activation=args.activation,
        tilt_params=dict(
            tdim=args.tdim,
            tlayers=args.tlayers,
            t_emb_dim=args.t_emb_dim,
            ntilts=args.ntilts,
        ),
    )
    config = dict(
        dataset_args=dataset_args, lattice_args=lattice_args, model_args=model_args
    )
    config["seed"] = args.seed
    cryodrgn.config.save(config, out_config)


def get_latest(args):
    # assumes args.num_epochs > latest checkpoint
    logger.info("Detecting latest checkpoint...")
    weights = [f"{args.outdir}/weights.{i}.pkl" for i in range(args.num_epochs)]
    weights = [f for f in weights if os.path.exists(f)]
    args.load = weights[-1]
    logger.info(f"Loading {args.load}")
    if args.do_pose_sgd:
        i = args.load.split(".")[-2]
        args.poses = f"{args.outdir}/pose.{i}.pkl"
        assert os.path.exists(args.poses)
        logger.info(f"Loading {args.poses}")
    return args

def main():
    pass
