"""
Train a VAE for heterogeneous reconstruction with known pose
"""
import argparse
import os
import pickle
import sys
import logging
from datetime import datetime as dt
from typing import Optional
import numpy as np
import torch
import torch.nn as nn
from torch.nn.parallel import DataParallel
import torch.nn.functional as F

try:
    import apex.amp as amp  # type: ignore  # PYR01
except ImportError:
    pass

import cryodrgn
from cryodrgn import __version__, ctf, dataset
from cryodrgn.beta_schedule import get_beta_schedule
from cryodrgn.lattice import Lattice
from cryodrgn.models import HetOnlyVAE, unparallelize
from cryodrgn.pose import PoseTracker
import cryodrgn.config

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser):
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
        help="Do not use mixed-precision training",
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
    return parser


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
        with torch.cuda.amp.autocast_mode.autocast():
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
    else:
        # print('AAA', y.shape, rot.shape)
        z_mu, z_logvar, z, y_recon, mask = run_batch(
            model, lattice, y, rot, ntilts, ctf_params, yr
        )
        loss, gen_loss, kld = loss_function(
            z_mu, z_logvar, y, ntilts, y_recon, mask, beta, beta_control, dose_filters
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


def preprocess_input(y, lattice, trans):
    # center the image
    B = y.size(0)
    D = lattice.D
    y = lattice.translate_ht(y.view(B, -1), trans.unsqueeze(1)).view(B, D, D)
    return y


def run_batch(model, lattice, y, rot, ntilts: Optional[int], ctf_params=None, yr=None):
    use_ctf = ctf_params is not None
    B = y.size(0)
    D = lattice.D
    c = None
    if use_ctf:
        freqs = lattice.freqs2d.unsqueeze(0).expand(
            B, *lattice.freqs2d.shape
        ) / ctf_params[:, 0].view(B, 1, 1)
        c = ctf.compute_ctf(freqs, *torch.split(ctf_params[:, 1:], 1, 1)).view(B, D, D)

    # encode
    if yr is not None:
        input_ = (yr,)
    else:
        input_ = (y,)
        if c is not None:
            input_ = (x * c.sign() for x in input_)  # phase flip by the ctf
    _model = unparallelize(model)
    assert isinstance(_model, HetOnlyVAE)
    z_mu, z_logvar = _model.encode(*input_)
    z = _model.reparameterize(z_mu, z_logvar)
    if ntilts is not None:
        z = torch.repeat_interleave(z, ntilts, dim=0)

    # decode
    mask = lattice.get_circular_mask(D // 2)  # restrict to circular mask
    y_recon = model(lattice.coords[mask] / lattice.extent / 2 @ rot, z).view(B, -1)
    if c is not None:
        y_recon *= c.view(B, -1)[:, mask]

    return z_mu, z_logvar, z, y_recon, mask


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
        loss = gen_loss + args.beta_control * (beta - kld) ** 2 / mask.sum().float()
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


def main(args):
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    t1 = dt.now()
    if args.outdir is not None and not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    logger.addHandler(logging.FileHandler(f"{args.outdir}/run.log"))

    if args.load == "latest":
        args = get_latest(args)
    logger.info(" ".join(sys.argv))
    logger.info(f"cryoDRGN {__version__}")
    logger.info(args)

    # set the random seed
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    # set the device
    use_cuda = torch.cuda.is_available()
    device = torch.device("cuda" if use_cuda else "cpu")
    logger.info("Use cuda {}".format(use_cuda))
    if not use_cuda:
        logger.warning("WARNING: No GPUs detected")

    # set beta schedule
    if args.beta is None:
        args.beta = 1.0 / args.ntilts
    try:
        args.beta = float(args.beta)
    except ValueError:
        assert (
            args.beta_control
        ), "Need to set beta control weight for schedule {}".format(args.beta)
    beta_schedule = get_beta_schedule(args.beta)

    # load index filter
    if args.ind is not None:
        logger.info("Filtering image dataset with {}".format(args.ind))
        if args.encode_mode == "tilt":
            particle_ind = pickle.load(open(args.ind, "rb"))
            pt, tp = dataset.TiltSeriesData.parse_particle_tilt(args.particles)
            ind = dataset.TiltSeriesData.particles_to_tilts(pt, particle_ind)
        else:
            ind = pickle.load(open(args.ind, "rb"))
    else:
        ind = None

    # load dataset
    logger.info(f"Loading dataset from {args.particles}")
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
    Nimg = data.N
    D = data.D

    if args.encode_mode == "conv":
        assert D - 1 == 64, "Image size must be 64x64 for convolutional encoder"

    # load poses
    pose_optimizer = None
    if args.do_pose_sgd:
        assert (
            args.domain == "hartley"
        ), "Need to use --domain hartley if doing pose SGD"
    do_pose_sgd = args.do_pose_sgd
    posetracker = PoseTracker.load(
        args.poses, Nimg, D, "s2s2" if do_pose_sgd else None, ind, device=device
    )
    pose_optimizer = (
        torch.optim.SparseAdam(list(posetracker.parameters()), lr=args.pose_lr)
        if do_pose_sgd
        else None
    )

    # load ctf
    if args.ctf is not None:
        if args.use_real:
            raise NotImplementedError(
                "Not implemented with real-space encoder. Use phase-flipped images instead"
            )
        logger.info("Loading ctf params from {}".format(args.ctf))
        ctf_params = ctf.load_ctf_for_training(D - 1, args.ctf)
        if args.ind is not None:
            ctf_params = ctf_params[ind, ...]
        assert ctf_params.shape == (Nimg, 8)
        if args.encode_mode == "tilt":  # TODO: Parse this in cryodrgn parse_ctf_star
            ctf_params = np.concatenate(
                (ctf_params, data.ctfscalefactor.reshape(-1, 1)), axis=1  # type: ignore
            )
            data.voltage = float(ctf_params[0, 4])
        ctf_params = torch.tensor(ctf_params, device=device)  # Nx8
    else:
        ctf_params = None

    # instantiate model
    lattice = Lattice(D, extent=0.5, device=device)
    if args.enc_mask is None:
        args.enc_mask = D // 2
    if args.enc_mask > 0:
        assert args.enc_mask <= D // 2
        enc_mask = lattice.get_circular_mask(args.enc_mask)
        in_dim = int(enc_mask.sum())
    elif args.enc_mask == -1:
        enc_mask = None
        in_dim = lattice.D**2 if not args.use_real else (lattice.D - 1) ** 2
    else:
        raise RuntimeError(
            "Invalid argument for encoder mask radius {}".format(args.enc_mask)
        )
    activation = {"relu": nn.ReLU, "leaky_relu": nn.LeakyReLU}[args.activation]
    tilt_params = {}
    if args.encode_mode == "tilt":
        tilt_params["t_emb_dim"] = args.t_emb_dim
        tilt_params["ntilts"] = args.ntilts
        tilt_params["tlayers"] = args.tlayers
        tilt_params["tdim"] = args.tdim
    model = HetOnlyVAE(
        lattice,
        args.qlayers,
        args.qdim,
        args.players,
        args.pdim,
        in_dim,
        args.zdim,
        encode_mode=args.encode_mode,
        enc_mask=enc_mask,
        enc_type=args.pe_type,
        enc_dim=args.pe_dim,
        domain=args.domain,
        activation=activation,
        feat_sigma=args.feat_sigma,
        tilt_params=tilt_params,
    )
    model.to(device)
    logger.info(model)
    logger.info(
        "{} parameters in model".format(
            sum(p.numel() for p in model.parameters() if p.requires_grad)
        )
    )
    logger.info(
        "{} parameters in encoder".format(
            sum(p.numel() for p in model.encoder.parameters() if p.requires_grad)
        )
    )
    logger.info(
        "{} parameters in decoder".format(
            sum(p.numel() for p in model.decoder.parameters() if p.requires_grad)
        )
    )

    # save configuration
    out_config = "{}/config.yaml".format(args.outdir)
    save_config(args, data, lattice, model, out_config)

    optim = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.wd)

    # Mixed precision training
    scaler = None
    if args.amp:
        assert (
            args.batch_size % 8 == 0
        ), "Batch size must be divisible by 8 for AMP training"
        assert (
            D - 1
        ) % 8 == 0, f"Image size must be divisible by 8 for AMP training: {D}"
        assert (
            args.pdim % 8 == 0
        ), "Decoder hidden layer dimension must be divisible by 8 for AMP training"
        assert (
            args.qdim % 8 == 0
        ), "Encoder hidden layer dimension must be divisible by 8 for AMP training"
        # Also check zdim, enc_mask dim? Add them as warnings for now.
        if args.zdim % 8 != 0:
            logger.warning(
                "Warning: z dimension is not a multiple of 8 -- AMP training speedup is not optimized"
            )
        if in_dim % 8 != 0:
            logger.warning(
                "Warning: Masked input image dimension is not a mutiple of 8 -- AMP training speedup is not optimized"
            )
        try:  # Mixed precision with apex.amp
            model, optim = amp.initialize(model, optim, opt_level="O1")
        except:  # noqa: E722
            # Mixed precision with pytorch (v1.6+)
            scaler = torch.cuda.amp.grad_scaler.GradScaler()

    # restart from checkpoint
    if args.load:
        logger.info("Loading checkpoint from {}".format(args.load))
        checkpoint = torch.load(args.load)
        model.load_state_dict(checkpoint["model_state_dict"])
        optim.load_state_dict(checkpoint["optimizer_state_dict"])
        start_epoch = checkpoint["epoch"] + 1
        model.train()
    else:
        start_epoch = 0

    # parallelize
    num_workers = args.num_workers
    if args.multigpu and torch.cuda.device_count() > 1:
        logger.info(f"Using {torch.cuda.device_count()} GPUs!")
        args.batch_size *= torch.cuda.device_count()
        logger.info(f"Increasing batch size to {args.batch_size}")
        model = DataParallel(model)
    elif args.multigpu:
        logger.warning(
            f"WARNING: --multigpu selected, but {torch.cuda.device_count()} GPUs detected"
        )

    cpu_count = os.cpu_count() or 1
    if num_workers > cpu_count:
        logger.warning(f"Reducing workers to {cpu_count} cpus")
        num_workers = cpu_count

    # training loop
    data_generator = dataset.make_dataloader(
        data,
        batch_size=args.batch_size,
        num_workers=num_workers,
        shuffler_size=args.shuffler_size,
    )

    num_epochs = args.num_epochs
    epoch = None
    Nparticles = Nimg if args.encode_mode != "tilt" else data.Np
    for epoch in range(start_epoch, num_epochs):
        t2 = dt.now()
        gen_loss_accum = 0
        loss_accum = 0
        kld_accum = 0
        batch_it = 0
        for i, minibatch in enumerate(data_generator):  # minibatch: [y, ind]
            ind = minibatch[-1].to(device)
            y = minibatch[0].to(device)
            B = len(ind)
            batch_it += B
            global_it = Nparticles * epoch + batch_it

            beta = beta_schedule(global_it)

            yr = None
            if args.use_real:
                assert hasattr(data, "particles_real")
                yr = torch.from_numpy(data.particles_real[ind.numpy()]).to(device)  # type: ignore  # PYR02
            if pose_optimizer is not None:
                pose_optimizer.zero_grad()

            dose_filters = None
            if args.encode_mode == "tilt":
                tilt_ind = minibatch[1].to(device)
                assert all(tilt_ind >= 0), tilt_ind
                rot, tran = posetracker.get_pose(tilt_ind.view(-1))
                ctf_param = (
                    ctf_params[tilt_ind.view(-1)] if ctf_params is not None else None
                )
                y = y.view(-1, D, D)
                Apix = ctf_params[0, 0] if ctf_params is not None else None
                if args.dose_per_tilt is not None:
                    dose_filters = data.get_dose_filters(tilt_ind, lattice, Apix)
            else:
                rot, tran = posetracker.get_pose(ind)
                ctf_param = ctf_params[ind] if ctf_params is not None else None

            loss, gen_loss, kld = train_batch(
                model,
                lattice,
                y,
                args.ntilts if args.encode_mode == "tilt" else None,
                rot,
                tran,
                optim,
                beta,
                args.beta_control,
                ctf_params=ctf_param,
                yr=yr,
                use_amp=args.amp,
                scaler=scaler,
                dose_filters=dose_filters,
            )
            if pose_optimizer is not None and epoch >= args.pretrain:
                pose_optimizer.step()

            # logging
            gen_loss_accum += gen_loss * B
            kld_accum += kld * B
            loss_accum += loss * B

            if batch_it % args.log_interval == 0:
                logger.info(
                    "# [Train Epoch: {}/{}] [{}/{} particles] gen loss={:.6f}, kld={:.6f}, beta={:.6f}, "
                    "loss={:.6f}".format(
                        epoch + 1,
                        num_epochs,
                        batch_it,
                        Nparticles,
                        gen_loss,
                        kld,
                        beta,
                        loss,
                    )
                )
        logger.info(
            "# =====> Epoch: {} Average gen loss = {:.6}, KLD = {:.6f}, total loss = {:.6f}; Finished in {}".format(
                epoch + 1,
                gen_loss_accum / Nparticles,
                kld_accum / Nparticles,
                loss_accum / Nparticles,
                dt.now() - t2,
            )
        )

        if args.checkpoint and epoch % args.checkpoint == 0:
            out_weights = "{}/weights.{}.pkl".format(args.outdir, epoch)
            out_z = "{}/z.{}.pkl".format(args.outdir, epoch)
            model.eval()
            with torch.no_grad():
                z_mu, z_logvar = eval_z(
                    model,
                    lattice,
                    data,
                    args.batch_size,
                    device,
                    trans=posetracker.trans,
                    use_tilt=args.encode_mode == "tilt",
                    ctf_params=ctf_params,
                    use_real=args.use_real,
                    shuffler_size=args.shuffler_size,
                )
                save_checkpoint(model, optim, epoch, z_mu, z_logvar, out_weights, out_z)
            if args.do_pose_sgd and epoch >= args.pretrain:
                out_pose = "{}/pose.{}.pkl".format(args.outdir, epoch)
                posetracker.save(out_pose)

    logger.info("Training complete")
    # save model weights, latent encoding, and evaluate the model on 3D lattice
    out_weights = "{}/weights.pkl".format(args.outdir)
    out_z = "{}/z.pkl".format(args.outdir)
    model.eval()
    with torch.no_grad():
        z_mu, z_logvar = eval_z(
            model,
            lattice,
            data,
            args.batch_size,
            device,
            posetracker.trans,
            args.encode_mode == "tilt",
            ctf_params,
            args.use_real,
        )
        save_checkpoint(model, optim, epoch, z_mu, z_logvar, out_weights, out_z)

    if args.do_pose_sgd and epoch >= args.pretrain:
        out_pose = "{}/pose.pkl".format(args.outdir)
        posetracker.save(out_pose)
    td = dt.now() - t1
    logger.info(
        "Finished in {} ({} per epoch)".format(td, td / (num_epochs - start_epoch))
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
