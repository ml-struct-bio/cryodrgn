"""
Heterogeneous NN reconstruction with hierarchical pose optimization
"""
import argparse
import os
import pickle
import sys
from datetime import datetime as dt

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader

import cryodrgn
from cryodrgn import ctf, dataset, lie_tools, utils
from cryodrgn.beta_schedule import LinearSchedule, get_beta_schedule
from cryodrgn.lattice import Lattice
from cryodrgn.losses import EquivarianceLoss
from cryodrgn.models import HetOnlyVAE, unparallelize
from cryodrgn.pose_search import PoseSearch

log = utils.log
vlog = utils.vlog


def add_args(parser):
    parser.add_argument(
        "particles", type=os.path.abspath, help="Input particles (.mrcs, .txt or .star)"
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
        "--ctf", metavar="pkl", type=os.path.abspath, help="CTF parameters (.pkl)"
    )
    parser.add_argument("--load", help="Initialize training from a checkpoint")
    parser.add_argument("--load-poses", help="Initialize training from a checkpoint")
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
        "-v", "--verbose", action="store_true", help="Increaes verbosity"
    )
    parser.add_argument(
        "--seed", type=int, default=np.random.randint(0, 100000), help="Random seed"
    )

    group = parser.add_argument_group("Dataset loading")
    group.add_argument(
        "--ind",
        type=os.path.abspath,
        metavar="PKL",
        help="Filter particle stack by these indices",
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
        "--lazy-single",
        action="store_true",
        help="Lazy loading if full dataset is too large to fit in memory",
    )
    group.add_argument(
        "--lazy",
        action="store_true",
        help="Memory efficient training by loading data in chunks",
    )
    group.add_argument(
        "--preprocessed",
        action="store_true",
        help="Skip preprocessing steps if input data is from cryodrgn preprocess_mrcs",
    )
    group.add_argument(
        "--max-threads",
        type=int,
        default=16,
        help="Maximum number of CPU cores for FFT parallelization (default: %(default)s)",
    )

    group = parser.add_argument_group("Tilt series")
    group.add_argument("--tilt", help="Particle stack file (.mrcs)")
    group.add_argument(
        "--tilt-deg",
        type=float,
        default=45,
        help="X-axis tilt offset in degrees (default: %(default)s)",
    )
    group.add_argument(
        "--enc-only",
        action="store_true",
        help="Use the tilt pair only in VAE and not in BNB search",
    )

    group = parser.add_argument_group("Training parameters")
    group.add_argument(
        "-n",
        "--num-epochs",
        type=int,
        default=30,
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
        default=1.0,
        help="Choice of beta schedule or a constant for KLD weight (default: %(default)s)",
    )
    group.add_argument(
        "--beta-control",
        type=float,
        help="KL-Controlled VAE gamma. Beta is KL target. (default: %(default)s)",
    )
    group.add_argument(
        "--equivariance",
        type=float,
        help="Strength of equivariance loss (default: %(default)s)",
    )
    group.add_argument(
        "--eq-start-it",
        type=int,
        default=100000,
        help="It at which equivariance turned on (default: %(default)s)",
    )
    group.add_argument(
        "--eq-end-it",
        type=int,
        default=200000,
        help="It at which equivariance max (default: %(default)s)",
    )
    group.add_argument(
        "--norm",
        type=float,
        nargs=2,
        default=None,
        help="Data normalization as shift, 1/scale (default: mean, std of dataset)",
    )
    group.add_argument(
        "--l-ramp-epochs", type=int, default=0, help="default: %(default)s"
    )
    group.add_argument(
        "--l-ramp-model",
        type=int,
        default=0,
        help="If 1, then during ramp only train the model up to l-max",
    )
    group.add_argument(
        "--reset-model-every", type=int, help="If set, reset the model every N epochs"
    )
    group.add_argument(
        "--reset-optim-every",
        type=int,
        help="If set, reset the optimizer every N epochs",
    )
    group.add_argument(
        "--reset-optim-after-pretrain",
        type=int,
        help="If set, reset the optimizer every N epochs",
    )
    group.add_argument(
        "--multigpu",
        action="store_true",
        help="Parallelize training across all detected GPUs",
    )

    group = parser.add_argument_group("Pose Search parameters")
    group.add_argument(
        "--l-start",
        type=int,
        default=12,
        help="Starting L radius (default: %(default)s)",
    )
    group.add_argument(
        "--l-end", type=int, default=32, help="End L radius (default: %(default)s)"
    )
    group.add_argument(
        "--niter", type=int, default=4, help="Number of iterations of grid subdivision"
    )
    group.add_argument(
        "--t-extent",
        type=float,
        default=10,
        help="+/- pixels to search over translations (default: %(default)s)",
    )
    group.add_argument(
        "--t-ngrid", type=float, default=7, help="Initial grid size for translations"
    )
    group.add_argument("--t-xshift", type=float, default=0)
    group.add_argument("--t-yshift", type=float, default=0)
    group.add_argument(
        "--pretrain",
        type=int,
        default=10000,
        help="Number of initial iterations with random poses (default: %(default)s)",
    )
    group.add_argument(
        "--ps-freq",
        type=int,
        default=5,
        help="Frequency of pose inference (default: every %(default)s epochs)",
    )
    group.add_argument(
        "--nkeptposes",
        type=int,
        default=8,
        help="Number of poses to keep at each refinement interation during branch and bound",
    )
    group.add_argument(
        "--base-healpy",
        type=int,
        default=2,
        help="Base healpy grid for pose search. Higher means exponentially higher resolution.",
    )
    group.add_argument(
        "--pose-model-update-freq",
        type=int,
        help="If set, only update the model used for pose search every N examples.",
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
        default=256,
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
        default=256,
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
        help="Num features in positional encoding (default: image D)",
    )
    group.add_argument(
        "--domain",
        choices=("hartley", "fourier"),
        default="fourier",
        help="Decoder representation domain (default: %(default)s)",
    )
    group.add_argument(
        "--activation",
        choices=("relu", "leaky_relu"),
        default="relu",
        help="Activation (default: %(default)s)",
    )
    return parser


def make_model(args, lattice, enc_mask, in_dim):
    return HetOnlyVAE(
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
        activation={"relu": nn.ReLU, "leaky_relu": nn.LeakyReLU}[args.activation],
        feat_sigma=args.feat_sigma,
    )


def pretrain(model, lattice, optim, minibatch, tilt, zdim):
    y, yt = minibatch
    use_tilt = yt is not None
    B = y.size(0)

    model.train()
    optim.zero_grad()

    rot = lie_tools.random_SO3(B, device=y.device)
    z = torch.randn((B, zdim), device=y.device)

    # reconstruct circle of pixels instead of whole image
    mask = lattice.get_circular_mask(lattice.D // 2)

    def gen_slice(R):
        return unparallelize(model).decode(lattice.coords[mask] @ R, z).view(B, -1)

    y = y.view(B, -1)[:, mask]
    if use_tilt:
        yt = yt.view(B, -1)[:, mask]
        gen_loss = 0.5 * F.mse_loss(gen_slice(rot), y) + 0.5 * F.mse_loss(
            gen_slice(tilt @ rot), yt
        )
    else:
        gen_loss = F.mse_loss(gen_slice(rot), y)

    gen_loss.backward()
    optim.step()
    return gen_loss.item()


def train(
    model,
    lattice,
    ps,
    optim,
    L,
    minibatch,
    beta,
    beta_control=None,
    equivariance=None,
    enc_only=False,
    poses=None,
    ctf_params=None,
):
    y, yt = minibatch
    use_tilt = yt is not None
    use_ctf = ctf_params is not None
    B = y.size(0)
    D = lattice.D

    if use_ctf:
        freqs = lattice.freqs2d.unsqueeze(0).expand(
            B, *lattice.freqs2d.shape
        ) / ctf_params[:, 0].view(B, 1, 1)
        ctf_i = ctf.compute_ctf(freqs, *torch.split(ctf_params[:, 1:], 1, 1)).view(
            B, D, D
        )
    else:
        ctf_i = None

    # TODO: Center image?
    # We do this in pose-supervised train_vae

    # VAE inference of z
    model.train()
    optim.zero_grad()
    input_ = (y, yt) if use_tilt else (y,)
    if use_ctf:
        input_ = (x * ctf_i.sign() for x in input_)  # phase flip by the ctf
    z_mu, z_logvar = unparallelize(model).encode(*input_)
    z = unparallelize(model).reparameterize(z_mu, z_logvar)

    if equivariance is not None:
        lamb, equivariance_loss = equivariance
        eq_loss = equivariance_loss(y, z_mu)

    # pose inference
    if poses is not None:  # use provided poses
        rot = poses[0]
        trans = poses[1]
    else:  # pose search
        model.eval()
        with torch.no_grad():
            rot, trans, _base_pose = ps.opt_theta_trans(
                y,
                z=z,
                images_tilt=None if enc_only else yt,
                ctf_i=ctf_i,
            )
        model.train()

    # reconstruct circle of pixels instead of whole image
    mask = lattice.get_circular_mask(L)

    def gen_slice(R):
        slice_ = model(lattice.coords[mask] @ R, z).view(B, -1)
        if use_ctf:
            slice_ *= ctf_i.view(B, -1)[:, mask]
        return slice_

    def translate(img):
        img = lattice.translate_ht(img, trans.unsqueeze(1), mask)
        return img.view(B, -1)

    y = y.view(B, -1)[:, mask]
    if use_tilt:
        yt = yt.view(B, -1)[:, mask]
    y = translate(y)
    if use_tilt:
        yt = translate(yt)

    if use_tilt:
        gen_loss = 0.5 * F.mse_loss(gen_slice(rot), y) + 0.5 * F.mse_loss(
            gen_slice(bnb.tilt @ rot), yt
        )  # noqa: F821
    else:
        gen_loss = F.mse_loss(gen_slice(rot), y)

    kld = -0.5 * torch.mean(1 + z_logvar - z_mu.pow(2) - z_logvar.exp())
    if torch.isnan(kld):
        log(z_mu[0])
        log(z_logvar[0])
        raise RuntimeError("KLD is nan")

    if beta_control is None:
        loss = gen_loss + beta * kld / mask.sum()
    else:
        loss = gen_loss + beta_control * (beta - kld) ** 2 / mask.sum()

    if equivariance is not None:
        loss += lamb * eq_loss

    loss.backward()

    optim.step()
    save_pose = [rot.detach().cpu().numpy()]
    save_pose.append(trans.detach().cpu().numpy())
    return (
        gen_loss.item(),
        kld.item(),
        loss.item(),
        eq_loss.item() if equivariance else None,
        save_pose,
    )


def eval_z(
    model,
    lattice,
    data,
    batch_size,
    device,
    trans=None,
    use_tilt=False,
    ctf_params=None,
):
    assert not model.training
    z_mu_all = []
    z_logvar_all = []
    data_generator = DataLoader(data, batch_size=batch_size, shuffle=False)
    for minibatch in data_generator:
        ind = minibatch[-1]
        y = minibatch[0].to(device)
        if use_tilt:
            yt = minibatch[1].to(device)
        B = len(ind)
        D = lattice.D
        if ctf_params is not None:
            freqs = lattice.freqs2d.unsqueeze(0).expand(
                B, *lattice.freqs2d.shape
            ) / ctf_params[ind, 0].view(B, 1, 1)
            c = ctf.compute_ctf(freqs, *torch.split(ctf_params[ind, 1:], 1, 1)).view(
                B, D, D
            )
        # if trans is not None:
        #    y = lattice.translate_ht(y.view(B,-1), trans[ind].unsqueeze(1)).view(B,D,D)
        #    if yt is not None: yt = lattice.translate_ht(yt.view(B,-1), trans[ind].unsqueeze(1)).view(B,D,D)
        input_ = (y, yt) if use_tilt else (y,)
        if ctf_params is not None:
            input_ = (x * c.sign() for x in input_)  # phase flip by the ctf
        z_mu, z_logvar = unparallelize(model).encode(*input_)
        z_mu_all.append(z_mu.detach().cpu().numpy())
        z_logvar_all.append(z_logvar.detach().cpu().numpy())
    z_mu_all = np.vstack(z_mu_all)
    z_logvar_all = np.vstack(z_logvar_all)
    return z_mu_all, z_logvar_all


def save_checkpoint(
    model,
    lattice,
    optim,
    epoch,
    norm,
    search_pose,
    z_mu,
    z_logvar,
    out_mrc_dir,
    out_weights,
    out_z,
    out_poses,
):
    """Save model weights, latent encoding z, and decoder volumes"""
    # save model weights
    torch.save(
        {
            "epoch": epoch,
            "model_state_dict": model.state_dict(),
            "optimizer_state_dict": optim.state_dict(),
            "search_pose": search_pose,
        },
        out_weights,
    )
    # save z
    with open(out_z, "wb") as f:
        pickle.dump(z_mu, f)
        pickle.dump(z_logvar, f)
    with open(out_poses, "wb") as f:
        pickle.dump(search_pose, f)


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
    )
    if args.tilt is not None:
        dataset_args["particles_tilt"] = args.tilt
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
    )
    config = dict(
        dataset_args=dataset_args, lattice_args=lattice_args, model_args=model_args
    )
    config["seed"] = args.seed
    with open(out_config, "wb") as f:
        pickle.dump(config, f)
        meta = dict(time=dt.now(), cmd=sys.argv, version=cryodrgn.__version__)
        pickle.dump(meta, f)


def sort_poses(poses):
    ind = [x[0] for x in poses]
    ind = np.concatenate(ind)
    rot = [x[1][0] for x in poses]
    rot = np.concatenate(rot)
    rot = rot[np.argsort(ind)]
    if len(poses[0][1]) == 2:
        trans = [x[1][1] for x in poses]
        trans = np.concatenate(trans)
        trans = trans[np.argsort(ind)]
        return (rot, trans)
    return (rot,)


def get_latest(args, flog):
    flog("Detecting latest checkpoint...")
    weights = [f"{args.outdir}/weights.{i}.pkl" for i in range(args.num_epochs)]
    weights = [f for f in weights if os.path.exists(f)]
    args.load = weights[-1]
    flog(f"Loading {args.load}")
    i = args.load.split(".")[-2]
    args.load_poses = f"{args.outdir}/pose.{i}.pkl"
    assert os.path.exists(args.load_poses)
    flog(f"Loading {args.load_poses}")
    return args


def main(args):
    t1 = dt.now()
    if args.outdir is not None and not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    LOG = f"{args.outdir}/run.log"

    def flog(msg):  # HACK: switch to logging module
        return utils.flog(msg, LOG)

    if args.load == "latest":
        args = get_latest(args, flog)
    flog(" ".join(sys.argv))
    flog(args)

    # set the random seed
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    # set the device
    use_cuda = torch.cuda.is_available()
    device = torch.device("cuda" if use_cuda else "cpu")
    flog("Use cuda {}".format(use_cuda))
    if not use_cuda:
        log("WARNING: No GPUs detected")

    # set beta schedule
    try:
        args.beta = float(args.beta)
    except ValueError:
        assert (
            args.beta_control
        ), "Need to set beta control weight for schedule {}".format(args.beta)
    beta_schedule = get_beta_schedule(args.beta)

    # load index filter
    if args.ind is not None:
        flog("Filtering image dataset with {}".format(args.ind))
        ind = pickle.load(open(args.ind, "rb"))
    else:
        ind = None

    # load dataset
    flog(f"Loading dataset from {args.particles}")
    if args.tilt is None:
        tilt = None
        args.use_real = args.encode_mode == "conv"

        if args.lazy:
            assert (
                args.preprocessed
            ), "Dataset must be preprocesed with `cryodrgn preprocess_mrcs` in order to use --lazy data loading"
            assert (
                not args.ind
            ), "For --lazy data loading, dataset must be filtered by `cryodrgn preprocess_mrcs`"
            # data = dataset.PreprocessedMRCData(args.particles, norm=args.norm)
            raise NotImplementedError("Use --lazy-single for on-the-fly image loading")
        elif args.lazy_single:
            data = dataset.LazyMRCData(
                args.particles,
                norm=args.norm,
                invert_data=args.invert_data,
                ind=ind,
                keepreal=args.use_real,
                window=args.window,
                datadir=args.datadir,
                window_r=args.window_r,
                flog=flog,
            )
        elif args.preprocessed:
            flog(
                "Using preprocessed inputs. Ignoring any --window/--invert-data options"
            )
            data = dataset.PreprocessedMRCData(
                args.particles, norm=args.norm, ind=ind, flog=flog
            )
        else:
            data = dataset.MRCData(
                args.particles,
                norm=args.norm,
                invert_data=args.invert_data,
                ind=ind,
                keepreal=args.use_real,
                window=args.window,
                datadir=args.datadir,
                max_threads=args.max_threads,
                window_r=args.window_r,
                flog=flog,
            )

    # Tilt series data -- lots of unsupported features
    else:
        assert args.encode_mode == "tilt"
        if args.lazy_single:
            raise NotImplementedError
        if args.lazy:
            raise NotImplementedError
        if args.preprocessed:
            raise NotImplementedError
        data = dataset.TiltMRCData(
            args.particles,
            args.tilt,
            norm=args.norm,
            invert_data=args.invert_data,
            ind=ind,
            window=args.window,
            keepreal=args.use_real,
            datadir=args.datadir,
            window_r=args.window_r,
            flog=flog,
        )
        tilt = torch.tensor(utils.xrot(args.tilt_deg).astype(np.float32), device=device)
    Nimg = data.N
    D = data.D

    if args.encode_mode == "conv":
        assert D - 1 == 64, "Image size must be 64x64 for convolutional encoder"

    # load ctf
    if args.ctf is not None:
        flog("Loading ctf params from {}".format(args.ctf))
        ctf_params = ctf.load_ctf_for_training(D - 1, args.ctf)
        if args.ind is not None:
            ctf_params = ctf_params[ind]
        assert ctf_params.shape == (Nimg, 8), ctf_params.shape
        ctf_params = torch.tensor(ctf_params, device=device)
    else:
        ctf_params = None

    lattice = Lattice(D, extent=0.5, device=device)
    if args.enc_mask is None:
        args.enc_mask = D // 2
    if args.enc_mask > 0:
        assert args.enc_mask <= D // 2
        enc_mask = lattice.get_circular_mask(args.enc_mask)
        in_dim = enc_mask.sum()
    elif args.enc_mask == -1:
        enc_mask = None
        in_dim = D**2
    else:
        raise RuntimeError(
            "Invalid argument for encoder mask radius {}".format(args.enc_mask)
        )

    model = make_model(args, lattice, enc_mask, in_dim)
    model.to(device)
    flog(model)
    flog(
        "{} parameters in model".format(
            sum(p.numel() for p in model.parameters() if p.requires_grad)
        )
    )

    if args.equivariance:
        assert args.equivariance > 0, "Regularization weight must be positive"
        equivariance_lambda = LinearSchedule(
            0, args.equivariance, args.eq_start_it, args.eq_end_it
        )
        equivariance_loss = EquivarianceLoss(model, D)

    optim = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.wd)

    if args.load == "latest":
        args = get_latest(args)

    if args.load:
        args.pretrain = 0
        flog("Loading checkpoint from {}".format(args.load))
        checkpoint = torch.load(args.load)
        model.load_state_dict(checkpoint["model_state_dict"])
        optim.load_state_dict(checkpoint["optimizer_state_dict"])
        start_epoch = checkpoint["epoch"] + 1
        model.train()
        if args.load_poses:
            sorted_poses = utils.load_pkl(args.load_poses)
    else:
        start_epoch = 0

    # parallelize
    if args.multigpu and torch.cuda.device_count() > 1:
        log(f"Using {torch.cuda.device_count()} GPUs!")
        args.batch_size *= torch.cuda.device_count()
        log(f"Increasing batch size to {args.batch_size}")
        model = nn.DataParallel(model)
    elif args.multigpu:
        log(
            f"WARNING: --multigpu selected, but {torch.cuda.device_count()} GPUs detected"
        )

    if args.pose_model_update_freq:
        assert not args.multigpu, "TODO"
        pose_model = make_model(args, lattice, enc_mask, in_dim)
        pose_model.to(device)
        pose_model.eval()
    else:
        pose_model = model

    # save configuration
    out_config = "{}/config.pkl".format(args.outdir)
    save_config(args, data, lattice, model, out_config)

    ps = PoseSearch(
        pose_model,
        lattice,
        args.l_start,
        args.l_end,
        tilt,
        t_extent=args.t_extent,
        t_ngrid=args.t_ngrid,
        niter=args.niter,
        nkeptposes=args.nkeptposes,
        base_healpy=args.base_healpy,
        t_xshift=args.t_xshift,
        t_yshift=args.t_yshift,
        device=device,
    )

    data_iterator = DataLoader(data, batch_size=args.batch_size, shuffle=True)

    # pretrain decoder with random poses
    global_it = 0
    flog("Using random poses for {} iterations".format(args.pretrain))
    while global_it < args.pretrain:
        for batch in data_iterator:
            global_it += len(batch[0])
            batch = (
                (batch[0].to(device), None)
                if tilt is None
                else (batch[0].to(device), batch[1].to(device))
            )
            loss = pretrain(model, lattice, optim, batch, tilt=ps.tilt, zdim=args.zdim)
            if global_it % args.log_interval == 0:
                flog(f"[Pretrain Iteration {global_it}] loss={loss:4f}")
            if global_it > args.pretrain:
                break

    # reset model after pretraining
    if args.reset_optim_after_pretrain:
        flog(">> Resetting optim after pretrain")
        optim = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.wd)

    # training loop
    num_epochs = args.num_epochs
    cc = 0
    if args.pose_model_update_freq:
        pose_model.load_state_dict(model.state_dict())
    for epoch in range(start_epoch, num_epochs):
        t2 = dt.now()
        kld_accum = 0
        gen_loss_accum = 0
        loss_accum = 0
        eq_loss_accum = 0
        batch_it = 0
        poses = []

        L_model = lattice.D // 2
        if args.l_ramp_epochs > 0:
            Lramp = args.l_start + int(
                epoch / args.l_ramp_epochs * (args.l_end - args.l_start)
            )
            ps.Lmin = min(Lramp, args.l_start)
            ps.Lmax = min(Lramp, args.l_end)
            if epoch < args.l_ramp_epochs and args.l_ramp_model:
                L_model = ps.Lmax

        if args.reset_model_every and (epoch - 1) % args.reset_model_every == 0:
            flog(">> Resetting model")
            model = make_model(args, lattice, enc_mask, in_dim)

        if args.reset_optim_every and (epoch - 1) % args.reset_optim_every == 0:
            flog(">> Resetting optim")
            optim = torch.optim.Adam(
                model.parameters(), lr=args.lr, weight_decay=args.wd
            )

        if epoch % args.ps_freq != 0:
            flog("Using previous iteration poses")
        for batch in data_iterator:
            ind = batch[-1]
            ind_np = ind.cpu().numpy()
            batch = (
                (batch[0].to(device), None)
                if tilt is None
                else (batch[0].to(device), batch[1].to(device))
            )
            batch_it += len(batch[0])
            global_it = Nimg * epoch + batch_it

            beta = beta_schedule(global_it)
            if args.equivariance:
                lamb = equivariance_lambda(global_it)
                equivariance_tuple = (lamb, equivariance_loss)
            else:
                equivariance_tuple = None

            # train the model
            if epoch % args.ps_freq != 0:
                p = [torch.tensor(x[ind_np], device=device) for x in sorted_poses]
            else:
                p = None

            cc += len(batch[0])
            if args.pose_model_update_freq and cc > args.pose_model_update_freq:
                pose_model.load_state_dict(model.state_dict())
                cc = 0

            ctf_i = ctf_params[ind] if ctf_params is not None else None
            gen_loss, kld, loss, eq_loss, pose = train(
                model,
                lattice,
                ps,
                optim,
                L_model,
                batch,
                beta,
                args.beta_control,
                equivariance_tuple,
                enc_only=args.enc_only,
                poses=p,
                ctf_params=ctf_i,
            )
            # logging
            poses.append((ind.cpu().numpy(), pose))
            kld_accum += kld * len(ind)
            gen_loss_accum += gen_loss * len(ind)
            if args.equivariance:
                eq_loss_accum += eq_loss * len(ind)

            loss_accum += loss * len(ind)
            if batch_it % args.log_interval == 0:
                eq_log = (
                    f"equivariance={eq_loss:.4f}, lambda={lamb:.4f}, "
                    if args.equivariance
                    else ""
                )
                log(
                    f"# [Train Epoch: {epoch+1}/{num_epochs}] [{batch_it}/{Nimg} images] gen loss={gen_loss:.4f}, "
                    f"kld={kld:.4f}, beta={beta:.4f}, {eq_log}loss={loss:.4f}"
                )

        eq_log = (
            "equivariance = {:.4f}, ".format(eq_loss_accum / Nimg)
            if args.equivariance
            else ""
        )
        flog(
            "# =====> Epoch: {} Average gen loss = {:.4}, KLD = {:.4f}, {}total loss = {:.4f}; Finished in {}".format(
                epoch + 1,
                gen_loss_accum / Nimg,
                kld_accum / Nimg,
                eq_log,
                loss_accum / Nimg,
                dt.now() - t2,
            )
        )

        sorted_poses = sort_poses(poses) if poses else None

        # save checkpoint
        if args.checkpoint and epoch % args.checkpoint == 0:
            out_mrc = "{}/reconstruct.{}.mrc".format(args.outdir, epoch)
            out_weights = "{}/weights.{}.pkl".format(args.outdir, epoch)
            out_poses = "{}/pose.{}.pkl".format(args.outdir, epoch)
            out_z = "{}/z.{}.pkl".format(args.outdir, epoch)
            model.eval()
            with torch.no_grad():
                z_mu, z_logvar = eval_z(
                    model,
                    lattice,
                    data,
                    args.batch_size,
                    device,
                    use_tilt=tilt is not None,
                    ctf_params=ctf_params,
                )
                save_checkpoint(
                    model,
                    lattice,
                    optim,
                    epoch,
                    data.norm,
                    sorted_poses,
                    z_mu,
                    z_logvar,
                    out_mrc,
                    out_weights,
                    out_z,
                    out_poses,
                )

    # save model weights and evaluate the model on 3D lattice
    model.eval()
    out_mrc = "{}/reconstruct".format(args.outdir)
    out_weights = "{}/weights.pkl".format(args.outdir)
    out_poses = "{}/pose.pkl".format(args.outdir)
    out_z = "{}/z.pkl".format(args.outdir)
    with torch.no_grad():
        z_mu, z_logvar = eval_z(
            model,
            lattice,
            data,
            args.batch_size,
            device,
            use_tilt=tilt is not None,
            ctf_params=ctf_params,
        )
        save_checkpoint(
            model,
            lattice,
            optim,
            epoch,
            data.norm,
            sorted_poses,
            z_mu,
            z_logvar,
            out_mrc,
            out_weights,
            out_z,
            out_poses,
        )

    td = dt.now() - t1
    flog("Finished in {} ({} per epoch)".format(td, td / (num_epochs - start_epoch)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    utils._verbose = args.verbose
    main(args)
