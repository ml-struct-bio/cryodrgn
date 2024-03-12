"""
Heterogeneous NN reconstruction with hierarchical pose optimization
"""
import os
import argparse
import numpy as np
import cryodrgn.utils
from cryodrgn.trainers.hps_trainer import HierarchicalPoseSearchTrainer


def add_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
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
        "--max-threads",
        type=int,
        default=16,
        help="Maximum number of CPU cores for data loading (default: %(default)s)",
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
        default=None,
        help="Choice of beta schedule or a constant for KLD weight (default: 1/zdim)",
    )
    group.add_argument(
        "--beta-control",
        type=float,
        help="KL-Controlled VAE gamma. Beta is KL target",
    )
    group.add_argument(
        "--equivariance",
        type=float,
        help="Strength of equivariance loss",
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
        "--l-ramp-epochs",
        type=int,
        default=0,
        help="Number of epochs to ramp up to --l-end (default: %(default)s)",
    )
    group.add_argument(
        "--l-ramp-model",
        type=int,
        default=0,
        help="If 1, then during ramp only train the model up to l-max (default: %(default)s)",
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
        "--niter",
        type=int,
        default=4,
        help="Number of iterations of grid subdivision (default: %(default)s)",
    )
    group.add_argument(
        "--t-extent",
        type=float,
        default=10,
        help="+/- pixels to search over translations (default: %(default)s)",
    )
    group.add_argument(
        "--t-ngrid",
        type=float,
        default=7,
        help="Initial grid size for translations (default: %(default)s)",
    )
    group.add_argument(
        "--t-xshift",
        type=float,
        default=0,
        help="X-axis translation shift (default: %(default)s)",
    )
    group.add_argument(
        "--t-yshift",
        type=float,
        default=0,
        help="Y-axis translation shift (default: %(default)s)",
    )
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
        help="Number of poses to keep at each refinement interation during branch and bound (default: %(default)s)",
    )
    group.add_argument(
        "--base-healpy",
        type=int,
        default=2,
        help="Base healpy grid for pose search. Higher means exponentially higher resolution (default: %(default)s)",
    )
    group.add_argument(
        "--pose-model-update-freq",
        type=int,
        help="If set, only update the model used for pose search every N examples",
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
        help="Num frequencies in positional encoding (default: image D/2)",
    )
    group.add_argument(
        "--domain",
        choices=("hartley", "fourier"),
        default="hartley",
        help="Volume decoder representation (default: %(default)s)",
    )
    group.add_argument(
        "--activation",
        choices=("relu", "leaky_relu"),
        default="relu",
        help="Activation (default: %(default)s)",
    )
    return parser


def main(args: argparse.Namespace):
    configs = {
        "model": "hps",
        "outdir": args.outdir,
        "particles": args.particles,
        "ctf": args.ctf,
        "pose": None,
        "dataset": None,
        "datadir": args.datadir,
        "ind": args.ind,
        "log_interval": args.log_interval,
        "verbose": args.verbose,
        "load": args.load,
        "load_poses": args.load_poses,
        "checkpoint": args.checkpoint,
        "z_dim": args.zdim,
        "use_gt_poses": False,
        "refine_gt_poses": False,
        "use_gt_trans": False,
        "invert_data": args.invert_data,
        "lazy": args.lazy,
        "window": args.window,
        "window_r": args.window_r,
        "shuffler_size": args.shuffler_size,
        "max_threads": args.max_threads,
        "num_workers": 0,
        "tilt": args.tilt,
        "tilt_deg": args.tilt_deg,
        "num_epochs": args.num_epochs,
        "batch_size": args.batch_size,
        "weight_decay": args.wd,
        "learning_rate": args.lr,
        "pose_learning_rate": args.lr,
        "lattice_extent": 0.5,
        "l_start": args.l_start,
        "l_end": args.l_end,
        "data_norm": args.norm,
        "multigpu": args.multigpu,
        "pretrain": args.pretrain,
        "t_extent": args.t_extent,
        "t_ngrid": args.t_ngrid,
        "t_xshift": args.t_xshift,
        "t_yshift": args.t_yshift,
        "hidden_layers": None,
        "hidden_dim": None,
        "encode_mode": args.encode_mode,
        "enc_mask": args.enc_mask,
        "use_real": args.use_real,
        "pe_type": args.pe_type,
        "pe_dim": args.pe_dim,
        "volume_domain": args.domain,
        "activation": args.activation,
        "feat_sigma": args.feat_sigma,
        "base_healpy": args.base_healpy,
        "subtomo_averaging": False,
        "volume_optim_type": "adam",
        "no_trans": False,
        "amp": False,
        "enc_only": args.enc_only,
        "beta": args.beta,
        "beta_control": args.beta_control,
        "equivariance": args.equivariance,
        "equivariance_start": args.eq_start_it,
        "equivariance_stop": args.eq_end_it,
        "l_ramp_epochs": args.l_ramp_epochs,
        "l_ramp_model": args.l_ramp_model,
        "reset_model_every": args.reset_model_every,
        "reset_optim_every": args.reset_optim_every,
        "reset_optim_after_pretrain": args.reset_optim_after_pretrain,
        "grid_niter": args.niter,
        "ps_freq": args.ps_freq,
        "n_kept_poses": args.nkeptposes,
        "pose_model_update_freq": args.pose_model_update_freq,
        "enc_layers": args.qlayers,
        "enc_dim": args.qdim,
        "dec_layers": args.players,
        "dec_dim": args.pdim,
    }

    cryodrgn.utils._verbose = False
    trainer = HierarchicalPoseSearchTrainer(configs)
    trainer.train()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    main(add_args(parser).parse_args())
