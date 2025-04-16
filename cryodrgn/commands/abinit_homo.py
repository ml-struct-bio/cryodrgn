"""Homogeneous neural net ab initio reconstruction with hierarchical pose optimization.

This command is an interface for the zdim=0 case of the ab initio pose reconstruction
method introduced in cryoDRGN v2. It creates an output directory and config file in the
style of cryoDRGN v4 while using a now-deprecated set of command-line arguments.

Example usages
--------------
$ cryodrgn abinit_homo particles.256.txt --ctf ctf.pkl --ind chosen-particles.pkl \
                                         -o cryodrn-out/256_abinit-homo

"""
import os
import argparse
import numpy as np
import cryodrgn.utils
from cryodrgn.trainers.hps_trainer import HierarchicalPoseSearchTrainer


def add_args(parser: argparse.ArgumentParser) -> None:
    """The command-line arguments for use with the command `cryodrgn abinit_homo`."""

    inputs_group = parser.add_argument_group("Listing Input Datasets")
    inputs_group.add_argument(
        "particles",
        type=os.path.abspath,
        help="Input particles (.mrcs, .star, .cs, or .txt)",
    )
    inputs_group.add_argument(
        "--ctf", metavar="pkl", type=os.path.abspath, help="CTF parameters (.pkl)"
    )
    inputs_group.add_argument("--ind", help="Filter particle stack by these indices")
    inputs_group.add_argument(
        "--datadir",
        type=os.path.abspath,
        help="Path prefix to particle stack if loading relative paths from a .star or .cs file",
    )
    inputs_group.add_argument("--load", help="Initialize training from a checkpoint")
    inputs_group.add_argument(
        "--load-poses",
        type=os.path.abspath,
        help="Initialize training from a checkpoint",
    )

    outputs_group = parser.add_argument_group("Managing Model Outputs and Logging")
    outputs_group.add_argument(
        "-o",
        "--outdir",
        type=os.path.abspath,
        required=True,
        help="Output directory to save model",
    )
    outputs_group.add_argument(
        "--checkpoint",
        type=int,
        default=1,
        help="Checkpointing interval in N_EPOCHS (default: %(default)s)",
    )
    outputs_group.add_argument(
        "--log-interval",
        type=int,
        default=1000,
        help="Logging interval in N_IMGS (default: %(default)s)",
    )
    outputs_group.add_argument(
        "-v", "--verbose", action="store_true", help="Increase verbosity"
    )

    loading_group = parser.add_argument_group("Managing Dataset Loading and Parsing")
    loading_group.add_argument(
        "--norm",
        type=float,
        nargs=2,
        default=None,
        help="Data normalization as shift, 1/scale (default: mean, std of dataset)",
    )
    loading_group.add_argument(
        "--uninvert-data",
        dest="invert_data",
        action="store_false",
        help="Do not invert data sign",
    )
    loading_group.add_argument(
        "--no-window",
        dest="window",
        action="store_false",
        help="Turn off real space windowing of dataset",
    )
    loading_group.add_argument(
        "--window-r",
        type=float,
        default=0.85,
        help="Windowing radius (default: %(default)s)",
    )
    loading_group.add_argument(
        "--lazy",
        action="store_true",
        help="Lazy loading if full dataset is too large to fit in memory",
    )
    loading_group.add_argument(
        "--shuffler-size",
        type=int,
        default=0,
        help="If non-zero, will use a data shuffler for faster lazy data loading.",
    )

    tilt_group = parser.add_argument_group("Tilt series")
    tilt_group.add_argument("--tilt", help="Particle stack file (.mrcs)")
    tilt_group.add_argument(
        "--tilt-deg",
        type=float,
        default=45,
        help="X-axis tilt offset in degrees (default: %(default)s)",
    )

    training_group = parser.add_argument_group("Reconstruction Training Parameters")
    training_group.add_argument(
        "--t-extent",
        type=float,
        default=10,
        help="+/- pixels to search over translations (default: %(default)s)",
    )
    training_group.add_argument(
        "--t-ngrid",
        type=float,
        default=7,
        help="Initial grid size for translations (default: %(default)s)",
    )
    training_group.add_argument(
        "--t-xshift",
        type=float,
        default=0,
        help="X-axis translation shift (default: %(default)s)",
    )
    training_group.add_argument(
        "--t-yshift",
        type=float,
        default=0,
        help="Y-axis translation shift (default: %(default)s)",
    )
    training_group.add_argument(
        "--no-trans", action="store_true", help="Don't search over translations"
    )
    training_group.add_argument(
        "--pretrain",
        type=int,
        default=10000,
        help="Number of initial iterations with random poses (default: %(default)s)",
    )
    training_group.add_argument(
        "--ps-freq",
        type=int,
        default=5,
        help="Frequency of pose inference (default: every %(default)s epochs)",
    )
    training_group.add_argument(
        "-n",
        "--num-epochs",
        type=int,
        default=30,
        help="Number of training epochs (default: %(default)s)",
    )
    training_group.add_argument(
        "-b",
        "--batch-size",
        type=int,
        default=8,
        help="Minibatch size (default: %(default)s)",
    )
    training_group.add_argument(
        "--wd",
        type=float,
        default=0.0,
        help="Weight decay in Adam optimizer (default: %(default)s)",
    )
    training_group.add_argument(
        "--lr",
        type=float,
        default=1e-4,
        help="Learning rate in Adam optimizer (default: %(default)s)",
    )
    training_group.add_argument(
        "--reset-model-every", type=int, help="If set, reset the model every N epochs"
    )
    training_group.add_argument(
        "--reset-optim-every",
        type=int,
        help="If set, reset the optimizer every N epochs",
    )
    training_group.add_argument(
        "--reset-optim-after-pretrain",
        type=int,
        help="If set, reset the optimizer every N epochs",
    )
    training_group.add_argument(
        "--no-amp",
        action="store_false",
        dest="amp",
        help="Do not use mixed-precision training",
    )
    training_group.add_argument(
        "--seed", type=int, default=np.random.randint(0, 100000), help="Random seed"
    )

    pose_group = parser.add_argument_group("Pose Search Parameters")
    pose_group.add_argument(
        "--l-start",
        type=int,
        default=12,
        help="Starting L radius (default: %(default)s)",
    )
    pose_group.add_argument(
        "--l-end", type=int, default=32, help="End L radius (default: %(default)s)"
    )
    pose_group.add_argument(
        "--niter",
        type=int,
        default=4,
        help="Number of iterations of grid subdivision (default: %(default)s)",
    )
    pose_group.add_argument(
        "--l-ramp-epochs",
        type=int,
        default=25,
        help="Number of epochs to ramp up to --l-end (default: %(default)s)",
    )
    pose_group.add_argument(
        "--probabilistic", action="store_true", help="Use probabilistic bound"
    )
    pose_group.add_argument(
        "--nkeptposes",
        type=int,
        default=8,
        help="Number of poses to keep at each refinement interation during branch and bound (default: %(default)s)",
    )
    pose_group.add_argument(
        "--base-healpy",
        type=int,
        default=2,
        help="Base healpy grid for pose search. Higher means exponentially higher resolution (default: %(default)s)",
    )
    pose_group.add_argument(
        "--pose-model-update-freq",
        type=int,
        help="If set, only update the model used for pose search every N examples",
    )

    archt_group = parser.add_argument_group("Network Architecture Parameters")
    archt_group.add_argument(
        "--layers",
        type=int,
        default=3,
        help="Number of hidden layers (default: %(default)s)",
    )
    archt_group.add_argument(
        "--dim",
        type=int,
        default=256,
        help="Number of nodes in hidden layers (default: %(default)s)",
    )
    archt_group.add_argument(
        "--l-extent",
        type=float,
        default=0.5,
        help="Coordinate lattice size (if not using positional encoding) (default: %(default)s)",
    )
    archt_group.add_argument(
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
    archt_group.add_argument(
        "--pe-dim",
        type=int,
        help="Num frequencies in positional encoding (default: D/2)",
    )
    archt_group.add_argument(
        "--domain",
        choices=("hartley", "fourier"),
        default="hartley",
        help="Volume decoder representation (default: %(default)s)",
    )
    archt_group.add_argument(
        "--activation",
        choices=("relu", "leaky_relu"),
        default="relu",
        help="Activation (default: %(default)s)",
    )
    archt_group.add_argument(
        "--feat-sigma",
        type=float,
        default=0.5,
        help="Scale for random Gaussian features (default: %(default)s)",
    )


def main(args: argparse.Namespace) -> None:
    """Passing command-line arguments as configurations to v3 model training engine."""

    print(
        "WARNING: "
        "This command is deprecated; use `cryodrgn train` as of cryoDRGN v4.0.0."
    )
    configs = {
        "model": "cryodrgn",
        "particles": args.particles,
        "ctf": args.ctf,
        "poses": None,
        "dataset": None,
        "datadir": args.datadir,
        "ind": args.ind,
        "pose_estimation": "abinit",
        "seed": args.seed,
        "log_interval": args.log_interval,
        "verbose": args.verbose,
        "load": args.load,
        "load_poses": args.load_poses,
        "checkpoint": args.checkpoint,
        "z_dim": 0,
        "invert_data": args.invert_data,
        "lazy": args.lazy,
        "window": args.window,
        "window_r": args.window_r,
        "shuffle": args.shuffler_size is not None and args.shuffler_size > 0,
        "shuffler_size": args.shuffler_size,
        "max_threads": 16,
        "num_workers": 0,
        "tilt": args.tilt,
        "tilt_deg": args.tilt_deg,
        "num_epochs": args.num_epochs,
        "batch_size": args.batch_size,
        "weight_decay": args.wd,
        "learning_rate": args.lr,
        "pose_learning_rate": args.lr,
        "l_extent": 0.5,
        "l_start": args.l_start,
        "l_end": args.l_end,
        "data_norm": args.norm,
        "multigpu": False,
        "pretrain": args.pretrain,
        "t_extent": args.t_extent,
        "t_ngrid": args.t_ngrid,
        "t_xshift": args.t_xshift,
        "t_yshift": args.t_yshift,
        "hidden_layers": args.layers,
        "hidden_dim": args.dim,
        "encode_mode": None,
        "enc_mask": None,
        "use_real": False,
        "pe_type": args.pe_type,
        "pe_dim": args.pe_dim,
        "volume_domain": args.domain,
        "activation": args.activation,
        "feat_sigma": args.feat_sigma,
        "base_healpy": args.base_healpy,
        "volume_optim_type": "adam",
        "no_trans": False,
        "amp": args.amp,
        "tilt_enc_only": False,
        "beta": None,
        "beta_control": None,
        "equivariance": None,
        "equivariance_start": None,
        "equivariance_stop": None,
        "l_ramp_epochs": args.l_ramp_epochs,
        "l_ramp_model": None,
        "reset_model_every": args.reset_model_every,
        "reset_optim_every": args.reset_optim_every,
        "reset_optim_after_pretrain": args.reset_optim_after_pretrain,
        "grid_niter": args.niter,
        "ps_freq": args.ps_freq,
        "n_kept_poses": args.nkeptposes,
        "pose_model_update_freq": args.pose_model_update_freq,
    }

    cryodrgn.utils._verbose = False
    trainer = HierarchicalPoseSearchTrainer(configs, args.outdir)
    trainer.train()
