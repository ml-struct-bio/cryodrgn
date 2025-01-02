"""Train a NN to model a 3D density map given 2D images with pose assignments.

This command is an interface for the zdim=0 case of the 3D reconstruction method
introduced in cryoDRGN v1. It creates an output directory and config file in the
style of cryoDRGN v4 while using a now-deprecated set of command-line arguments.

Example usage
-------------
$ cryodrgn train_nn projections.mrcs --poses angles.pkl --ctf ctf.pkl \
                                     -o output/train_nn -n 10

# Run with more epochs
$ cryodrgn train_nn projections.star --poses angles.pkl --ctf.pkl \
                                     -o outs/003_train-nn --num-epochs 30 --lr 0.01

# Restart after already running the same command with some epochs completed
$ cryodrgn train_nn projections.star --poses angles.pkl --ctf.pkl \
                                     -o outs/003_train-nn --num-epochs 75 --lr 0.01 \
                                     --load latest

"""
import os
import argparse
import numpy as np
import cryodrgn.utils
from cryodrgn.trainers.hps_trainer import HierarchicalPoseSearchTrainer


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
        "--shuffler-size",
        type=int,
        default=0,
        help="If non-zero, will use a data shuffler for faster lazy data loading.",
    )
    group.add_argument(
        "--datadir",
        type=os.path.abspath,
        help="Path prefix to particle stack if loading relative paths from a .star or .cs file",
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
        default=5,
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
        default=1e-4,
        help="Learning rate for pose optimizer (default: %(default)s)",
    )

    group = parser.add_argument_group("Network Architecture")
    group.add_argument(
        "--layers",
        type=int,
        default=3,
        help="Number of hidden layers (default: %(default)s)",
    )
    group.add_argument(
        "--dim",
        type=int,
        default=1024,
        help="Number of nodes in hidden layers (default: %(default)s)",
    )
    group.add_argument(
        "--l-extent",
        type=float,
        default=0.5,
        help="Coordinate lattice size (if not using positional encoding) (default: %(default)s)",
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
        "--pe-dim",
        type=int,
        help="Num frequencies in positional encoding (default: D/2)",
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
    group.add_argument(
        "--feat-sigma",
        type=float,
        default=0.5,
        help="Scale for random Gaussian features (default: %(default)s)",
    )


def main(args: argparse.Namespace) -> None:
    print(
        "WARNING: "
        "This command is deprecated; use `cryodrgn train` as of cryoDRGN v4.0.0."
    )
    configs = {
        "model": "hps",
        "outdir": args.outdir,
        "particles": args.particles,
        "ctf": args.ctf,
        "poses": args.poses,
        "dataset": None,
        "datadir": args.datadir,
        "ind": args.ind,
        "log_interval": args.log_interval,
        "verbose": args.verbose,
        "load": args.load,
        "checkpoint": args.checkpoint,
        "z_dim": 0,
        "use_gt_poses": True,
        "refine_gt_poses": args.do_pose_sgd,
        "use_gt_trans": False,
        "invert_data": args.invert_data,
        "lazy": args.lazy,
        "window": args.window,
        "window_r": args.window_r,
        "shuffler_size": args.shuffler_size,
        "max_threads": None,
        "num_workers": 0,
        "num_epochs": args.num_epochs,
        "batch_size": args.batch_size,
        "weight_decay": args.wd,
        "learning_rate": args.lr,
        "pose_learning_rate": args.pose_lr,
        "l_extent": args.l_extent,
        "data_norm": args.norm,
        "multigpu": args.multigpu,
        "pretrain": args.pretrain,
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
        "subtomo_averaging": False,
        "volume_optim_type": "adam",
        "no_trans": False,
        "amp": args.amp,
        "tilt_enc_only": False,
        "reset_optim_after_pretrain": False,
        "pose_sgd_emb_type": args.emb_type,
    }

    cryodrgn.utils._verbose = False
    trainer = HierarchicalPoseSearchTrainer(configs)
    trainer.train()
