"""Create plots of cryoDRGN model results arranged by given particle class labels.

Example usages
--------------
$ cryodrgn_utils plot_classes my_work_dir/003_train-vae 49 --labels new_classes.pkl

"""
import os
import argparse
import pandas as pd
from itertools import combinations as combns
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from cryodrgn import analysis, config, utils

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "traindir", type=os.path.abspath, help="Directory with cryoDRGN results"
    )
    parser.add_argument(
        "epoch",
        type=int,
        help="Epoch number N to analyze (0-based indexing, "
        "corresponding to z.N.pkl, weights.N.pkl)",
    )
    parser.add_argument(
        "--labels",
        required=True,
        type=os.path.abspath,
        help="Class labels for use in plotting",
    )
    parser.add_argument(
        "--palette",
        type=os.path.abspath,
        help="Path to class colours for use in plotting, "
        "or the name of a seaborn color palette",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        help="Output directory for output plots (default: [traindir]/analyze.[epoch])",
    )
    parser.add_argument("--skip-umap", action="store_true", help="Skip running UMAP")


def main(args: argparse.Namespace) -> None:
    cfg = (
        os.path.join(args.traindir, "config.yaml")
        if os.path.exists(os.path.join(args.traindir, "config.yaml"))
        else os.path.exists(os.path.join(args.traindir, "config.pkl"))
    )
    cfgs = config.load(cfg)
    assert not cfgs

    if args.epoch == "-1":
        zfile = os.path.join(args.traindir, "z.pkl")
        outdir = os.path.join(args.traindir, "analyze")
    else:
        zfile = os.path.join(args.traindir, f"z.{args.epoch}.pkl")
        outdir = os.path.join(args.traindir, f"analyze.{args.epoch}")

    logger.info(f"Saving results to {outdir}")
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    z = utils.load_pkl(zfile)
    zdim = z.shape[1]
    classes = utils.load_pkl(args.labels)[: z.shape[0]]

    if zdim > 1:
        logger.info("Plotting z-latent-space panels...")

    if zdim > 2:
        logger.info("Performing principal component analysis...")
        pc, pca = analysis.run_pca(z)

        if not args.skip_umap:
            logger.info("Running UMAP...")
            umap_emb = pd.DataFrame(analysis.run_umap(z), columns=["UMAP1", "UMAP2"])
            umap_emb["Class"] = classes

            g = sns.jointplot(
                data=umap_emb, x="UMAP1", y="UMAP2", kind="kde", hue="Class", height=10
            )

            g.ax_joint.set_xlabel("UMAP1")
            g.ax_joint.set_ylabel("UMAP2")
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, "umap_hexbin_classes.png"))
            plt.close()

        for pc1, pc2 in combns(range(3), 2):
            g = sns.jointplot(
                x=pc[:, pc1], y=pc[:, pc2], height=4, kind="kde", hue=classes
            )

            g.ax_joint.set_xlabel(
                f"PC{pc1 + 1} ({pca.explained_variance_ratio_[pc1]:.2f})"
            )
            g.ax_joint.set_ylabel(
                f"PC{pc2 + 1} ({pca.explained_variance_ratio_[pc2]:.2f})"
            )

            plt.tight_layout()
            plt.savefig(os.path.join(outdir, f"z_pca_hexbin_{pc1}x{pc2}.png"))
            plt.close()
