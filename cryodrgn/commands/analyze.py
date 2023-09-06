"""
Visualize latent space and generate volumes
"""

import argparse
import os
import os.path
import shutil
from datetime import datetime as dt
import logging
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import cryodrgn
from cryodrgn import analysis, utils, config

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument(
        "workdir", type=os.path.abspath, help="Directory with cryoDRGN results"
    )
    parser.add_argument(
        "epoch",
        type=int,
        help="Epoch number N to analyze (0-based indexing, corresponding to z.N.pkl, weights.N.pkl)",
    )
    parser.add_argument("--device", type=int, help="Optionally specify CUDA device")
    parser.add_argument(
        "-o",
        "--outdir",
        help="Output directory for analysis results (default: [workdir]/analyze.[epoch])",
    )
    parser.add_argument(
        "--skip-vol", action="store_true", help="Skip generation of volumes"
    )
    parser.add_argument("--skip-umap", action="store_true", help="Skip running UMAP")

    group = parser.add_argument_group("Extra arguments for volume generation")
    group.add_argument(
        "--Apix",
        type=float,
        default=1,
        help="Pixel size to add to .mrc header (default: %(default)s A/pix)",
    )
    group.add_argument(
        "--flip", action="store_true", help="Flip handedness of output volumes"
    )
    group.add_argument(
        "--invert", action="store_true", help="Invert contrast of output volumes"
    )
    group.add_argument(
        "-d",
        "--downsample",
        type=int,
        help="Downsample volumes to this box size (pixels)",
    )
    group.add_argument(
        "--pc",
        type=int,
        default=2,
        help="Number of principal component traversals to generate (default: %(default)s)",
    )
    group.add_argument(
        "--ksample",
        type=int,
        default=20,
        help="Number of kmeans samples to generate (default: %(default)s)",
    )
    group.add_argument(
        "--vol-start-index",
        type=int,
        default=0,
        help="Default value of start index for volume generation (default: %(default)s)",
    )

    return parser


def analyze_z1(z, outdir, vg):
    """Plotting and volume generation for 1D z"""
    assert z.shape[1] == 1
    z = z.reshape(-1)
    N = len(z)

    plt.figure(1)
    plt.scatter(np.arange(N), z, alpha=0.1, s=2)
    plt.xlabel("particle")
    plt.ylabel("z")
    plt.savefig(f"{outdir}/z.png")

    plt.figure(2)
    sns.distplot(z)
    plt.xlabel("z")
    plt.savefig(f"{outdir}/z_hist.png")

    ztraj = np.percentile(z, np.linspace(5, 95, 10))
    vg.gen_volumes(outdir, ztraj)


def analyze_zN(z, outdir, vg, skip_umap=False, num_pcs=2, num_ksamples=20):
    zdim = z.shape[1]

    # Principal component analysis
    logger.info("Performing principal component analysis...")
    pc, pca = analysis.run_pca(z)
    logger.info("Generating volumes...")
    for i in range(num_pcs):
        start, end = np.percentile(pc[:, i], (5, 95))
        z_pc = analysis.get_pc_traj(pca, z.shape[1], 10, i + 1, start, end)
        vg.gen_volumes(f"{outdir}/pc{i+1}", z_pc)

    # kmeans clustering
    logger.info("K-means clustering...")
    K = num_ksamples
    kmeans_labels, centers = analysis.cluster_kmeans(z, K)
    centers, centers_ind = analysis.get_nearest_point(z, centers)
    if not os.path.exists(f"{outdir}/kmeans{K}"):
        os.mkdir(f"{outdir}/kmeans{K}")
    utils.save_pkl(kmeans_labels, f"{outdir}/kmeans{K}/labels.pkl")
    np.savetxt(f"{outdir}/kmeans{K}/centers.txt", centers)
    np.savetxt(f"{outdir}/kmeans{K}/centers_ind.txt", centers_ind, fmt="%d")
    logger.info("Generating volumes...")
    vg.gen_volumes(f"{outdir}/kmeans{K}", centers)

    # UMAP -- slow step
    umap_emb = None
    if zdim > 2 and not skip_umap:
        logger.info("Running UMAP...")
        umap_emb = analysis.run_umap(z)
        utils.save_pkl(umap_emb, f"{outdir}/umap.pkl")

    # Make some plots
    logger.info("Generating plots...")

    def plt_pc_labels(x=0, y=1):
        plt.xlabel(f"PC{x+1} ({pca.explained_variance_ratio_[x]:.2f})")
        plt.ylabel(f"PC{y+1} ({pca.explained_variance_ratio_[y]:.2f})")

    def plt_pc_labels_jointplot(g, x=0, y=1):
        g.ax_joint.set_xlabel(f"PC{x+1} ({pca.explained_variance_ratio_[x]:.2f})")
        g.ax_joint.set_ylabel(f"PC{y+1} ({pca.explained_variance_ratio_[y]:.2f})")

    def plt_umap_labels():
        plt.xticks([])
        plt.yticks([])
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")

    def plt_umap_labels_jointplot(g):
        g.ax_joint.set_xlabel("UMAP1")
        g.ax_joint.set_ylabel("UMAP2")

    # PCA -- Style 1 -- Scatter
    plt.figure(figsize=(4, 4))
    plt.scatter(pc[:, 0], pc[:, 1], alpha=0.1, s=1, rasterized=True)
    plt_pc_labels()
    plt.tight_layout()
    plt.savefig(f"{outdir}/z_pca.png")

    # PCA -- Style 2 -- Scatter, with marginals
    g = sns.jointplot(x=pc[:, 0], y=pc[:, 1], alpha=0.1, s=1, rasterized=True, height=4)
    plt_pc_labels_jointplot(g)
    plt.tight_layout()
    plt.savefig(f"{outdir}/z_pca_marginals.png")

    # PCA -- Style 3 -- Hexbin
    g = sns.jointplot(x=pc[:, 0], y=pc[:, 1], height=4, kind="hex")
    plt_pc_labels_jointplot(g)
    plt.tight_layout()
    plt.savefig(f"{outdir}/z_pca_hexbin.png")

    if umap_emb is not None:
        # Style 1 -- Scatter
        plt.figure(figsize=(4, 4))
        plt.scatter(umap_emb[:, 0], umap_emb[:, 1], alpha=0.1, s=1, rasterized=True)
        plt_umap_labels()
        plt.tight_layout()
        plt.savefig(f"{outdir}/umap.png")

        # Style 2 -- Scatter with marginal distributions
        g = sns.jointplot(
            x=umap_emb[:, 0],
            y=umap_emb[:, 1],
            alpha=0.1,
            s=1,
            rasterized=True,
            height=4,
        )
        plt_umap_labels_jointplot(g)
        plt.tight_layout()
        plt.savefig(f"{outdir}/umap_marginals.png")

        # Style 3 -- Hexbin / heatmap
        g = sns.jointplot(x=umap_emb[:, 0], y=umap_emb[:, 1], kind="hex", height=4)
        plt_umap_labels_jointplot(g)
        plt.tight_layout()
        plt.savefig(f"{outdir}/umap_hexbin.png")

    # Plot kmeans sample points
    colors = analysis._get_chimerax_colors(K)
    analysis.scatter_annotate(
        pc[:, 0],
        pc[:, 1],
        centers_ind=centers_ind,
        annotate=True,
        colors=colors,
    )
    plt_pc_labels()
    plt.tight_layout()
    plt.savefig(f"{outdir}/kmeans{K}/z_pca.png")

    g = analysis.scatter_annotate_hex(
        pc[:, 0],
        pc[:, 1],
        centers_ind=centers_ind,
        annotate=True,
        colors=colors,
    )
    plt_pc_labels_jointplot(g)
    plt.tight_layout()
    plt.savefig(f"{outdir}/kmeans{K}/z_pca_hex.png")

    if umap_emb is not None:
        analysis.scatter_annotate(
            umap_emb[:, 0],
            umap_emb[:, 1],
            centers_ind=centers_ind,
            annotate=True,
            colors=colors,
        )
        plt_umap_labels()
        plt.tight_layout()
        plt.savefig(f"{outdir}/kmeans{K}/umap.png")

        g = analysis.scatter_annotate_hex(
            umap_emb[:, 0],
            umap_emb[:, 1],
            centers_ind=centers_ind,
            annotate=True,
            colors=colors,
        )
        plt_umap_labels_jointplot(g)
        plt.tight_layout()
        plt.savefig(f"{outdir}/kmeans{K}/umap_hex.png")

    # Plot PC trajectories
    for i in range(num_pcs):
        start, end = np.percentile(pc[:, i], (5, 95))
        z_pc = analysis.get_pc_traj(pca, z.shape[1], 10, i + 1, start, end)
        if umap_emb is not None:
            # UMAP, colored by PCX
            analysis.scatter_color(
                umap_emb[:, 0],
                umap_emb[:, 1],
                pc[:, i],
                label=f"PC{i+1}",
            )
            plt_umap_labels()
            plt.tight_layout()
            plt.savefig(f"{outdir}/pc{i+1}/umap.png")

            # UMAP, with PC traversal
            z_pc_on_data, pc_ind = analysis.get_nearest_point(z, z_pc)
            dists = ((z_pc_on_data - z_pc) ** 2).sum(axis=1) ** 0.5
            if np.any(dists > 2):
                logger.warn(
                    f"Warning: PC{i+1} point locations in UMAP plot may be inaccurate"
                )
            plt.figure(figsize=(4, 4))
            plt.scatter(
                umap_emb[:, 0], umap_emb[:, 1], alpha=0.05, s=1, rasterized=True
            )
            plt.scatter(
                umap_emb[pc_ind, 0],
                umap_emb[pc_ind, 1],
                c="cornflowerblue",
                edgecolor="black",
            )
            plt_umap_labels()
            plt.tight_layout()
            plt.savefig(f"{outdir}/pc{i+1}/umap_traversal.png")

            # UMAP, with PC traversal, connected
            plt.figure(figsize=(4, 4))
            plt.scatter(
                umap_emb[:, 0], umap_emb[:, 1], alpha=0.05, s=1, rasterized=True
            )
            plt.plot(umap_emb[pc_ind, 0], umap_emb[pc_ind, 1], "--", c="k")
            plt.scatter(
                umap_emb[pc_ind, 0],
                umap_emb[pc_ind, 1],
                c="cornflowerblue",
                edgecolor="black",
            )
            plt_umap_labels()
            plt.tight_layout()
            plt.savefig(f"{outdir}/pc{i+1}/umap_traversal_connected.png")

        # 10 points, from 5th to 95th percentile of PC1 values
        t = np.linspace(start, end, 10, endpoint=True)
        plt.figure(figsize=(4, 4))
        if i > 0 and i == num_pcs - 1:
            plt.scatter(pc[:, i - 1], pc[:, i], alpha=0.1, s=1, rasterized=True)
            plt.scatter(np.zeros(10), t, c="cornflowerblue", edgecolor="white")
            plt_pc_labels(i - 1, i)
        else:
            plt.scatter(pc[:, i], pc[:, i + 1], alpha=0.1, s=1, rasterized=True)
            plt.scatter(t, np.zeros(10), c="cornflowerblue", edgecolor="white")
            plt_pc_labels(i, i + 1)
        plt.tight_layout()
        plt.savefig(f"{outdir}/pc{i+1}/pca_traversal.png")

        if i > 0 and i == num_pcs - 1:
            g = sns.jointplot(
                pc[:, i - 1], pc[:, i], alpha=0.1, s=1, rasterized=True, height=4
            )
            g.ax_joint.scatter(np.zeros(10), t, c="cornflowerblue", edgecolor="white")
            plt_pc_labels_jointplot(g, i - 1, i)
        else:
            g = sns.jointplot(
                x=pc[:, i], y=pc[:, i + 1], alpha=0.1, s=1, rasterized=True, height=4
            )
            g.ax_joint.scatter(t, np.zeros(10), c="cornflowerblue", edgecolor="white")
            plt_pc_labels_jointplot(g)
        plt.tight_layout()
        plt.savefig(f"{outdir}/pc{i+1}/pca_traversal_hex.png")


class VolumeGenerator:
    """Helper class to call analysis.gen_volumes"""

    def __init__(self, weights, config, vol_args={}, skip_vol=False):
        self.weights = weights
        self.config = config
        self.vol_args = vol_args
        self.skip_vol = skip_vol

    def gen_volumes(self, outdir, z_values):
        if self.skip_vol:
            return
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        zfile = f"{outdir}/z_values.txt"
        np.savetxt(zfile, z_values)
        analysis.gen_volumes(self.weights, self.config, zfile, outdir, **self.vol_args)


def main(args):
    t1 = dt.now()
    E = args.epoch
    workdir = args.workdir
    zfile = f"{workdir}/z.{E}.pkl"
    weights = f"{workdir}/weights.{E}.pkl"
    cfg = (
        f"{workdir}/config.yaml"
        if os.path.exists(f"{workdir}/config.yaml")
        else f"{workdir}/config.pkl"
    )
    outdir = f"{workdir}/analyze.{E}"
    if E == -1:
        zfile = f"{workdir}/z.pkl"
        weights = f"{workdir}/weights.pkl"
        outdir = f"{workdir}/analyze"

    if args.outdir:
        outdir = args.outdir
    logger.info(f"Saving results to {outdir}")
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    z = utils.load_pkl(zfile)
    zdim = z.shape[1]

    vol_args = dict(
        Apix=args.Apix,
        downsample=args.downsample,
        flip=args.flip,
        device=args.device,
        invert=args.invert,
        vol_start_index=args.vol_start_index,
    )
    vg = VolumeGenerator(weights, cfg, vol_args, skip_vol=args.skip_vol)

    if zdim == 1:
        analyze_z1(z, outdir, vg)
    else:
        analyze_zN(
            z,
            outdir,
            vg,
            skip_umap=args.skip_umap,
            num_pcs=args.pc,
            num_ksamples=args.ksample,
        )

    # copy over template if file doesn't exist
    cfg = config.load(cfg)
    if cfg["model_args"]["encode_mode"] == "tilt":
        out_ipynb = f"{outdir}/cryoDRGN_ET_viz.ipynb"
        if not os.path.exists(out_ipynb):
            logger.info("Creating jupyter notebook...")
            ipynb = f"{cryodrgn._ROOT}/templates/cryoDRGN_ET_viz_template.ipynb"
            shutil.copyfile(ipynb, out_ipynb)
        else:
            logger.info(f"{out_ipynb} already exists. Skipping")
        logger.info(out_ipynb)
    else:
        out_ipynb = f"{outdir}/cryoDRGN_viz.ipynb"
        if not os.path.exists(out_ipynb):
            logger.info("Creating jupyter notebook...")
            ipynb = f"{cryodrgn._ROOT}/templates/cryoDRGN_viz_template.ipynb"
            shutil.copyfile(ipynb, out_ipynb)
        else:
            logger.info(f"{out_ipynb} already exists. Skipping")
        logger.info(out_ipynb)

        # copy over template if file doesn't exist
        out_ipynb = f"{outdir}/cryoDRGN_filtering.ipynb"
        if not os.path.exists(out_ipynb):
            logger.info("Creating jupyter notebook...")
            ipynb = f"{cryodrgn._ROOT}/templates/cryoDRGN_filtering_template.ipynb"
            shutil.copyfile(ipynb, out_ipynb)
        else:
            logger.info(f"{out_ipynb} already exists. Skipping")
        logger.info(out_ipynb)

    # copy over template if file doesn't exist
    out_ipynb = f"{outdir}/cryoDRGN_figures.ipynb"
    if not os.path.exists(out_ipynb):
        logger.info("Creating jupyter notebook...")
        ipynb = f"{cryodrgn._ROOT}/templates/cryoDRGN_figures_template.ipynb"
        shutil.copyfile(ipynb, out_ipynb)
    else:
        logger.info(f"{out_ipynb} already exists. Skipping")
    logger.info(out_ipynb)

    logger.info(f"Finished in {dt.now()-t1}")


if __name__ == "__main__":
    matplotlib.use("Agg")  # non-interactive backend
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
