"""Describe the latent space produced by a cryoDRGN model by directly comparing volumes.

Example usage
-------------
$ cryodrgn analyze_landscape 003_abinit-het/ 49

# Sample more volumes from k-means centroids generated from the latent space; use a
# larger box size for the sampled volumes instead of downsampling to 128x128
$ cryodrgn analyze_landscape 005_train-vae/ 39 -N 5000 -d 256

"""
import argparse
import os
import shutil
from collections import Counter
from datetime import datetime as dt
import logging
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import torch
import seaborn as sns
from matplotlib.colors import ListedColormap
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from cryodrgn import analysis, utils
from cryodrgn.commands.analyze import VolumeGenerator
from cryodrgn.mrcfile import parse_mrc, write_mrc
from cryodrgn.masking import cosine_dilation_mask
import cryodrgn.models.config

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    """The command-line arguments for use with `cryodrgn analyze_landscape`."""

    parser.add_argument(
        "workdir", type=os.path.abspath, help="Directory with cryoDRGN results"
    )
    parser.add_argument(
        "epoch",
        type=str,
        help="Epoch number N to analyze "
        "(0-based indexing, corresponding to z.N.pkl, weights.N.pkl)",
    )
    parser.add_argument("--device", type=int, help="Optionally specify CUDA device")
    parser.add_argument(
        "-o",
        "--outdir",
        type=os.path.abspath,
        help="Output directory for landscape analysis results (default: "
        "[workdir]/landscape.[epoch])",
    )
    parser.add_argument("--skip-umap", action="store_true", help="Skip running UMAP")
    parser.add_argument(
        "--vol-ind", type=os.path.abspath, help="Index .pkl for filtering volumes"
    )

    group = parser.add_argument_group("Extra arguments for volume generation")
    group.add_argument(
        "-N",
        "--sketch-size",
        type=int,
        default=1000,
        help="Number of volumes to generate for analysis (default: %(default)s)",
    )
    group.add_argument(
        "--Apix",
        type=float,
        default=1,
        help="Pixel size to add to .mrc header (default: %(default)s A/pix)",
    )
    group.add_argument(
        "--flip", action="store_true", help="Flip handedness of output volume"
    )
    group.add_argument(
        "-d",
        "--downsample",
        type=int,
        default=128,
        help="Downsample volumes to this box size (pixels) (default: %(default)s)",
    )
    group.add_argument(
        "--skip-vol", action="store_true", help="Skip generation of volumes"
    )
    group.add_argument(
        "--vol-start-index",
        type=int,
        default=0,
        help="Default value of start index for volume generation "
        "(default: %(default)s)",
    )

    group = parser.add_argument_group("Extra arguments for mask generation")
    group.add_argument(
        "--thresh",
        type=float,
        help="Density value to threshold for masking (default: "
        "half of max density value)",
    )
    group.add_argument(
        "--dilate",
        type=int,
        default=5,
        help="Dilate initial mask by this amount (default: %(default)s pixels)",
    )
    group.add_argument(
        "--cosine-edge",
        type=int,
        default=0,
        help="Apply a cosine edge to the mask (default: %(default)s pixels)",
    )
    group.add_argument(
        "--mask",
        metavar="MRC",
        type=os.path.abspath,
        help="Path to a custom mask. Must be same box size as generated volumes.",
    )

    group = parser.add_argument_group("Extra arguments for clustering")
    group.add_argument(
        "--linkage",
        default="average",
        help="Linkage for agglomerative clustering (e.g. average, "
        "ward) (default: %(default)s)",
    )
    group.add_argument(
        "-M", type=int, default=10, help="Number of clusters (default: %(default)s)"
    )

    group = parser.add_argument_group("Extra arguments for landscape visualization")
    group.add_argument(
        "--pc-dim",
        type=int,
        default=20,
        help="PCA dimensionality reduction (default: %(default)s)",
    )
    group.add_argument(
        "--plot-dim",
        type=int,
        default=5,
        help="Number of dimensions to plot (default: %(default)s)",
    )


def generate_volumes(z, outdir, vg, K):
    # kmeans clustering
    logger.info("Sketching distribution...")
    kmeans_labels, centers = analysis.cluster_kmeans(z, K, on_data=True, reorder=True)
    centers, centers_ind = analysis.get_nearest_point(z, centers)
    if not os.path.exists(f"{outdir}/kmeans{K}"):
        os.mkdir(f"{outdir}/kmeans{K}")
    utils.save_pkl(kmeans_labels, f"{outdir}/kmeans{K}/labels.pkl")
    np.savetxt(f"{outdir}/kmeans{K}/centers.txt", centers)
    np.savetxt(f"{outdir}/kmeans{K}/centers_ind.txt", centers_ind, fmt="%d")
    logger.info("Generating volumes...")
    vg.gen_volumes(f"{outdir}/kmeans{K}", centers)


def make_mask(
    outdir: str,
    K: int,
    dilate: int,
    cosine_edge: int,
    thresh: float,
    in_mrc: Optional[str] = None,
    Apix: float = 1.0,
    vol_start_index: int = 0,
) -> None:
    """Create a masking filter for use across all volumes produced by k-means."""
    kmean_dir = os.path.join(outdir, f"kmeans{K}")

    if in_mrc is None:
        if thresh is None:
            thresh = np.mean(
                [
                    np.percentile(
                        parse_mrc(
                            os.path.join(kmean_dir, f"vol_{vol_start_index+i:03d}.mrc")
                        )[0],
                        99.99,
                    )
                    / 2.0
                    for i in range(K)
                ]
            )

        def mask_fx(vol):
            return cosine_dilation_mask(
                vol, thresh, dilate, edge_dist=cosine_edge, apix=Apix, verbose=False
            )

        # Combine all masks by taking their union (without loading all into memory)
        vol, _ = parse_mrc(os.path.join(kmean_dir, f"vol_{vol_start_index:03d}.mrc"))
        mask = 1.0 - mask_fx(vol)
        for i in range(1, K):
            vol, _ = parse_mrc(
                os.path.join(kmean_dir, f"vol_{vol_start_index+i:03d}.mrc")
            )
            mask *= 1.0 - mask_fx(vol)
        mask = 1.0 - mask
    else:
        mask = np.array(parse_mrc(in_mrc)[0])

    # Save mask, and then view its slices from three acxes as saved plots
    out_mrc = os.path.join(outdir, "mask.mrc")
    logger.info(f"Saving {out_mrc}")
    write_mrc(out_mrc, mask.astype(np.float32), Apix=Apix)
    view_slices(mask, out_png=os.path.join(outdir, "mask_slices.png"))


def view_slices(y: np.array, out_png: str, D: Optional[int] = None) -> None:
    if D is None:
        D = y.shape[0]

    fig, ax = plt.subplots(1, 3, figsize=(10, 8))
    ax[0].imshow(y[D // 2, :, :])
    ax[1].imshow(y[:, D // 2, :])
    ax[2].imshow(y[:, :, D // 2])

    plt.savefig(out_png)


def choose_cmap(M):
    if M <= 10:
        cmap = "tab10"
    elif M <= 20:
        cmap = "tab20"
    else:
        cmap = ListedColormap(sns.color_palette("husl").as_hex())
    return cmap


def get_colors_for_cmap(cmap, M):
    if M <= 20:
        colors = plt.cm.get_cmap(cmap)(np.arange(M) / (np.ceil(M / 10) * 10))
    else:
        colors = plt.cm.get_cmap(cmap)(np.linspace(0, 1, M))
    return colors


def analyze_volumes(
    outdir,
    K,
    dim,
    M,
    linkage,
    vol_ind=None,
    plot_dim=5,
    particle_ind_orig=None,
    Apix=1,
    vol_start_index=0,
):
    kmean_dir = os.path.join(outdir, f"kmeans{K}")
    cmap = choose_cmap(M)

    # Load mean volume; compute it instead if it does not exist
    vol_mean_fl = os.path.join(kmean_dir, "vol_mean.mrc")
    if not os.path.exists(vol_mean_fl):
        volm = torch.stack(
            [
                torch.Tensor(
                    parse_mrc(
                        os.path.join(kmean_dir, f"vol_{vol_start_index+i:03d}.mrc")
                    )[0]
                )
                for i in range(K)
            ]
        ).mean(dim=0)
        write_mrc(vol_mean_fl, np.array(volm).astype(np.float32), Apix=Apix)
    else:
        volm = torch.Tensor(parse_mrc(vol_mean_fl)[0])

    assert isinstance(volm, torch.Tensor)

    # Load the mask; we will only use the non-zero co-ordinates for analyses
    mask = torch.Tensor(parse_mrc(os.path.join(outdir, "mask.mrc"))[0])
    mask_inds = mask > 0.0
    logger.info(f"{mask_inds.sum()} voxels in mask")

    # load volumes
    vols = torch.stack(
        [
            (
                mask[mask_inds]
                * torch.Tensor(
                    parse_mrc(
                        os.path.join(kmean_dir, f"vol_{vol_start_index+i:03d}.mrc")
                    )[0]
                )[mask_inds]
            ).flatten()
            for i in range(K)
        ]
    )
    vols[vols < 0] = 0

    # load umap
    umap = utils.load_pkl(f"{outdir}/umap.pkl")
    ind = np.loadtxt(f"{outdir}/kmeans{K}/centers_ind.txt").astype(int)

    if vol_ind is not None:
        logger.info(f"Filtering to {len(vol_ind)} volumes")
        vols = vols[vol_ind]
        ind = ind[vol_ind]

    # compute PCA
    pca = PCA(dim)
    pca.fit(vols)
    pc = pca.transform(vols)
    utils.save_pkl(pc, os.path.join(outdir, f"vol_pca_{K}.pkl"))
    utils.save_pkl(pca, os.path.join(outdir, "vol_pca_obj.pkl"))
    logger.info("Explained variance ratio:")
    logger.info(pca.explained_variance_ratio_)

    # save rxn coordinates
    for i in range(plot_dim):
        subdir = f"{outdir}/vol_pcs/pc{i+1}"
        if not os.path.exists(subdir):
            os.makedirs(subdir)
        min_, max_ = pc[:, i].min(), pc[:, i].max()
        logger.info((min_, max_))
        for j, val in enumerate(
            np.linspace(min_, max_, 10, endpoint=True), start=vol_start_index
        ):
            v = volm.clone()
            v[mask_inds] += torch.Tensor(pca.components_[i]) * val
            write_mrc(f"{subdir}/{j}.mrc", np.array(v).astype(np.float32), Apix=Apix)

    # which plots to show???
    def plot(i, j):
        plt.figure()
        plt.scatter(pc[:, i], pc[:, j])
        plt.xlabel(f"Volume PC{i+1} (EV: {pca.explained_variance_ratio_[i]:03f})")
        plt.ylabel(f"Volume PC{j+1} (EV: {pca.explained_variance_ratio_[j]:03f})")
        plt.savefig(os.path.join(outdir, f"vol_pca_{K}_{i+1}_{j+1}.png"))

    for i in range(plot_dim - 1):
        plot(i, i + 1)

    # clustering
    subdir = os.path.join(outdir, f"clustering_L2_{linkage}_{M}")
    os.makedirs(subdir, exist_ok=True)
    cluster = AgglomerativeClustering(n_clusters=M, linkage=linkage)
    labels = cluster.fit_predict(vols)
    utils.save_pkl(labels, os.path.join(subdir, "state_labels.pkl"))

    kmeans_labels = utils.load_pkl(os.path.join(outdir, f"kmeans{K}", "labels.pkl"))
    kmeans_counts = Counter(kmeans_labels)
    for i in range(M):
        vol_i = np.where(labels == i)[0]
        logger.info(f"State {i}: {len(vol_i)} volumes")
        if vol_ind is not None:
            vol_i = np.arange(K)[vol_ind][vol_i]

        vol_fl = os.path.join(kmean_dir, f"vol_{vol_start_index+i:03d}.mrc")
        vol_i_all = torch.stack([torch.Tensor(parse_mrc(vol_fl)[0]) for i in vol_i])
        nparticles = np.array([kmeans_counts[i] for i in vol_i])
        vol_i_mean = np.average(vol_i_all, axis=0, weights=nparticles)
        vol_i_std = (
            np.average((vol_i_all - vol_i_mean) ** 2, axis=0, weights=nparticles) ** 0.5
        )
        write_mrc(
            os.path.join(subdir, f"state_{i}_mean.mrc"),
            vol_i_mean.astype(np.float32),
            Apix=Apix,
        )
        write_mrc(
            os.path.join(subdir, f"state_{i}_std.mrc"),
            vol_i_std.astype(np.float32),
            Apix=Apix,
        )

        os.makedirs(os.path.join(subdir, f"state_{i}"), exist_ok=True)
        for v in vol_i:
            os.symlink(
                os.path.join(kmean_dir, f"vol_{vol_start_index+v:03d}.mrc"),
                os.path.join(subdir, f"state_{i}", f"vol_{vol_start_index+v:03d}.mrc"),
            )

        particle_ind = analysis.get_ind_for_cluster(kmeans_labels, vol_i)
        logger.info(f"State {i}: {len(particle_ind)} particles")
        if particle_ind_orig is not None:
            utils.save_pkl(
                particle_ind_orig[particle_ind],
                os.path.join(subdir, f"state_{i}_particle_ind.pkl"),
            )
        else:
            utils.save_pkl(
                particle_ind, os.path.join(subdir, f"state_{i}_particle_ind.pkl")
            )

    # plot clustering results
    def hack_barplot(counts_):
        if M <= 20:  # HACK TO GET COLORS
            with sns.color_palette(cmap):
                g = sns.barplot(x=np.arange(M), y=counts_)
        else:  # default is husl
            g = sns.barplot(x=np.arange(M), y=counts_)
        return g

    plt.figure()
    counts = Counter(labels)
    g = hack_barplot([counts[i] for i in range(M)])  # type: ignore  (bug in Counter type-checking?)
    for i in range(M):
        g.text(i - 0.1, counts[i] + 2, counts[i])  # type: ignore  (bug in Counter type-checking?)
    plt.xlabel("State")
    plt.ylabel("Count")
    plt.savefig(os.path.join(subdir, "state_volume_counts.png"))

    plt.figure()
    particle_counts = [
        np.sum([kmeans_counts[ii] for ii in np.where(labels == i)[0]]) for i in range(M)
    ]
    g = hack_barplot(particle_counts)
    for i in range(M):
        g.text(i - 0.1, particle_counts[i] + 2, particle_counts[i])
    plt.xlabel("State")
    plt.ylabel("Count")
    plt.savefig(os.path.join(subdir, "state_particle_counts.png"))

    def plot_w_labels(i, j):
        plt.figure()
        plt.scatter(pc[:, i], pc[:, j], c=labels, cmap=cmap)
        plt.xlabel(f"Volume PC{i+1} (EV: {pca.explained_variance_ratio_[i]:03f})")
        plt.ylabel(f"Volume PC{j+1} (EV: {pca.explained_variance_ratio_[j]:03f})")
        plt.savefig(os.path.join(subdir, f"vol_pca_{K}_{i+1}_{j+1}.png"))

    for i in range(plot_dim - 1):
        plot_w_labels(i, i + 1)

    def plot_w_labels_annotated(i, j):
        fig, ax = plt.subplots(figsize=(16, 16))
        plt.scatter(pc[:, i], pc[:, j], c=labels, cmap=cmap)
        annots = np.arange(K)
        if vol_ind is not None:
            annots = annots[vol_ind]
        for ii, k in enumerate(annots):
            ax.annotate(str(k), pc[ii, [i, j]] + np.array([0.1, 0.1]))
        plt.xlabel(f"Volume PC{i+1} (EV: {pca.explained_variance_ratio_[i]:03f})")
        plt.ylabel(f"Volume PC{j+1} (EV: {pca.explained_variance_ratio_[j]:03f})")
        plt.savefig(os.path.join(subdir, f"vol_pca_{K}_annotated_{i+1}_{j+1}.png"))

    for i in range(plot_dim - 1):
        plot_w_labels_annotated(i, i + 1)

    # plot clusters on UMAP
    umap_i = umap[ind]
    fig, ax = plt.subplots(figsize=(8, 8))
    plt.scatter(
        umap[:, 0], umap[:, 1], alpha=0.1, s=1, rasterized=True, color="lightgrey"
    )
    colors = get_colors_for_cmap(cmap, M)
    for i in range(M):
        c = umap_i[np.where(labels == i)]
        plt.scatter(c[:, 0], c[:, 1], label=i, color=colors[i])
    plt.legend()
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.savefig(f"{subdir}/umap.png")

    fig, ax = plt.subplots(figsize=(16, 16))
    plt.scatter(
        umap[:, 0], umap[:, 1], alpha=0.1, s=1, rasterized=True, color="lightgrey"
    )
    plt.scatter(umap_i[:, 0], umap_i[:, 1], c=labels, cmap=cmap)
    annots = np.arange(K)
    if vol_ind is not None:
        annots = annots[vol_ind]
    for i, k in enumerate(annots):
        ax.annotate(str(k), umap_i[i] + np.array([0.1, 0.1]))

    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.savefig(os.path.join(subdir, "umap_annotated.png"))


def main(args: argparse.Namespace) -> None:
    t1 = dt.now()
    logger.info(args)

    E = args.epoch
    workdir = args.workdir
    zfile = os.path.join(workdir, f"z.{E}.pkl")
    weights = os.path.join(workdir, f"weights.{E}.pkl")
    cfg_file = os.path.join(workdir, "config.yaml")
    cfg_file = (
        cfg_file if os.path.exists(cfg_file) else os.path.join(workdir, "config.pkl")
    )
    outdir = args.outdir or os.path.join(workdir, f"landscape.{E}")
    logger.info(f"Saving results to {outdir}")
    os.makedirs(outdir, exist_ok=True)
    z = utils.load_pkl(zfile)
    K = args.sketch_size

    vol_args = dict(
        Apix=args.Apix,
        downsample=args.downsample,
        flip=args.flip,
        device=args.device,
        vol_start_index=args.vol_start_index,
    )
    vg = VolumeGenerator(weights, cfg_file, vol_args, skip_vol=args.skip_vol)

    if args.vol_ind is not None:
        args.vol_ind = utils.load_pkl(args.vol_ind)

    if not args.skip_vol:
        generate_volumes(z, outdir, vg, K)
    else:
        logger.info("Skipping volume generation")

    if args.skip_umap:
        assert os.path.exists(os.path.join(outdir, "umap.pkl"))
        logger.info("Skipping UMAP")
    else:
        umap_fl = os.path.join(workdir, f"analyze.{E}", "umap.pkl")
        logger.info(f"Copying UMAP from {umap_fl}")
        if os.path.exists(umap_fl):
            shutil.copyfile(umap_fl, os.path.join(outdir, "umap.pkl"))
        else:
            raise NotImplementedError

    if args.mask:
        logger.info(f"Using custom mask {args.mask}")
    make_mask(
        outdir,
        K,
        args.dilate,
        args.cosine_edge,
        args.thresh,
        args.mask,
        Apix=args.Apix,
        vol_start_index=args.vol_start_index,
    )

    logger.info("Analyzing volumes...")
    # get particle indices if the dataset was originally filtered
    cfgs = cryodrgn.models.config.load(cfg_file)
    particle_ind = (
        utils.load_pkl(cfgs["dataset_args"]["ind"])
        if cfgs["dataset_args"]["ind"] is not None
        else None
    )
    analyze_volumes(
        outdir,
        K,
        args.pc_dim,
        args.M,
        args.linkage,
        vol_ind=args.vol_ind,
        plot_dim=args.plot_dim,
        particle_ind_orig=particle_ind,
        Apix=args.Apix,
        vol_start_index=args.vol_start_index,
    )

    td = dt.now() - t1
    logger.info(f"Finished in {td}")
