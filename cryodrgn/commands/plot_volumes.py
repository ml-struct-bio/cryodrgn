"""Visualize reconstructed volumes against other model outputs such as latent z-coords.

Example usage
-------------
$ cryodrgn analyze_volumes 003_abinit-het/ 49

"""
import argparse
import os
import os.path
from datetime import datetime as dt
import logging
import subprocess

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sknetwork as skn
from sknetwork.embedding import Spring
from itertools import combinations as combns
from cryodrgn import analysis, utils, config

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
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
        help="Output directory for analysis results (default: [workdir]/analyze_volumes.[epoch])",
    )
    parser.add_argument(
        "--skip-vol", action="store_true", help="Skip generation of volumes"
    )
    parser.add_argument("--skip-umap", action="store_true", help="Skip running UMAP")

    group = parser.add_argument_group("Extra arguments for volume generation")
    group.add_argument(
        "--Apix",
        type=float,
        help="Pixel size to add to .mrc header (default is to infer from ctf.pkl file else 1)",
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
        default=50,
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


def analyze_zN(z, outdir, vg, skip_umap=False, num_ksamples=20):
    zdim = z.shape[1]

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
    umap_fl = f"{outdir}/umap.pkl"

    if zdim > 2 and not skip_umap:
        if os.path.exists(umap_fl):
            umap_emb = utils.load_pkl(umap_fl)
        else:
            logger.info("Running UMAP...")
            umap_emb = analysis.run_umap(z)
            utils.save_pkl(umap_emb, umap_fl)

    # Make some plots
    logger.info("Generating plots...")
    colors = analysis._get_chimerax_colors(K)

    # Style 1 -- Scatter
    fig, ax = plt.subplots(figsize=(7, 7))
    # ax.scatter(umap_emb[:, 0], umap_emb[:, 1], alpha=0.1, s=1, rasterized=True)

    edge_list = list()
    for i, j in combns(range(K), 2):
        dist_val = ((umap_emb[centers_ind[i]] - umap_emb[centers_ind[j]]) ** 2).sum()
        edge_list.append((i, j, dist_val**-2))

    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    xrng, yrng = xlim[1] - xlim[0], ylim[1] - ylim[0]
    img_xscl, img_yscl = 2 * xrng * K**-0.83, 2 * yrng * K**-0.83

    gph = skn.data.from_edge_list(edge_list, matrix_only=False)
    spring = Spring(2, approx_radius=0.6 * (img_xscl * img_yscl) ** 0.5, n_iter=100000)
    spring_emb = spring.fit_transform(
        gph.adjacency, position_init=umap_emb[centers_ind]
    )
    smin0, smax0 = spring_emb[:, 0].min(), spring_emb[:, 0].max()
    smin1, smax1 = spring_emb[:, 1].min(), spring_emb[:, 1].max()
    spring_emb[:, 0] = (spring_emb[:, 0] - smin0) / (smax0 - smin0)
    spring_emb[:, 1] = (spring_emb[:, 1] - smin1) / (smax1 - smin1)

    subprocess.run("export DISPLAY=:0.0", shell=True)
    for i, (emb_x, emb_y) in enumerate(spring_emb):
        vol_fl = os.path.join(outdir, f"kmeans{K}", f"vol_{i:03d}.mrc")
        vol_lbl = os.path.splitext(vol_fl)[0]
        vol_png = vol_lbl + ".png"

        job_msg = subprocess.run(
            f'chimerax --cmd "open {vol_fl}; '
            f'save {vol_png} supersample 3 height 1000 width 1000; exit";',
            shell=True,
            capture_output=True,
        )
        if not os.path.exists(vol_png):
            print(job_msg.stderr.decode("utf-8").strip())
            exit(0)

        vol_img = plt.imread(vol_png)
        img_x = xlim[0] + emb_x * xrng - img_xscl / 2
        img_y = ylim[0] + emb_y * yrng - img_yscl / 2
        ax.imshow(
            vol_img,
            extent=(img_x, img_x + img_xscl, img_y, img_y + img_yscl),
            zorder=100,
            transform=ax.transData,
            aspect="auto",
        )

        lbl_size = 5 + 180 / K
        ax.text(
            img_x, img_y, i, c=colors[i], size=lbl_size, weight="semibold", zorder=200
        )

    ax.set_xlim((xlim[0] - img_xscl, xlim[1] + img_xscl))
    ax.set_ylim((ylim[0] - img_yscl, ylim[1] + img_yscl))
    ax.axis("off")

    fig.savefig(f"{outdir}/kmeans{K}/umap.png", dpi=300, bbox_inches="tight")
    plt.close()


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

        os.makedirs(outdir, exist_ok=True)
        zfile = f"{outdir}/z_values.txt"
        np.savetxt(zfile, z_values)
        analysis.gen_volumes(self.weights, self.config, zfile, outdir, **self.vol_args)


def main(args: argparse.Namespace) -> None:
    matplotlib.use("Agg")  # non-interactive backend
    t1 = dt.now()
    workdir = args.workdir

    zfile = f"{workdir}/z.{args.epoch}.pkl"
    weights = f"{workdir}/weights.{args.epoch}.pkl"
    cfg = (
        f"{workdir}/config.yaml"
        if os.path.exists(f"{workdir}/config.yaml")
        else f"{workdir}/config.pkl"
    )

    configs = config.load(cfg)
    outdir = f"{workdir}/analyze_volumes.{args.epoch}"

    if args.Apix:
        use_apix = args.Apix

    # find A/px from CTF if not given
    else:
        if configs["dataset_args"]["ctf"]:
            ctf_params = utils.load_pkl(configs["dataset_args"]["ctf"])
            orig_apixs = set(ctf_params[:, 1])

            # TODO: add support for multiple optics groups
            if len(orig_apixs) > 1:
                use_apix = 1.0
                logger.info(
                    "cannot find unique A/px in CTF parameters, "
                    "defaulting to A/px=1.0"
                )

            else:
                orig_apix = tuple(orig_apixs)[0]
                orig_sizes = set(ctf_params[:, 0])
                orig_size = tuple(orig_sizes)[0]

                if len(orig_sizes) > 1:
                    logger.info(
                        "cannot find unique original box size in CTF "
                        f"parameters, defaulting to first found: {orig_size}"
                    )

                cur_size = configs["lattice_args"]["D"] - 1
                use_apix = round(orig_apix * orig_size / cur_size, 6)
                logger.info(f"using A/px={use_apix} as per CTF parameters")

        else:
            use_apix = 1.0
            logger.info("cannot find A/px in CTF parameters, " "defaulting to A/px=1.0")

    if args.epoch == -1:
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
        Apix=use_apix,
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
            num_ksamples=args.ksample,
        )

    logger.info(f"Finished in {dt.now() - t1}")
