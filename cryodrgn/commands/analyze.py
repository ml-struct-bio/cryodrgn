"""Visualize latent space and generate volumes."""

import argparse
import os
import os.path
import shutil
from datetime import datetime as dt
import logging
import nbformat

import numpy as np
import torch
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from cryodrgn import _ROOT, analysis, utils
from cryodrgn.models.config import load_configs
from cryodrgn.commands.eval_vol import VolumeEvaluator

logger = logging.getLogger(__name__)
TEMPLATE_DIR = os.path.join(_ROOT, "templates")


def add_args(parser):
    parser.add_argument(
        "outdir", type=os.path.abspath, help="Directory with cryoDRGN results"
    )
    parser.add_argument(
        "epoch",
        type=int,
        help="Epoch number N to analyze (0-based indexing, "
        "corresponding to z.N.pkl, weights.N.pkl)",
    )

    parser.add_argument("--device", type=int, help="Optionally specify CUDA device")
    parser.add_argument(
        "-o",
        "--outdir",
        help="Output directory for analysis results "
        "(default: [workdir]/analyze.[epoch])",
    )
    parser.add_argument(
        "--skip-vol", action="store_true", help="Skip generation of volumes"
    )
    parser.add_argument("--skip-umap", action="store_true", help="Skip running UMAP")

    group = parser.add_argument_group("Extra arguments for volume generation")
    group.add_argument(
        "--Apix",
        type=float,
        help="Pixel size to add to .mrc header "
        "(default is to infer from ctf.pkl file else 1)",
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
        help="Number of principal component traversals "
        "to generate (default: %(default)s)",
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
        help="Default value of start index for volume "
        "generation (default: %(default)s)",
    )


class ModelAnalyzer:
    """An engine for analyzing the output of a reconstruction model.

    Attributes
    ----------
    train_configs (TrainingConfigurations): Parameters that were used when
                                            the model was trained.
    """

    def __init__(
        self,
        train_configs: dict,
        epoch: int = -1,
        skip_vol: bool = False,
        apix: float = None,
        downsample: int = None,
        n_per_pc: int = 10,
        ksample: int = 20,
        vol_start_index: int = 0,
    ) -> None:
        self.train_configs = train_configs
        self.traindir = self.train_configs.outdir
        self.z_dim = self.train_configs.z_dim

        if "run.log" in os.listdir(self.traindir):
            self.log = os.path.join(self.traindir, "run.log")
        elif "training.log" in os.listdir(self.traindir):
            self.log = os.path.join(self.traindir, "training.log")

        else:
            raise ValueError(f"No training log file found in `{self.traindir}`!")

        self.use_cuda = torch.cuda.is_available()
        self.device = torch.device("cuda:0" if self.use_cuda else "cpu")
        logger.info(f"Use cuda {self.use_cuda}")

        if apix:
            self.apix = apix

        # find A/px from CTF if not given
        else:
            if self.configs["dataset_args"]["ctf"]:
                ctf_params = utils.load_pkl(self.configs["dataset_args"]["ctf"])
                orig_apixs = set(ctf_params[:, 1])

                # TODO: add support for multiple optics groups
                if len(orig_apixs) > 1:
                    self.apix = 1.0
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

                    cur_size = self.configs["lattice_args"]["D"] - 1
                    self.apix = round(orig_apix * orig_size / cur_size, 6)
                    logger.info(f"using A/px={self.apix} as per CTF parameters")

            else:
                self.apix = 1.0
                logger.info(
                    "cannot find A/px in CTF parameters, " "defaulting to A/px=1.0"
                )

        # use last completed epoch if no epoch given
        if epoch == -1:
            self.epoch = max(
                int(fl.split(".")[1])
                for fl in os.listdir(self.traindir)
                if fl[:8] == "weights."
            )

        else:
            self.epoch = epoch

        # load model data
        self.weights_file = os.path.join(self.traindir, f"weights.{self.epoch}.pkl")
        if self.z_dim > 0:
            self.z = utils.load_pkl(
                os.path.join(self.traindir, f"conf.{self.epoch}.pkl")
            )
        else:
            self.z = None

        # create an output directory for these analyses
        self.outdir = os.path.join(self.traindir, f"analyze.{self.epoch}")
        os.makedirs(self.outdir, exist_ok=True)
        logger.info(f"Saving results to {self.outdir}")

        if skip_vol:
            self.volume_generator = None
        else:
            self.volume_generator = VolumeEvaluator(
                self.weights_file,
                self.train_configs,
                self.device,
                verbose=False,
                apix=self.apix,
            )

    def generate_volumes(
        self, z_values, voldir, prefix="vol_", suffix=None, start_index=0
    ):
        if not self.configs.skip_vol:
            os.makedirs(voldir, exist_ok=True)
            np.savetxt(os.path.join(voldir, "z_values.txt"), z_values)

            self.volume_generator.produce_volumes(
                z_values, voldir, prefix, suffix, start_index
            )

    def analyze(self):
        if self.z_dim == 0:
            logger.warning("No analyses available for homogeneous reconstruction!")

        if self.z_dim == 1:
            self.analyze_z1()
        else:
            self.analyze_zN()

        # create Jupyter notebooks for data analysis and visualization by
        # copying them over from the template directory
        if self.train_configs.quick_config["capture_setup"] == "spa":
            out_ipynb = os.path.join(self.outdir, "cryoDRGN-analysis.ipynb")

            if not os.path.exists(out_ipynb):
                logger.info("Creating analysis+visualization notebook...")
                ipynb = os.path.join(TEMPLATE_DIR, "analysis-template.ipynb")
                shutil.copyfile(ipynb, out_ipynb)

            else:
                logger.info(f"{out_ipynb} already exists. Skipping")

            # edit the notebook with the epoch to analyze
            with open(out_ipynb, "r") as f:
                viz_ntbook = nbformat.read(f, as_version=nbformat.NO_CONVERT)

            viz_ntbook["cells"][3]["source"] = viz_ntbook["cells"][3]["source"].replace(
                "EPOCH = None", "EPOCH = {self.epoch}"
            )

            with open(out_ipynb, "w") as f:
                nbformat.write(viz_ntbook, f)

        if self.configs.sample_z_idx is not None:
            sampledir = os.path.join(self.outdir, "samples")
            os.makedirs(sampledir, exist_ok=True)

            for z_idx in self.sample_z_idx:
                logger.info(f"Sampling index {z_idx}")
                self.generate_volumes(self.z[z_idx][None], sampledir, suffix=z_idx)

        if self.configs.trajectory_1d is not None:
            logger.info(
                "Generating 1d linear trajectory from "
                f"{self.trajectory_1d[0]} to {self.trajectory_1d[1]} "
                f"({self.trajectory_1d[2]} samples)"
            )

            z_0 = self.z[self.trajectory_1d[0]]
            z_1 = self.z[self.trajectory_1d[1]]
            n_zs = self.trajectory_1d[2]
            z_list = self.linear_interpolation(z_0, z_1, n_zs)

            trajdir = os.path.join(
                self.outdir,
                "trajectories",
                f"1d_{self.trajectory_1d[0]}"
                f"_{self.trajectory_1d[self.trajectory_1d[2]]}",
            )

            os.makedirs(trajdir, exist_ok=True)
            self.generate_volumes(z_list, trajdir)

        if self.configs.direct_traversal_txt is not None:
            dir_traversal_vertices_ind = np.loadtxt(self.configs.direct_traversal_txt)
            travdir = os.path.join(self.outdir, "direct_traversal")
            z_values = np.zeros((0, self.z_dim))

            for i, ind in enumerate(dir_traversal_vertices_ind[:-1]):
                z_0 = self.z[int(int)]
                z_1 = self.z[int(dir_traversal_vertices_ind[i + 1])]
                z_values = np.concatenate(
                    [
                        z_values,
                        self.linear_interpolation(z_0, z_1, 10, exclude_last=True),
                    ],
                    0,
                )

            self.generate_volumes(z_values, travdir)

        if self.configs.z_values_txt is not None:
            z_values = np.loadtxt(self.configs.z_values_txt)
            zvaldir = os.path.join(self.outdir, "trajectory")
            self.generate_volumes(z_values, zvaldir)

        logger.info("Done")

    def analyze_z1(self) -> None:
        """Plotting and volume generation for 1D z"""
        assert self.z.shape[1] == 1
        z = self.z.reshape(-1)
        n = len(z)

        plt.figure(1)
        plt.scatter(np.arange(n), z, alpha=0.1, s=2)
        plt.xlabel("particle")
        plt.ylabel("z")
        plt.savefig(os.path.join(self.outdir, "z.png"))
        plt.close()

        plt.figure(2)
        sns.distplot(z)
        plt.xlabel("z")
        plt.savefig(os.path.join(self.outdir, "z_hist.png"))
        plt.close()

        ztraj = np.percentile(z, np.linspace(5, 95, 10))
        self.generate_volumes(ztraj, self.outdir)

        kmeans_labels, centers = analysis.cluster_kmeans(
            z[..., None], self.ksample, reorder=False
        )
        centers, centers_ind = analysis.get_nearest_point(z[:, None], centers)

        volpath = os.path.join(self.outdir, f"kmeans{self.configs.ksample}")
        self.generate_volumes(centers, volpath)

    def analyze_zN(self) -> None:
        zdim = self.z.shape[1]

        # Principal component analysis
        logger.info("Performing principal component analysis...")
        pc, pca = analysis.run_pca(self.z)
        logger.info("Generating volumes...")

        for i in range(self.configs.pc):
            start, end = np.percentile(pc[:, i], (5, 95))
            z_pc = analysis.get_pc_traj(
                pca, self.z.shape[1], self.configs.n_per_pc, i + 1, start, end
            )

            volpath = os.path.join(self.outdir, f"pc{i + 1}_{self.configs.n_per_pc}")
            self.generate_volumes(z_pc, volpath)

        # kmeans clustering
        logger.info("K-means clustering...")
        k = self.configs.ksample
        kmeans_labels, centers = analysis.cluster_kmeans(self.z, k)
        centers, centers_ind = analysis.get_nearest_point(self.z, centers)
        kmean_path = os.path.join(self.outdir, f"kmeans{k}")
        os.makedirs(kmean_path, exist_ok=True)

        utils.save_pkl(kmeans_labels, os.path.join(kmean_path, "labels.pkl"))
        np.savetxt(os.path.join(kmean_path, "centers.txt"), centers)
        np.savetxt(os.path.join(kmean_path, "centers_ind.txt"), centers_ind, fmt="%d")

        logger.info("Generating volumes...")
        self.generate_volumes(centers, kmean_path)

        # UMAP -- slow step
        umap_emb = None
        if zdim > 2 and not self.configs.skip_umap:
            logger.info("Running UMAP...")
            umap_emb = analysis.run_umap(self.z)
            utils.save_pkl(umap_emb, os.path.join(self.outdir, "umap.pkl"))

        # Make some plots
        logger.info("Generating plots...")

        # Plot learning curve
        loss = analysis.parse_loss(self.log)
        plt.figure(figsize=(4, 4))
        plt.plot(loss)
        plt.xlabel("Epoch")
        plt.ylabel("Loss")
        plt.axvline(
            x=self.epoch, linestyle="--", color="black", label=f"Epoch {self.epoch}"
        )
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(self.outdir, f"learning_curve_epoch{self.epoch}.png"))

        def plt_pc_labels(pc1=0, pc2=1):
            plt.xlabel(f"PC{pc1 + 1} " f"({pca.explained_variance_ratio_[pc1]:.2f})")
            plt.ylabel(f"PC{pc2 + 1} " f"({pca.explained_variance_ratio_[pc2]:.2f})")

        def plt_pc_labels_jointplot(g, pc1=0, pc2=1):
            g.ax_joint.set_xlabel(
                f"PC{pc1 + 1} ({pca.explained_variance_ratio_[pc1]:.2f})"
            )
            g.ax_joint.set_ylabel(
                f"PC{pc2 + 1} ({pca.explained_variance_ratio_[pc2]:.2f})"
            )

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
        plt.savefig(os.path.join(self.outdir, "z_pca.png"))
        plt.close()

        # PCA -- Style 2 -- Scatter, with marginals
        g = sns.jointplot(
            x=pc[:, 0], y=pc[:, 1], alpha=0.1, s=1, rasterized=True, height=4
        )
        plt_pc_labels_jointplot(g)
        plt.tight_layout()
        plt.savefig(os.path.join(self.outdir, "z_pca_marginals.png"))
        plt.close()

        # PCA -- Style 3 -- Hexbin
        g = sns.jointplot(x=pc[:, 0], y=pc[:, 1], height=4, kind="hex")
        plt_pc_labels_jointplot(g)
        plt.tight_layout()
        plt.savefig(os.path.join(self.outdir, "z_pca_hexbin.png"))
        plt.close()

        if umap_emb is not None:
            # Style 1 -- Scatter
            plt.figure(figsize=(4, 4))
            plt.scatter(umap_emb[:, 0], umap_emb[:, 1], alpha=0.1, s=1, rasterized=True)
            plt_umap_labels()
            plt.tight_layout()
            plt.savefig(os.path.join(self.outdir, "umap.png"))
            plt.close()

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
            plt.savefig(os.path.join(self.outdir, "umap_marginals.png"))
            plt.close()

            # Style 3 -- Hexbin / heatmap
            g = sns.jointplot(x=umap_emb[:, 0], y=umap_emb[:, 1], kind="hex", height=4)
            plt_umap_labels_jointplot(g)
            plt.tight_layout()
            plt.savefig(os.path.join(self.outdir, "umap_hexbin.png"))
            plt.close()

        # Plot kmeans sample points
        colors = analysis._get_chimerax_colors(k)
        analysis.scatter_annotate(
            pc[:, 0],
            pc[:, 1],
            centers_ind=centers_ind,
            annotate=True,
            colors=colors,
        )
        plt_pc_labels()
        plt.tight_layout()
        plt.savefig(os.path.join(kmean_path, "z_pca.png"))
        plt.close()

        g = analysis.scatter_annotate_hex(
            pc[:, 0],
            pc[:, 1],
            centers_ind=centers_ind,
            annotate=True,
            colors=colors,
        )
        plt_pc_labels_jointplot(g)
        plt.tight_layout()
        plt.savefig(os.path.join(kmean_path, "z_pca_hex.png"))
        plt.close()

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
            plt.savefig(os.path.join(kmean_path, "umap.png"))
            plt.close()

            g = analysis.scatter_annotate_hex(
                umap_emb[:, 0],
                umap_emb[:, 1],
                centers_ind=centers_ind,
                annotate=True,
                colors=colors,
            )
            plt_umap_labels_jointplot(g)
            plt.tight_layout()
            plt.savefig(os.path.join(kmean_path, "umap_hex.png"))
            plt.close()

        # Plot PC trajectories
        for i in range(self.configs.pc):
            start, end = np.percentile(pc[:, i], (5, 95))
            pc_path = os.path.join(self.outdir, f"pc{i + 1}_{self.configs.n_per_pc}")
            z_pc = analysis.get_pc_traj(
                pca, self.z.shape[1], self.configs.n_per_pc, i + 1, start, end
            )

            if umap_emb is not None:
                # UMAP, colored by PCX
                analysis.scatter_color(
                    umap_emb[:, 0],
                    umap_emb[:, 1],
                    pc[:, i],
                    label=f"PC{i + 1}",
                )
                plt_umap_labels()
                plt.tight_layout()
                plt.savefig(os.path.join(pc_path, "umap.png"))
                plt.close()

                # UMAP, with PC traversal
                z_pc_on_data, pc_ind = analysis.get_nearest_point(self.z, z_pc)
                dists = ((z_pc_on_data - z_pc) ** 2).sum(axis=1) ** 0.5

                if np.any(dists > 2):
                    logger.warning(
                        f"Warning: PC{i + 1} point locations "
                        "in UMAP plot may be inaccurate"
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
                plt.savefig(os.path.join(pc_path, "umap_traversal.png"))
                plt.close()

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
                plt.savefig(os.path.join(pc_path, "umap_traversal_connected.png"))
                plt.close()

            # 10 points, from 5th to 95th percentile of PC1 values
            t = np.linspace(start, end, self.configs.n_per_pc, endpoint=True)
            plt.figure(figsize=(4, 4))

            if i > 0 and i == self.configs.pc - 1:
                plt.scatter(pc[:, i - 1], pc[:, i], alpha=0.1, s=1, rasterized=True)
                plt.scatter(
                    np.zeros(self.configs.n_per_pc),
                    t,
                    c="cornflowerblue",
                    edgecolor="white",
                )
                plt_pc_labels(i - 1, i)

            else:
                plt.scatter(pc[:, i], pc[:, i + 1], alpha=0.1, s=1, rasterized=True)
                plt.scatter(
                    t,
                    np.zeros(self.configs.n_per_pc),
                    c="cornflowerblue",
                    edgecolor="white",
                )
                plt_pc_labels(i, i + 1)

            plt.tight_layout()
            plt.savefig(os.path.join(pc_path, "pca_traversal.png"))
            plt.close()

            if i > 0 and i == self.configs.pc - 1:
                g = sns.jointplot(
                    x=pc[:, i - 1],
                    y=pc[:, i],
                    alpha=0.1,
                    s=1,
                    rasterized=True,
                    height=4,
                )
                g.ax_joint.scatter(
                    np.zeros(self.configs.n_per_pc),
                    t,
                    c="cornflowerblue",
                    edgecolor="white",
                )
                plt_pc_labels_jointplot(g, i - 1, i)

            else:
                g = sns.jointplot(
                    x=pc[:, i],
                    y=pc[:, i + 1],
                    alpha=0.1,
                    s=1,
                    rasterized=True,
                    height=4,
                )
                g.ax_joint.scatter(
                    t,
                    np.zeros(self.configs.n_per_pc),
                    c="cornflowerblue",
                    edgecolor="white",
                )
                plt_pc_labels_jointplot(g)

            plt.tight_layout()
            plt.savefig(os.path.join(pc_path, "pca_traversal_hex.png"))
            plt.close()


def main(args):
    matplotlib.use("Agg")  # non-interactive backend
    t0 = dt.now()

    train_configs = load_configs(args.outdir)
    analyzer = ModelAnalyzer(train_configs, args.epoch, args.downsample)
    analyzer.analyze()

    logger.info(f"Finished in {dt.now() - t0}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
