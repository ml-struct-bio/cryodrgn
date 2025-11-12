"""Visualizing latent space and generating volumes for trained models."""

import os
import argparse
import shutil
import logging
from types import SimpleNamespace
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import torch
import nbformat

import cryodrgn
from cryodrgn import analysis, utils
from cryodrgn import models_ai as models
from cryodrgn.mrcfile import write_mrc
from cryodrgn.lattice import Lattice

TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")
matplotlib.use("Agg")


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "workdir", type=os.path.abspath, help="Directory with cryoDRGN results"
    )
    parser.add_argument(
        "epoch",
        type=int,
        help="Epoch number N to analyze (1-based indexing, corresponding to z.N.pkl, weights.N.pkl)",
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

    group = parser.add_argument_group("Modify number of volumes to generate")
    group.add_argument(
        "--pc",
        type=int,
        default=2,
        help="Number of principal component traversals to generate (default: %(default)s)",
    )
    parser.add_argument(
        "--n-per-pc",
        type=int,
        default=10,
        help="Number of samples of the latent reconstruction space to take along each principal component axis (default: %(default)s)",
    )
    group.add_argument(
        "--ksample",
        type=int,
        default=20,
        help="Number of kmeans samples to generate (default: %(default)s)",
    )

    group = parser.add_argument_group("Volume post-processing")
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
        "--low-pass", type=float, help="Low-pass filter resolution in Angstroms"
    )
    group.add_argument(
        "--crop",
        type=int,
        help="crop volume to this box size after downsampling or low-pass filtering (pixels)",
    )
    group.add_argument(
        "--vol-start-index",
        type=int,
        default=1,
        help="Default value of start index for volume generation (default: %(default)s)",
    )


class VolumeGenerator:
    """Helper class to call analysis.gen_volumes"""

    def __init__(
        self,
        hypervolume,
        lattice,
        zdim,
        invert,
        radius_mask,
        data_norm=(0, 1),
        vol_start_index=1,
    ):
        self.hypervolume = hypervolume
        self.lattice = lattice
        self.zdim = zdim
        self.invert = invert
        self.radius_mask = radius_mask
        self.data_norm = data_norm
        self.vol_start_index = vol_start_index

    def gen_volumes(self, outdir, z_values, suffix=None):
        """
        z_values: [nz, zdim]
        """
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        zfile = f"{outdir}/z_values.txt"
        np.savetxt(zfile, z_values)

        for i, z in enumerate(z_values):
            if suffix is None:
                out_mrc = "{}/{}{:03d}.mrc".format(
                    outdir, "vol_", i + self.vol_start_index
                )
            else:
                out_mrc = "{}/{}{:03d}.mrc".format(outdir, "vol_", suffix)

            vol = models.eval_volume_method(
                self.hypervolume,
                self.lattice,
                self.zdim,
                self.data_norm,
                zval=z,
                radius=self.radius_mask,
            )
            if self.invert:
                vol *= -1

            write_mrc(out_mrc, vol.cpu().numpy().astype(np.float32))


class ModelAnalyzer:
    """An engine for analyzing the output of a reconstruction model.

    Attributes
    ----------
    configs (AnalysisConfigurations):   Values of all parameters that can be
                                        set by the user.
    train_configs (TrainingConfigurations): Parameters that were used when
                                            the model was trained.

    epoch (int): Which epoch will be analyzed.

    skip_umap (bool):   UMAP clustering is relatively computationally intense
                        so sometimes we choose not to do it
    n_per_pc (int):     How many samples of the latent reconstruction space
                        will be taken along each principal component axis.
    """

    @classmethod
    def get_last_cached_epoch(cls, traindir: str) -> int:
        chkpnt_files = [fl for fl in os.listdir(traindir) if fl[:8] == "weights."]

        epoch = (
            -2
            if not chkpnt_files
            else max(
                int(fl.split(".")[1])
                for fl in os.listdir(traindir)
                if fl[:8] == "weights."
            )
        )

        return epoch

    def __init__(
        self, traindir: str, config_vals: dict, train_config_vals: dict
    ) -> None:
        self.logger = logging.getLogger(__name__)

        self.configs = SimpleNamespace(**config_vals)
        self.train_configs = SimpleNamespace(**train_config_vals["training"])
        self.traindir = traindir

        # Find how input data was normalized for training
        self.out_cfgs = {k: v for k, v in train_config_vals.items() if k != "training"}
        if "data_norm_mean" not in self.out_cfgs:
            self.out_cfgs["data_norm_mean"] = 0.0
        if "data_norm_std" not in self.out_cfgs:
            self.out_cfgs["data_norm_std"] = 1.0

        # Use last completed epoch if no epoch given specified by the user
        if self.configs.epoch == -1:
            self.epoch = self.get_last_cached_epoch(traindir)
        else:
            self.epoch = self.configs.epoch

        if self.epoch == -2:
            raise ValueError(
                f"Cannot perform any analyses for output directory `{self.traindir}` "
                f"which does not contain any saved drngai training checkpoints!"
            )

        self.use_cuda = torch.cuda.is_available()
        self.device = torch.device("cuda:0" if self.use_cuda else "cpu")
        self.logger.info(f"Use cuda {self.use_cuda}")

        # Load reconstruction model from the saved checkpoint file
        checkpoint_path = os.path.join(self.traindir, f"weights.{self.epoch}.pkl")
        self.logger.info(f"Loading model from {checkpoint_path}")
        checkpoint = torch.load(checkpoint_path, weights_only=False)

        hypervolume_params = checkpoint["hypervolume_params"]
        hypervolume = models.HyperVolume(**hypervolume_params)
        hypervolume.load_state_dict(checkpoint["hypervolume_state_dict"])
        hypervolume.eval()
        hypervolume.to(self.device)

        lattice = Lattice(
            checkpoint["hypervolume_params"]["resolution"],
            extent=0.5,
            device=self.device,
        )

        self.zdim = checkpoint["hypervolume_params"]["z_dim"]
        radius_mask = (
            checkpoint["output_mask_radius"]
            if "output_mask_radius" in checkpoint
            else None
        )
        self.vg = VolumeGenerator(
            hypervolume,
            lattice,
            self.zdim,
            self.configs.invert,
            radius_mask,
            data_norm=(self.out_cfgs["data_norm_mean"], self.out_cfgs["data_norm_std"]),
        )

        # Load the conformations if the using a heterogeneous reconstruction model
        if self.train_configs.zdim > 0:
            self.z = utils.load_pkl(os.path.join(self.traindir, f"z.{self.epoch}.pkl"))
            self.n_samples = self.z.shape[0]
        else:
            self.z = None
            self.n_samples = None

        # Create an output directory for these analyses
        self.outdir = os.path.join(self.traindir, f"analysis_{self.epoch}")
        os.makedirs(self.outdir, exist_ok=True)

    @staticmethod
    def linear_interpolation(z_0, z_1, n, exclude_last=False):
        delta = 0 if not exclude_last else 1.0 / n
        t = np.linspace(0, 1 - delta, n)[..., None]

        return z_0[None] * (1.0 - t) + z_1[None] * t

    def analyze(self):
        if self.zdim == 0:
            self.logger.info("No analyses available for homogeneous reconstruction!")
            return

        if self.zdim == 1:
            self.analyze_z1()
        else:
            self.analyze_zN()

        # Create Jupyter notebooks for data analysis and visualization by
        # copying them over from the template codebase directory
        out_ipynb = os.path.join(self.outdir, "cryoDRGN-analysis.ipynb")

        if not os.path.exists(out_ipynb):
            self.logger.info("Creating analysis+visualization notebook...")
            ipynb = os.path.join(
                cryodrgn._ROOT, "templates", "cryoDRGN_filtering_template.ipynb"
            )
            shutil.copyfile(ipynb, out_ipynb)
        else:
            self.logger.info(f"{out_ipynb} already exists. Skipping")

        # Edit the notebook with the epoch to analyze
        with open(out_ipynb, "r") as f:
            filter_ntbook = nbformat.read(f, as_version=nbformat.NO_CONVERT)

        filter_ntbook["cells"][3]["source"] = filter_ntbook["cells"][3][
            "source"
        ].replace("EPOCH = None", f"EPOCH = {self.epoch}")

        with open(out_ipynb, "w") as f:
            nbformat.write(filter_ntbook, f)

        self.logger.info("Done")

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
        self.vg.gen_volumes(self.outdir, ztraj)

        kmeans_labels, centers = analysis.cluster_kmeans(
            z[..., None], self.ksample, reorder=False
        )
        centers, centers_ind = analysis.get_nearest_point(z[:, None], centers)

        volpath = os.path.join(self.outdir, f"kmeans{self.configs.ksample}")
        self.vg.gen_volumes(volpath, centers)

    def analyze_zN(self) -> None:
        zdim = self.z.shape[1]

        # Principal component analysis
        self.logger.info("Performing principal component analysis...")
        pc, pca = analysis.run_pca(self.z)
        self.logger.info("Generating volumes...")

        for i in range(self.configs.pc):
            start, end = np.percentile(pc[:, i], (5, 95))
            z_pc = analysis.get_pc_traj(
                pca, self.z.shape[1], self.configs.n_per_pc, i + 1, start, end
            )

            volpath = os.path.join(self.outdir, f"pc{i + 1}_{self.configs.n_per_pc}")
            self.vg.gen_volumes(volpath, z_pc)

        # Kmeans clustering
        self.logger.info("K-means clustering...")
        k = min(self.configs.ksample, self.n_samples)
        if self.n_samples < self.configs.ksample:
            self.logger.warning(f"Changing ksample to # of samples: {self.n_samples}")

        kmeans_labels, centers = analysis.cluster_kmeans(self.z, k)
        centers, centers_ind = analysis.get_nearest_point(self.z, centers)
        kmean_path = os.path.join(self.outdir, f"kmeans{k}")
        os.makedirs(kmean_path, exist_ok=True)

        utils.save_pkl(kmeans_labels, os.path.join(kmean_path, "labels.pkl"))
        np.savetxt(os.path.join(kmean_path, "centers.txt"), centers)
        np.savetxt(os.path.join(kmean_path, "centers_ind.txt"), centers_ind, fmt="%d")

        self.logger.info("Generating volumes...")
        self.vg.gen_volumes(kmean_path, centers)

        # UMAP -- slow step
        umap_emb = None
        if zdim > 2 and not self.configs.skip_umap:
            self.logger.info("Running UMAP...")

            if self.n_samples and self.n_samples < 15:
                n_neighbours = self.n_samples - 1
            else:
                n_neighbours = 15

            umap_emb = analysis.run_umap(self.z, n_neighbors=n_neighbours)
            utils.save_pkl(umap_emb, os.path.join(self.outdir, "umap.pkl"))

        # Make some plots
        self.logger.info("Generating plots...")

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
            z_pc = analysis.get_pc_traj(pca, self.z.shape[1], 10, i + 1, start, end)

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
                    self.logger.warning(
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
            t = np.linspace(start, end, 10, endpoint=True)
            plt.figure(figsize=(4, 4))

            if i > 0 and i == self.configs.pc - 1:
                plt.scatter(pc[:, i - 1], pc[:, i], alpha=0.1, s=1, rasterized=True)
                plt.scatter(np.zeros(10), t, c="cornflowerblue", edgecolor="white")
                plt_pc_labels(i - 1, i)

            else:
                plt.scatter(pc[:, i], pc[:, i + 1], alpha=0.1, s=1, rasterized=True)
                plt.scatter(t, np.zeros(10), c="cornflowerblue", edgecolor="white")
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
                    np.zeros(10), t, c="cornflowerblue", edgecolor="white"
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
                    t, np.zeros(10), c="cornflowerblue", edgecolor="white"
                )
                plt_pc_labels_jointplot(g)

            plt.tight_layout()
            plt.savefig(os.path.join(pc_path, "pca_traversal_hex.png"))
            plt.close()


def main(args: argparse.Namespace) -> None:
    cfg = dict(
        epoch=args.epoch,
        device=args.device,
        outdir=args.outdir,
        skip_vol=args.skip_vol,
        skip_umap=args.skip_umap,
        pc=args.pc,
        n_per_pc=args.n_per_pc,
        ksample=args.ksample,
        apix=args.Apix,
        flip=args.flip,
        invert=args.invert,
        downsample=args.downsample,
    )

    cfg_file = os.path.join(args.workdir, "config.yaml")
    if not os.path.exists(cfg_file):
        raise ValueError(
            f"Config file {cfg_file} does not exist! "
            f"Has the model for {args.workdir} been trained yet?"
        )

    analyzer = ModelAnalyzer(args.workdir, cfg, utils.load_yaml(cfg_file))
    analyzer.analyze()
