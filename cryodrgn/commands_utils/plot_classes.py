"""Create plots of cryoDRGN model results arranged by given particle class labels.

Class labels are expected to be saved as a pickled .pkl file containing a single
one-dimensional `np.array` object containing categorical data (e.g. string labels
or integers) with length equal to the number of particles in the input dataset.
You can also use a 1D array of indices in {0 ... nparticles - 1) in which case two
classes will be plotted defined by particle membership in this subset index.

For example, to create labels representing three classes chosen at random for
a dataset of 10k particles:
> cryodrgn.utils.save_pkl(np.random.choice(3, 10000), "random_labels.pkl")

Colour palettes can be specified using known seaborn palette names:
seaborn.pydata.org/tutorial/color_palettes.html
Or alternatively a .pkl file storing a palette dictionary:
> cryodrgn.utils.save_pkl(
>     {"Class_A": "blue", "Class_B": "red", "Class_C": "#000FFF"},
>     "my_palette.pkl"
> )

This command will create a new folder [traindir]/analyze.[epoch] in the given training
directory for saving output plots, unless `cryodrgn analyze` has done so already.

Example usages
--------------
$ cryodrgn_utils plot_classes my_work_dir/003_train-vae 49 --labels new_classes.pkl

# Use your own color palette saved to a file
$ cryodrgn_utils plot_classes my_work_dir/003_train-vae 49 --labels new_classes.pkl \
                              --palette my_colours.pkl

# Use a colour palette from the seaborn plotting package
$ cryodrgn_utils plot_classes my_work_dir/003_train-vae 49 --labels new_classes.pkl \
                              --palette rocket

# Save plots to .svg files instead of .pngs, which will preserve resolution with scaling
$ cryodrgn_utils plot_classes 005_train-vae 39 --labels new_classes.pkl --svg

"""
import os
import argparse
from itertools import combinations as combns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from cryodrgn import analysis, config, utils

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    """The command-line arguments added for use with `cryodrgn_utils plot_classes`."""
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
        help="Class labels for use in plotting (.pkl); a pickled numpy array ",
    )

    parser.add_argument(
        "--palette",
        help="Path to class colours for use in plotting, "
        "or the name of a seaborn color palette",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        help="Output directory for output plots (default: [traindir]/analyze.[epoch])",
    )
    parser.add_argument(
        "--skip-umap",
        action="store_true",
        help="Skip running computationally-intensive UMAP step",
    )
    parser.add_argument(
        "--svg",
        action="store_true",
        help="Save figures as rasterized .svg image files as opposed to .png",
    )


def main(args: argparse.Namespace) -> None:
    """Plot reconstruction outputs by given class labels (see `add_args()` above)."""

    if args.epoch == -1:
        z_file = os.path.join(args.traindir, "z.pkl")
    else:
        z_file = os.path.join(args.traindir, f"z.{args.epoch}.pkl")

    if args.outdir is not None:
        outdir = str(args.outdir)
    elif args.epoch == -1:
        outdir = os.path.join(args.traindir, "analyze")
    else:
        outdir = os.path.join(args.traindir, f"analyze.{args.epoch}")

    if os.path.exists(os.path.join(args.traindir, "config.yaml")):
        cfg_file = os.path.join(args.traindir, "config.yaml")
    elif os.path.exists(os.path.join(args.traindir, "config.pkl")):
        cfg_file = os.path.join(args.traindir, "config.pkl")
    else:
        raise ValueError(
            f"Given directory `{args.traindir}` does not appear to be a cryoDRGN "
            "output folder as it does not contain a `config.yaml` or `config.pkl` file!"
        )
    cfgs = config.load(cfg_file)

    if not os.path.exists(z_file):
        if cfgs["cmd"][1] not in {"train_vae", "abinit_het"}:
            logger.warning(
                f"Given cryoDRGN output folder `{args.traindir}` is associated with a "
                f"homogeneous reconstruction experiment (`cryodrgn {cfgs['cmd'][1]}`), "
                "for which no class-based plots are currently available!"
            )
            exit(0)
        else:
            raise ValueError(f"Cannot find saved latent space embeddings `{z_file}`!")

    logger.info(f"Saving results to {outdir}")
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    z_mat = utils.load_pkl(z_file)
    classes = utils.load_pkl(args.labels)
    if classes.ndim != 1:
        raise ValueError("Class labels must be a 1D array!")

    n_imgs = z_mat.shape[0]
    if n_imgs < classes.shape[0]:
        raise ValueError(
            f"Given class labels have more entries ({classes.shape[0]}) than there "
            f"were images in the reconstruction experiment ({n_imgs})!"
        )
    elif n_imgs > classes.shape[0]:
        if not all(
            isinstance(i, (int, np.integer)) and 0 <= i < n_imgs for i in classes
        ):
            raise ValueError(
                f"Given class labels have fewer entries ({classes.shape[0]}) than "
                f"there were images in the reconstruction experiment "
                f"({n_imgs}), so they are being interpreted as indices "
                f"defining a subet of the images, but not all of its entries are "
                f"in the range [0, ..., {n_imgs - 1}]!"
            )
        inds = set(classes)
        classes = np.array([int(i) in inds for i in range(n_imgs)])
        nclasses = 2
    else:
        nclasses = len(set(classes))

    if args.palette:
        if os.path.exists(args.palette):
            palette = utils.load_pkl(args.palette)
        else:
            try:
                palette = sns.color_palette(args.palette, nclasses)
            except ValueError as e:
                raise ValueError(
                    f"Palette {args.palette} is not available in seaborn\n\t({e})!"
                )
    else:
        palette = sns.color_palette("Set1", nclasses)

    def save_figure(out_lbl: str) -> None:
        """Utility for saving the current figure according to the given format."""
        plt.tight_layout()

        if args.svg:
            plt.savefig(os.path.join(outdir, f"{out_lbl}.svg"), format="svg")
        else:
            plt.savefig(os.path.join(outdir, f"{out_lbl}.png"), format="png")

        plt.close()

    z_dim = z_mat.shape[1]
    if z_dim > 2:
        logger.info("Performing principal component analysis...")
        pc, pca = analysis.run_pca(z_mat)

        if not args.skip_umap:
            logger.info("Running UMAP...")
            umap_emb = pd.DataFrame(
                analysis.run_umap(z_mat), columns=["UMAP1", "UMAP2"]
            )
            umap_emb["Class"] = classes

            logger.info("Plotting UMAP clustering densities...")
            g = sns.jointplot(
                data=umap_emb,
                x="UMAP1",
                y="UMAP2",
                kind="kde",
                hue="Class",
                palette=palette,
                height=10,
            )
            g.ax_joint.set_xlabel("UMAP1", size=23, weight="semibold")
            g.ax_joint.set_ylabel("UMAP2", size=23, weight="semibold")
            save_figure("umap_kde_classes")

        logger.info("Plotting PCA clustering densities...")
        for pc1, pc2 in combns(range(3), 2):
            g = sns.jointplot(
                x=pc[:, pc1],
                y=pc[:, pc2],
                kind="kde",
                hue=classes,
                palette=palette,
                height=10,
            )
            g.ax_joint.set_xlabel(
                f"PC{pc1 + 1} ({pca.explained_variance_ratio_[pc1]:.2f})",
                size=23,
                weight="semibold",
            )
            g.ax_joint.set_ylabel(
                f"PC{pc2 + 1} ({pca.explained_variance_ratio_[pc2]:.2f})",
                size=23,
                weight="semibold",
            )
            save_figure(f"z_pca_kde_classes_pc.{pc1 + 1}xpc.{pc2 + 1}")
