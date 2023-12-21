"""Interactive filtering of mapped particles."""

import os
import pickle
import argparse

import pandas as pd
import yaml
import re
import numpy as np
import logging

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import colors
from matplotlib.backend_bases import Event, MouseButton
from matplotlib.widgets import LassoSelector, RadioButtons
from matplotlib.path import Path as PlotPath
from scipy.spatial.transform import Rotation as RR

from cryodrgn import analysis
from cryodrgn import utils

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument(
        "outdir", type=os.path.abspath, help="experiment output directory"
    )

    parser.add_argument(
        "--epoch",
        "-e",
        type=int,
        default=-1,
        help="which train epoch to use for filtering",
    )
    parser.add_argument(
        "--kmeans",
        "-k",
        type=int,
        default=-1,
        help="which set of k-means clusters to use for filtering",
    )

    parser.add_argument(
        "--plot-inds",
        type=str,
        help="path to a file containing previously selected indices "
             "that will be plotted at the beginning",
        dest='plot_inds'
    )


def main(args) -> None:
    workdir = args.outdir
    epoch = args.epoch
    kmeans = args.kmeans
    plot_inds = args.plot_inds

    train_configs_file = os.path.join(workdir, "config.yaml")
    if not os.path.exists(train_configs_file):
        raise ValueError("Missing config.yaml file " "in given output folder!")

    conf_fls = [fl for fl in os.listdir(workdir) if re.fullmatch(r"z\.[0-9]+\.pkl", fl)]

    if not conf_fls:
        raise NotImplementedError(
            "Filtering is not available for the output "
            "of homogeneous reconstruction!"
        )

    with open(train_configs_file, "r") as f:
        train_configs = yaml.safe_load(f)

    if epoch == -1:
        epoch = max(int(x.split(".")[1]) for x in conf_fls)
        logger.info(f"Using epoch {epoch} for filtering...")

    anlzdir = os.path.join(workdir, f"analyze.{epoch}")
    z = utils.load_pkl(os.path.join(workdir, f"z.{epoch}.pkl"))

    # load poses
    if train_configs["dataset_args"]["do_pose_sgd"]:
        pose_pkl = os.path.join(workdir, f"pose.{epoch}.pkl")

        with open(pose_pkl, "rb") as f:
            rot, trans = pickle.load(f)

    else:
        pose_pkl = train_configs["dataset_args"]["poses"]
        rot, trans = utils.load_pkl(pose_pkl)

    ctf_params = utils.load_pkl(train_configs["dataset_args"]["ctf"])
    all_indices = np.array(range(ctf_params.shape[0]))

    if isinstance(train_configs["dataset_args"]["ind"], int):
        ctf_params = ctf_params[:train_configs["dataset_args"]["ind"], :]
        all_indices = all_indices[:train_configs["dataset_args"]["ind"]]

    elif isinstance(train_configs["dataset_args"]["ind"], str):
        inds = utils.load_pkl(train_configs["dataset_args"]["ind"])
        ctf_params = ctf_params[inds, :]
        rot = rot[inds, :, :]
        trans = trans[inds, :]
        all_indices = all_indices[inds]

    pc, pca = analysis.run_pca(z)
    umap = utils.load_pkl(os.path.join(anlzdir, "umap.pkl"))

    # load preselected indices if they have been specified
    if plot_inds:
        with open(plot_inds, "rb") as file:
            pre_indices = pickle.load(file)
    else:
        pre_indices = None

    if kmeans == -1:
        kmeans_dirs = [
            d
            for d in os.listdir(anlzdir)
            if os.path.isdir(os.path.join(anlzdir, d))
            and re.match(r"^kmeans[0-9]+$", d)
        ]

        if len(kmeans_dirs) == 0:
            raise RuntimeError(
                "Did not find any k-means clustering result "
                "outputs for this experiment!"
            )

        kmeans_dir = os.path.join(anlzdir, kmeans_dirs[0])
        if len(kmeans_dirs) > 1:
            print(
                "Found more than one set of k-means results but no "
                "particular k-means set specified, "
                f"defaulting to {kmeans_dir}"
            )

    else:
        kmeans_dir = os.path.join(anlzdir, f"kmeans{kmeans}")

        if not os.path.exists(kmeans_dir):
            raise ValueError(
                "This experiment does not contain results for "
                f"k-means clustering using k={kmeans}!"
            )

    kmeans_lbls = utils.load_pkl(os.path.join(kmeans_dir, "labels.pkl"))

    plot_df = analysis.load_dataframe(
        z=z,
        pc=pc,
        euler=RR.from_matrix(rot).as_euler("zyz", degrees=True),
        trans=trans,
        labels=kmeans_lbls,
        umap=umap,
        df1=ctf_params[:, 2],
        df2=ctf_params[:, 3],
        dfang=ctf_params[:, 4],
        phase=ctf_params[:, 8],
        znorm=np.sum(z**2, axis=1) ** 0.5,
    )

    selector = SelectFromScatter(plot_df, pre_indices)
    input("Press Enter after making your selection...")
    selected_indices = [all_indices[i] for i in selector.indices]
    plt.close()  # Close the figure to avoid interference with other plots

    select_str = " ... ".join(
        [','.join([str(i) for i in selected_indices[:6]]),
         ','.join([str(i) for i in selected_indices[-6:]])]
        )

    print(f"Selected {len(selected_indices)} particles from original list of "
          f"{len(all_indices)} "
          f"particles numbered [{min(all_indices)}, ... , {max(all_indices)}]:\n{select_str}")

    save_option = (
        input("Do you want to save the selection? (yes/no): ").strip().lower())

    if save_option == "yes":
        filename = input(
            "Enter filename to save selection (absolute, without extension): "
        ).strip()

        # saving the selected indices
        if filename:
            selected_full_path = filename + ".pkl"

            with open(selected_full_path, "wb") as file:
                pickle.dump(selected_indices, file)
            print(f"Selection saved to {selected_full_path}")

            # Saving the inverse selection
            inverse_filename = filename + "_inverse.pkl"
            inverse_indices = np.setdiff1d(all_indices, selected_indices)

            with open(inverse_filename, "wb") as file:
                pickle.dump(inverse_indices, file)

            print(f"Inverse selection saved to {inverse_filename}")

    else:
        print("Exiting without saving selection.")


class SelectFromScatter:
    def __init__(self, data_table: pd.DataFrame, pre_indices: list[int]) -> None:
        self.data_table = data_table
        self.scatter = None

        self.fig = plt.figure(constrained_layout=True)
        gs = self.fig.add_gridspec(2, 3, width_ratios=[1, 7, 1])
        self.main_ax = self.fig.add_subplot(gs[:, 1])

        self.select_cols = [
            col for col in data_table.select_dtypes("number").columns if col != "index"
        ]

        self.xcol, self.ycol = "UMAP1", "UMAP2"
        self.color_col = "None"
        self.pnt_colors = None

        lax = self.fig.add_subplot(gs[0, 0])
        lax.axis("off")
        lax.set_title("choose\nx-axis", size=13)
        self.menu_x = RadioButtons(lax, labels=self.select_cols, active=0)
        self.menu_x.on_clicked(self.choose_xaxis)

        rax = self.fig.add_subplot(gs[1, 0])
        rax.axis("off")
        rax.set_title("choose\ny-axis", size=13)
        self.menu_y = RadioButtons(rax, labels=self.select_cols, active=1)
        self.menu_y.on_clicked(self.choose_yaxis)

        cax = self.fig.add_subplot(gs[:, 2])
        cax.axis("off")
        cax.set_title("choose\ncolors", size=13)
        self.menu_col = RadioButtons(cax, labels=["None"] + self.select_cols, active=0)
        self.menu_col.on_clicked(self.choose_colors)

        self.lasso = LassoSelector(self.main_ax, onselect=self.choose_points)
        self.indices = pre_indices if pre_indices else list()
        self.annot = None
        self.pik_txt = None

        self.handl_id = self.fig.canvas.mpl_connect(
            "motion_notify_event", self.hover_points
        )
        self.fig.canvas.mpl_connect("button_press_event", self.on_click)
        self.fig.canvas.mpl_connect("button_release_event", self.on_release)

        self.plot()

    def plot(self) -> None:
        self.main_ax.clear()
        pnt_colors = ["gray" for _ in range(self.data_table.shape[0])]

        if len(self.indices):
            for idx in self.indices:
                pnt_colors[idx] = "goldenrod"

        elif self.color_col != "None":
            clr_vals = self.data_table[self.color_col]

            if clr_vals.dtype == np.int64:
                use_cmap = sns.color_palette("tab10", as_cmap=True)

                def use_norm(x):
                    return x

            elif clr_vals.min() < 0 and clr_vals.max() > 0:
                use_max = max(abs(clr_vals))
                use_norm = colors.Normalize(vmin=-use_max, vmax=use_max)
                use_cmap = sns.color_palette("Spectral", as_cmap=True)

            elif clr_vals.max() < 0:
                use_norm = colors.Normalize(vmin=clr_vals.min(), vmax=0)
                use_cmap = sns.color_palette("ch:s=1.25,rot=-.7", as_cmap=True)
            else:
                use_norm = colors.Normalize(vmin=0, vmax=clr_vals.max())
                use_cmap = sns.color_palette("ch:s=0.25,rot=-.7", as_cmap=True)

            pnt_colors = use_cmap(use_norm(self.data_table[self.color_col]))

        self.scatter = self.main_ax.scatter(
            x=self.data_table[self.xcol],
            y=self.data_table[self.ycol],
            s=1,
            alpha=0.5,
            color=pnt_colors,
        )

        self.main_ax.set_xlabel(self.xcol, size=17, weight="semibold")
        self.main_ax.set_ylabel(self.ycol, size=17, weight="semibold")
        self.main_ax.set_title("Select Points Manually", size=23, weight="semibold")

        self.annot = self.main_ax.annotate(
            "",
            xy=(0, 0),
            xytext=(20, 20),
            textcoords="offset points",
            bbox=dict(boxstyle="round", fc="w"),
            arrowprops=dict(arrowstyle="->"),
            annotation_clip=False,
        )
        self.annot.set_visible(False)

        self.pik_txt = self.main_ax.text(
            0.99,
            0.01,
            f"# of selected particles: {len(self.indices)}",
            size=11,
            fontstyle="italic",
            ha="right",
            va="bottom",
            transform=self.main_ax.transAxes,
        )

        plt.show(block=False)
        plt.draw()

    def choose_xaxis(self, xlbl: str) -> None:
        self.xcol = xlbl
        self.plot()

    def choose_yaxis(self, ylbl: str) -> None:
        self.ycol = ylbl
        self.plot()

    def choose_colors(self, colors: str) -> None:
        self.color_col = colors

        if self.color_col != "None":
            self.indices = list()

        self.plot()

    def choose_points(self, verts: np.array) -> None:
        self.indices = np.where(
            PlotPath(verts).contains_points(self.scatter.get_offsets())
        )[0]

        self.color_col = "None"
        self.menu_col.set_active(0)
        self.plot()

    def hover_points(self, event: Event):
        vis = self.annot.get_visible()

        if event.inaxes == self.main_ax:
            cont, ix = self.scatter.contains(event)

            if cont:
                pos = self.scatter.get_offsets()[ix["ind"][0]]

                self.annot.xy = pos
                ant_lbls = [
                    str(int(self.data_table.iloc[x]["index"])) for x in ix["ind"]
                ]

                if len(ant_lbls) > 5:
                    ant_lbls = ant_lbls[:5]
                    ant_txt = ",".join(ant_lbls) + ",..."
                else:
                    ant_txt = ",".join(ant_lbls)

                self.annot.set_text(ant_txt)
                self.annot.get_bbox_patch().set_facecolor("0.5")
                self.annot.get_bbox_patch().set_alpha(0.4)

                self.annot.set_visible(True)
                self.fig.canvas.draw_idle()

            else:
                if vis:
                    self.annot.set_visible(False)
                    self.fig.canvas.draw_idle()

    def on_click(self, event: Event) -> None:
        if hasattr(event, "button") and event.button is MouseButton.LEFT:
            self.fig.canvas.mpl_disconnect(self.handl_id)

    def on_release(self, event: Event) -> None:
        if hasattr(event, "button") and event.button is MouseButton.LEFT:
            self.handl_id = self.fig.canvas.mpl_connect(
                "motion_notify_event", self.hover_points
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
