"""Interactive filtering of particles plotted using various model variables.

This command opens an interactive interface in which the particles (or tilts) used in
the given reconstruction experiment are plotted as a 2-d scatter plot with trained
model variables as the axes values. This interface allows for lasso-type selection
of particles within this space which can then be saved upon closing the window;
you can also colour the points by a third model variable, or select another pair of
model variables to use as axes.

Note that `cryodrgn analyze` must be run first in the given workdir using the
epoch to filter on!

Example usage
-------------
# If no epoch is given, the default is to find the last available epoch in the workdir
$ cryodrgn filter 00_trainvae

# Choose an epoch yourself; save final selection to `indices.pkl` without prompting
$ cryodrgn filter my_outdir --epoch 30 -f

# Choose another epoch; this time choose file name but pre-select directory to save in
$ cryodrgn filter my_outdir --epoch 30 --sel-dir /data/my_indices/

# If you have done multiple k-means clusterings, you can pick which one to use
$ cryodrgn filter my_outdir/ -k 25

# If you already have indices you can start by plotting them
$ cryodrgn filter my_outdir/01_trainvae --plot-inds candidate-particles.pkl

"""
from __future__ import annotations

import argparse
import logging
import os
import pickle
import re
from typing import Any, Mapping, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from matplotlib import colors
from matplotlib.backend_bases import MouseButton, MouseEvent
from matplotlib.gridspec import GridSpec
from matplotlib.path import Path as PlotPath
from matplotlib.widgets import Button, LassoSelector, RadioButtons
from scipy.spatial import transform

from cryodrgn import analysis, utils
from cryodrgn.dataset import ImageDataset, TiltSeriesData

logger = logging.getLogger(__name__)


def prepare_filter_workspace(
    workdir: str,
    *,
    epoch: int = -1,
    kmeans: int = -1,
    plot_inds: Optional[str] = None,
) -> Tuple[pd.DataFrame, np.ndarray, Optional[np.ndarray]]:
    """Load PCA / UMAP / k-means / poses for filtering; shared by CLI and GIF recorder."""

    train_configs_file = os.path.join(workdir, "config.yaml")
    if not os.path.exists(train_configs_file):
        raise ValueError("Missing config.yaml file in given output folder!")

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
    if not os.path.isdir(anlzdir):
        raise ValueError(
            f"No analysis available for epoch {epoch} "
            f"— first run `cryodrgn analyze {workdir} {epoch}`"
        )
    z = utils.load_pkl(os.path.join(workdir, f"z.{epoch}.pkl"))

    if "poses" in train_configs["dataset_args"]:
        pose_pkl = train_configs["dataset_args"]["poses"]
    else:
        pose_pkl = os.path.join(workdir, f"pose.{epoch}.pkl")

    rot, trans = utils.load_pkl(pose_pkl)
    if train_configs["dataset_args"]["ctf"] is not None:
        ctf_params = utils.load_pkl(train_configs["dataset_args"]["ctf"])
    else:
        ctf_params = None

    pre_indices = None if plot_inds is None else utils.load_pkl(plot_inds)
    if "encode_mode" in train_configs["model_args"]:
        enc_mode = train_configs["model_args"]["encode_mode"]
    else:
        enc_mode = "autodec"

    if isinstance(train_configs["dataset_args"]["ind"], int):
        indices = slice(train_configs["dataset_args"]["ind"])
    elif train_configs["dataset_args"]["ind"] is not None:
        indices = utils.load_pkl(train_configs["dataset_args"]["ind"])
    else:
        indices = None

    imgs_fl = train_configs["dataset_args"]["particles"]
    if ctf_params is not None and enc_mode != "tilt":
        all_indices = np.array(range(ctf_params.shape[0]))
    elif enc_mode == "tilt":
        pt, tp = TiltSeriesData.parse_particle_tilt(imgs_fl)
        all_indices = np.array(range(len(pt)))
    else:
        all_indices = np.array(range(len(ImageDataset(mrcfile=imgs_fl, lazy=True))))

    if indices is not None:
        ctf_params = ctf_params[indices, :] if ctf_params is not None else None
        all_indices = all_indices[indices]

        if "poses" in train_configs["dataset_args"]:
            rot = rot[indices, :, :]
            trans = trans[indices, :]

    pc, _pca = analysis.run_pca(z)
    umap = utils.load_pkl(os.path.join(anlzdir, "umap.pkl"))

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
    znorm = np.sum(z**2, axis=1) ** 0.5

    if rot.shape[0] == z.shape[0]:
        plot_df = analysis.load_dataframe(
            z=z,
            pc=pc,
            euler=transform.Rotation.from_matrix(rot).as_euler("zyz", degrees=True),
            trans=trans,
            labels=kmeans_lbls,
            umap=umap,
            znorm=znorm,
        )
        if ctf_params is not None:
            plot_df["df1"] = ctf_params[:, 2]
            plot_df["df2"] = ctf_params[:, 3]
            plot_df["dfang"] = ctf_params[:, 4]
            plot_df["phase"] = ctf_params[:, 8]

    else:
        plot_df = analysis.load_dataframe(
            z=z, pc=pc, labels=kmeans_lbls, umap=umap, znorm=znorm
        )

    return plot_df, all_indices, pre_indices


def add_args(parser: argparse.ArgumentParser) -> None:
    """Specifies the command-line interface used by `cryodrgn filter`."""
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
        "--force",
        "-f",
        action="store_true",
        help="save selection to `indices.pkl` without prompting",
    )
    parser.add_argument(
        "--plot-inds",
        type=str,
        help="path to a file containing previously selected indices "
        "that will be plotted at the beginning",
    )
    parser.add_argument(
        "--sel-dir",
        type=str,
        help="directory to save the particle selection into",
    )


def main(args: argparse.Namespace) -> None:
    """Launching the interactive interface for filtering particles from command-line."""

    if args.sel_dir is not None:
        if not os.path.exists(args.sel_dir):
            raise ValueError(f"Directory {args.sel_dir} does not exist!")
        elif not os.path.isdir(args.sel_dir):
            raise ValueError(f"Path {args.sel_dir} is not a directory!")

    workdir = args.outdir
    epoch = args.epoch
    kmeans = args.kmeans
    plot_inds = args.plot_inds

    plot_df, all_indices, pre_indices = prepare_filter_workspace(
        workdir,
        epoch=epoch,
        kmeans=kmeans,
        plot_inds=plot_inds,
    )

    # Launch the plot and the interactive command-line prompt; once points are selected,
    # close the figure to avoid interference with other plots
    selector = SelectFromScatter(plot_df, pre_indices)
    selected_indices = [all_indices[i] for i in selector.indices]
    plt.close()

    if selected_indices:
        select_str = " ... ".join(
            [
                ",".join([str(i) for i in selected_indices[:6]]),
                ",".join([str(i) for i in selected_indices[-6:]]),
            ]
        )
        print(
            f"Selected {len(selected_indices)} particles from original list of "
            f"{len(all_indices)} particles numbered "
            f"[{min(all_indices)}, ... , {max(all_indices)}]:\n{select_str}"
        )

        if args.force:
            filename = "indices"
        else:
            if args.sel_dir:
                sel_msg = f"Enter filename to save selection under {args.sel_dir} "
            else:
                sel_msg = "Enter filename to save selection "
            filename = input(sel_msg + "(absolute, without extension):").strip()
        if args.sel_dir:
            filename = os.path.join(args.sel_dir, filename)

        # Saving the selected indices
        if filename:
            selected_dir = os.path.dirname(filename)
            if selected_dir:
                os.makedirs(selected_dir, exist_ok=True)

            selected_full_path = filename + ".pkl"
            print(f"Saving selection to `{selected_full_path}`")
            with open(selected_full_path, "wb") as file:
                pickle.dump(np.array(selected_indices, dtype=int), file)

            # Saving the inverse selection
            inverse_filename = filename + "_inverse.pkl"
            inverse_indices = np.setdiff1d(all_indices, selected_indices)
            print(f"Saving inverse selection to `{inverse_filename}`")
            with open(inverse_filename, "wb") as file:
                pickle.dump(np.array(inverse_indices, dtype=int), file)

    else:
        print("Exiting without having made a selection.")


class SelectFromScatter:
    """An interactive scatterplot for choosing particles using a lasso tool."""

    def __init__(
        self,
        data_table: pd.DataFrame,
        pre_indices: Optional[Sequence[int]] = None,
        *,
        interactive: bool = True,
        figure_kwargs: Optional[Mapping[str, Any]] | None = None,
    ) -> None:
        self.interactive = interactive
        self.data_table = data_table
        self.scatter = None

        # Create a plotting region subdivided into three parts vertically, the middle
        # big part being used for the scatterplot and the thin sides used for legends
        fig_kw: dict[str, Any] = {"constrained_layout": True}
        if figure_kwargs:
            fig_kw.update(dict(figure_kwargs))
        self.fig = plt.figure(**fig_kw)
        gs = self.gridspec()
        self.main_ax = self.fig.add_subplot(gs[:, 1])

        # Find the columns in the given data frame that can be used as plotting
        # covariates based on being a numeric non-index column
        self.select_cols = [
            col for col in data_table.select_dtypes("number").columns if col != "index"
        ]

        self.xcol, self.ycol = "UMAP1", "UMAP2"
        self.color_col = "None"
        self.pnt_colors = None

        # Create user interfaces for selecting the covariates to plot and each axis
        lax = self.fig.add_subplot(gs[0, 0])
        lax.axis("off")
        lax.set_title("choose\nx-axis", size=13)
        self.menu_x = RadioButtons(lax, labels=self.select_cols, active=0)
        self.menu_x.on_clicked(self.update_xaxis)

        rax = self.fig.add_subplot(gs[1, 0])
        rax.axis("off")
        rax.set_title("choose\ny-axis", size=13)
        self.menu_y = RadioButtons(rax, labels=self.select_cols, active=1)
        self.menu_y.on_clicked(self.update_yaxis)

        # Create interface for coloring the plotted points by the values of a covariate
        cax = self.fig.add_subplot(gs[:, 2])
        cax.axis("off")
        cax.set_title("choose\ncolors", size=13)
        self.menu_col = RadioButtons(cax, labels=["None"] + self.select_cols, active=0)
        self.menu_col.on_clicked(self.choose_colors)

        # Add the interface buttons for saving a selection to file and exiting without
        # saving; the save button is only made visible when a selection is made
        self.save_ax = self.fig.add_subplot(gs[2, 0])
        self.exit_ax = self.fig.add_subplot(gs[3, 0])
        self.save_btn = Button(
            self.save_ax,
            "Save Selection",
            color="#164316",
            hovercolor="#01BC01",
        )
        self.exit_btn = Button(
            self.exit_ax,
            "Exit Without Saving",
            color="#601515",
            hovercolor="#BA0B0B",
        )
        self.save_btn.label.set_color("white")
        self.exit_btn.label.set_color("white")
        self.save_btn.on_clicked(self.save_click)
        self.exit_btn.on_clicked(self.exit_click)
        self.save_ax.set_visible(False)
        self.exit_ax.set_visible(True)

        # Create and initialize user interface for selecting points in the scatterplot
        self.lasso = LassoSelector(self.main_ax, onselect=self.choose_points)
        self.indices = pre_indices if pre_indices is not None else list()
        self.pik_txt = None
        self.hover_txt = None

        self.handl_id = self.fig.canvas.mpl_connect(
            "motion_notify_event", self.hover_points
        )
        self.fig.canvas.mpl_connect("button_press_event", self.on_click)
        self.fig.canvas.mpl_connect("button_release_event", self.on_release)

        self.plot()

    def gridspec(self) -> GridSpec:
        """Defines the layout of the plots and menus in the interactive interface."""
        return self.fig.add_gridspec(
            4, 3, width_ratios=[1, 7, 1], height_ratios=[7, 7, 1, 1]
        )

    def plot(self) -> None:
        """Redraw the plot using the current plot info upon e.g. input from user."""
        self.main_ax.clear()
        pnt_colors = ["gray" for _ in range(self.data_table.shape[0])]

        # With a non-empty selection, change the color of the selected points
        # and make the save button visible; otherwise, remove the save button again
        if len(self.indices) > 0:
            for idx in self.indices:
                pnt_colors[idx] = "goldenrod"
            self.save_ax.set_visible(True)
        else:
            self.save_ax.set_visible(False)

        if self.color_col != "None" and len(self.indices) == 0:
            clr_vals = self.data_table[self.color_col]

            if clr_vals.dtype == np.int64:
                use_cmap = sns.color_palette("tab10", as_cmap=True)

                def use_norm(x):
                    return x

            elif clr_vals.min() < 0 < clr_vals.max():
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

        # Plot the points with the chosen colors and axes on a labelled scatterplot
        self.scatter = self.main_ax.scatter(
            x=self.data_table[self.xcol],
            y=self.data_table[self.ycol],
            s=1,
            alpha=0.5,
            c=pnt_colors,
        )
        self.main_ax.set_xlabel(self.xcol, size=17, weight="semibold")
        self.main_ax.set_ylabel(self.ycol, size=17, weight="semibold")
        self.main_ax.set_title("Select Points Manually", size=23, weight="semibold")

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
        if self.interactive:
            plt.show()
            plt.draw()
        else:
            self.fig.canvas.draw()

    def update_xaxis(self, xlbl: str) -> None:
        """When we choose a new x-axis label, we remake the plot with the new axes."""
        self.xcol = xlbl
        self.plot()

    def update_yaxis(self, ylbl: str) -> None:
        """When we choose a new y-axis label, we remake the plot with the new axes."""
        self.ycol = ylbl
        self.plot()

    def choose_colors(self, chosen_colors: str) -> None:
        """User selecting new colors from menu necessitate updating the plot."""
        self.color_col = chosen_colors

        if self.color_col != "None":
            self.indices = list()

        self.plot()

    def choose_points(self, verts: np.array) -> None:
        """Update the chosen points and the plot when points are circled by the user."""
        self.indices = np.where(
            PlotPath(verts).contains_points(self.scatter.get_offsets())
        )[0]

        self.color_col = "None"
        self.menu_col.set_active(0)
        self.plot()

    def hover_points(self, event: MouseEvent) -> None:
        """Update the plot label listing points the mouse is currently hovering over."""

        # Erase any existing annotation for points hovered over
        if self.hover_txt is not None and self.hover_txt.get_text():
            self.hover_txt.set_text("")
            self.fig.canvas.draw_idle()

        # If hovering over the plotting region, find if we are hovering over any points
        if event.inaxes == self.main_ax:
            cont, ix = self.scatter.contains(event)

            if cont:
                ant_lbls = [
                    str(int(self.data_table.iloc[x]["index"])) for x in ix["ind"]
                ]

                # If there are a lot of points we are hovering over, shorten the label
                if len(ant_lbls) > 4:
                    ant_lbls = ant_lbls[:4]
                    ant_txt = ",".join(ant_lbls) + ",..."
                else:
                    ant_txt = ",".join(ant_lbls)

                # Add the new label to the bottom-left corner and redraw the plot
                self.hover_txt = self.main_ax.text(
                    0.01,
                    0.01,
                    f"hovering over particles:\n{ant_txt}",
                    size=11,
                    fontstyle="italic",
                    ha="left",
                    va="bottom",
                    transform=self.main_ax.transAxes,
                )
                self.fig.canvas.draw_idle()

    def on_click(self, event: MouseEvent) -> None:
        """When we click the mouse button to make a selection, we disable hover-text."""
        if hasattr(event, "button") and event.button is MouseButton.LEFT:
            self.fig.canvas.mpl_disconnect(self.handl_id)

    def on_release(self, event: MouseEvent) -> None:
        """When the mouse is released after making a selection, re-enable hover-text."""
        if hasattr(event, "button") and event.button is MouseButton.LEFT:
            self.handl_id = self.fig.canvas.mpl_connect(
                "motion_notify_event", self.hover_points
            )

    def save_click(self, event: MouseEvent) -> None:
        """When the save button is clicked, we close display."""
        if hasattr(event, "button") and event.button is MouseButton.LEFT:
            # Close the plot so we can move on to saving the selection
            plt.close("all")

    def exit_click(self, event: MouseEvent) -> None:
        """When the exit button is clicked, we clear the selection and close display."""
        if hasattr(event, "button") and event.button is MouseButton.LEFT:
            self.indices = list()
            plt.close("all")
