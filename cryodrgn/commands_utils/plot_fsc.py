"""Create a plot of one or more sets of computed Fourier shell correlations.

Example usage
-------------
$ cryodrgn_utils plot_fsc vol1-fsc.txt vol2-fsc.txt
$ cryodrgn_utils plot_fsc vol1-fsc.txt -o vol1-fsc.png
$ cryodrgn_utils plot_fsc fsc-a.txt fsc-b.txt fsc-c.txt -o fsc.png --Apix 2.75

"""
import os
import argparse
from typing import Union, Optional
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import logging

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input", nargs="+", help="input cryoDRGN fsc text files")
    parser.add_argument(
        "-a", "--Apix", type=float, default=1.0, help="physical pixel size in angstrom"
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        default="fsc-plot.png",
        help="output plot file name (.png)",
    )


def plot_fsc_vals(fsc_arr: pd.DataFrame, label: str, **plot_args) -> None:
    """Add this set of FSC curves to the current plot using the given aesthetics."""
    plotting_args = dict(linewidth=3.1, alpha=0.81)
    plotting_args.update(**plot_args)

    # if one of the columns is the pixel resolution, use that as the x-axis...
    if "pixres" in fsc_arr.columns:
        for col in set(fsc_arr.columns) - {"pixres"}:
            plt.plot(fsc_arr.pixres, fsc_arr[col], label=label, **plotting_args)

    # ...otherwise just plot the values sequentially
    else:
        for col in set(fsc_arr.columns):
            plt.plot(fsc_arr[col], label=label, **plotting_args)


def create_fsc_plot(
    fsc_vals: Union[pd.DataFrame, dict[str, pd.DataFrame]],
    outfile: Optional[str] = None,
    Apix: float = 1.0,
    title: Optional[str] = None,
) -> None:
    # Create a subplot
    fig, ax = plt.subplots(figsize=(10, 5))

    if isinstance(fsc_vals, dict):
        for plt_lbl, fsc_array in fsc_vals.items():
            plot_fsc_vals(fsc_array, plt_lbl, linewidth=0.9 + 3.5 / len(fsc_vals))

    elif isinstance(fsc_vals, (np.ndarray, pd.DataFrame)):
        plot_fsc_vals(fsc_vals, "")

    else:
        raise TypeError(f"Unrecognized type for `fsc_vals`: {type(fsc_vals).__name__}!")

    res_given = isinstance(fsc_vals, pd.DataFrame) and fsc_vals.shape[1] == 2
    res_given |= isinstance(fsc_vals, dict) and all(
        fsc_arr.shape[1] == 2 for fsc_arr in fsc_vals.values()
    )

    if res_given:
        use_xticks = np.arange(0.1, 0.6, 0.1)
        xtick_lbls = [f"1/{val:.1f}Å" for val in ((1 / use_xticks) * Apix)]
        plt.xticks(use_xticks, xtick_lbls)

    # titles for axes
    plt.xlabel("Spatial frequency", size=14)
    plt.ylabel("Fourier shell correlation", size=14)

    plt.axhline(y=0.143, color="b", linewidth=1.4)
    plt.grid(True, linewidth=0.7, color="0.5", alpha=0.33)
    plt.xticks(size=10)
    plt.yticks(size=10)
    plt.ylim(0, 1.0)
    plt.xlim(0, 0.5)
    ax.set_aspect(0.3)  # Set the aspect ratio on the plot specifically
    plt.tight_layout()
    plt.subplots_adjust(right=0.8)

    if title:
        plt.title(title)

    # Create the legend on the figure, not the plot
    if isinstance(fsc_vals, dict):
        plt.legend(loc="best", prop={"size": 12})

    if outfile:
        if os.path.dirname(outfile):
            os.makedirs(os.path.dirname(outfile), exist_ok=True)
        plt.savefig(outfile, dpi=300, bbox_inches="tight")
    else:
        plt.show()


def main(args):
    # Load and plot FSC curves from each user-given file (see `add_args()` above)."""
    fsc_arrays = dict()

    for file in args.input:
        fsc_arr = pd.read_csv(file, sep=" ")
        if fsc_arr.shape[1] > 2:
            fsc_arrays = fsc_arr

            if len(args.input) > 1:
                logger.info(
                    f"Only using the first FSC file `{file}` which "
                    f"already contains multiple FSC curves!"
                )
                break
        else:
            fsc_arrays[os.path.splitext(file)[0]] = fsc_arr

    if isinstance(fsc_arrays, dict) and len(fsc_arrays) == 1:
        fsc_arrays = tuple(fsc_arrays.values())[0]

    create_fsc_plot(fsc_arrays, args.outfile, args.Apix)
