"""Create a plot of one or more sets of computed Fourier shell correlations.

Example usage
-------------
# Plot two curves on the same plot and save it to `fsc-plot.png`
$ cryodrgn_utils plot_fsc vol1-fsc.txt vol2-fsc.txt

# Plot one curve and save it to `vol1-fsc.png`
$ cryodrgn_utils plot_fsc vol1-fsc.txt -o vol1-fsc.png

# Plot three curves at `fsc.png` and use an A/px value other than 1.0
$ cryodrgn_utils plot_fsc fsc-a.txt fsc-b.txt fsc-c.txt -o fsc.png --Apix 2.75

"""
import os
import argparse
import logging
from typing import Union, Optional
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def add_args(parser: argparse.ArgumentParser) -> None:
    """The command-line arguments available for use with `cryodrgn_utils plot_fsc`."""

    parser.add_argument(
        "input",
        nargs="+",
        help="input cryoDRGN FSC text files, with either one column containing FSCs or "
        "two space-delimited columns containing pixres and FSCs",
    )
    parser.add_argument(
        "-a",
        "--Apix",
        type=float,
        help="physical pixel size in angstroms for proper frequency x-axis labels",
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
    fsc_vals: Union[np.ndarray, pd.DataFrame, dict[str, pd.DataFrame]],
    outfile: Optional[str] = None,
    apix: Optional[float] = None,
    title: Optional[str] = None,
) -> None:
    """Plot a given set of Fourier shell correlation values on a single canvas.

    Arguments
    ---------
    fsc_vals:   An array or DataFrame of FSC values, in which case each column will be
                treated as an FSC curve, or a dictionary of FSC curves expressed as
                DataFrames with an optional `pixres` columns.
    outfile:    Where to save the plot. If not given, plot will be displayed on screen.
    apix:       Supply an A/px value for creating proper x-axis frequency labels.
    title:      Optionally add this title to the plot.

    """
    fig, ax = plt.subplots(figsize=(10, 5))

    if isinstance(fsc_vals, dict):
        for plt_lbl, fsc_array in fsc_vals.items():
            plot_fsc_vals(fsc_array, plt_lbl, linewidth=0.9 + 3.5 / len(fsc_vals))

    elif isinstance(fsc_vals, (np.ndarray, pd.DataFrame)):
        plot_fsc_vals(fsc_vals, "")

    else:
        raise TypeError(f"Unrecognized type for `fsc_vals`: {type(fsc_vals).__name__}!")

    res_given = isinstance(fsc_vals, pd.DataFrame) and fsc_vals.shape[1] == 2
    if isinstance(fsc_vals, dict):
        res_given |= all(fsc_arr.shape[1] == 2 for fsc_arr in fsc_vals.values())

    if res_given:
        use_xticks = np.arange(0.1, 0.6, 0.1)
        xtick_lbls = [f"1/{val:.1f}Å" for val in ((1 / use_xticks) * apix)]
        plt.xticks(use_xticks, xtick_lbls)
    elif apix is not None:
        logger.warning(
            f"Supplied A/px={apix} but can't produce frequency x-axis labels if "
            f"input arrays don't have `pixres` columns!"
        )

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
    """Load and plot FSC curves from each user-given file (see `add_args()` above)."""
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
