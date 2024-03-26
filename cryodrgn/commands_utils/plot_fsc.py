"""Create a plot of one or more sets of computed Fourier shell correlations.

Example usages
--------------
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


def plot_fsc_vals(fsc_arr: pd.DataFrame, label):
    plot_args = dict(linewidth=3.1, alpha=0.81)

    if "pixres" in fsc_arr.columns:
        plt.plot(fsc_arr.pixres, fsc_arr.fsc, label=label, **plot_args)
    elif fsc_arr.shape[1] == 1:
        plt.plot(fsc_arr.iloc[:, 0], label=label, **plot_args)
    else:
        raise ValueError(f"Unrecognized format for\n{fsc_arr=}!")


def create_fsc_plot(
    fsc_vals: Union[pd.DataFrame, dict[str, pd.DataFrame]],
    outfile: Optional[str] = None,
    Apix: float = 1.0,
) -> None:
    # Create a subplot
    fig, ax = plt.subplots(figsize=(10, 5))

    if isinstance(fsc_vals, dict):
        for plt_lbl, fsc_array in fsc_vals.items():
            plot_fsc_vals(fsc_array, plt_lbl)

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
        plt.xlabel("1/resolution", size=14)
        plt.ylabel("Fourier shell correlation", size=14)
    else:
        plt.xlabel("Spatial Frequency", size=14)
        plt.ylabel("Fourier shell correlation", size=14)

    plt.grid(True)
    plt.xticks(size=10)
    plt.yticks(size=10)
    plt.ylim(0, 1.0)
    plt.xlim(0, 0.5)
    ax.set_aspect(0.3)  # Set the aspect ratio on the plot specifically
    plt.tight_layout()
    plt.subplots_adjust(right=0.8)

    # Create the legend on the figure, not the plot
    if isinstance(fsc_vals, dict):
        plt.legend(loc="best", bbox_to_anchor=(0.5, 0.5, 0.5, 0.5), prop={"size": 12})

    if outfile:
        if os.path.dirname(outfile):
            os.makedirs(os.path.dirname(outfile), exist_ok=True)
        plt.savefig(outfile, dpi=300, bbox_inches="tight")
    else:
        plt.show()


def main(args):
    # Load and plot data from each file
    fsc_arrays = dict()
    for file in args.input:
        fsc_arr = pd.read_csv(file, sep=" ", header=None)

        if fsc_arr.shape[1] > 2:
            raise ValueError(
                "FSC text files must either have two columns labelled `pixres` "
                "and `fsc` or one unlabelled column containing FSCs; "
                "see cryodrgn_utils fsc for the former"
            )
        elif fsc_arr.shape[1] == 2:
            fsc_arr = pd.read_csv(file, sep=" ")

        fsc_arrays[os.path.splitext(file)[0]] = fsc_arr

    if len(fsc_arrays) == 1:
        fsc_arrays = tuple(fsc_arrays.values())[0]

    create_fsc_plot(fsc_arrays, args.outfile, args.Apix)
