"""Create a plot of Fourier shell correlations using computed values.

Example usages
--------------
$ cryodrgn_utils plot_fsc vol1-fsc.txt vol2-fsc.txt
$ cryodrgn_utils plot_fsc vol1-fsc.txt -o vol1-fsc.png

"""
import os
import argparse
from typing import Union
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Load data from file
def load_data(file):
    data = np.loadtxt(file)
    x = data[:, 0]
    y = data[:, 1]
    return x, y


# Plot data
def plot_data(x, y, label):
    plt.plot(x, y, label=label)


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
        help="output file name (.png)",
    )


def create_fsc_plot(
    fsc_vals: Union[pd.Series, dict[str, pd.Series]], outfile: str, Apix: float = 1.0
) -> None:
    # Create a subplot
    fig, ax = plt.subplots(figsize=(10, 5))

    if isinstance(fsc_vals, dict):
        for plt_lbl, fsc_array in fsc_vals.items():
            plot_data(*fsc_array, plt_lbl)
    elif isinstance(fsc_vals, pd.Series):
        plot_data(*fsc_vals, "")
    else:
        raise TypeError(f"Unrecognized type for `fsc_vals`: {type(fsc_vals).__name__}!")

    if Apix != 0:
        freq = np.arange(1, 6) * 0.1
        res = ["1/{:.1f}".format(val) for val in ((1 / freq) * Apix)]
        print(res)
        res_text = res
        plt.xticks(np.arange(1, 6) * 0.1, res_text)
        plt.xlabel("1/resolution (1/Å)")
        plt.ylabel("Fourier shell correlation")
    else:
        plt.xlabel("Spatial Frequency")
        plt.ylabel("Fourier shell correlation")

    plt.ylim(0, 1.0)
    plt.xlim(0, 0.5)
    ax.set_aspect(0.3)  # Set the aspect ratio on the plot specifically

    # Create the legend on the figure, not the plot
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", prop={"size": 6})
    plt.grid(True)
    plt.tight_layout()
    plt.subplots_adjust(right=0.8)

    if outfile:
        plt.savefig(outfile, dpi=300, bbox_inches="tight")
    else:
        plt.show()


def main(args):
    # Load and plot data from each file
    fsc_arrays = dict()
    for file in args.input:
        fsc_arrays[os.path.basename(file)] = load_data(file)

    create_fsc_plot(fsc_arrays, args.Apix)
