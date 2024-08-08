"""Utilities for reading/loading data from .star files."""

import os
import numpy as np
import pandas as pd
from typing import Tuple, List, Optional, Union
from datetime import datetime as dt


def parse_star(starfile: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    blocks = dict()
    cur_block = None

    with open(starfile, "r") as f:
        while line := f.readline():
            if line.startswith("data_"):
                if not line.startswith("data_optics"):
                    if "data_" in blocks:
                        raise ValueError("Multiple data blocks detected!")
                    cur_block = "data_"
                else:
                    cur_block = "data_optics"

                blocks[cur_block] = {"headers": list(), "body": list()}

            elif line.startswith("_"):
                blocks[cur_block]["headers"].append(line.split()[0])

            elif not line.startswith("#") and not line.startswith("loop_"):
                vals = line.strip().split()
                if len(vals):
                    blocks[cur_block]["body"].append(vals)

    for block in blocks:
        blocks[block]["body"] = np.array(blocks[block]["body"])
        assert blocks[block]["body"].ndim == 2, (
            f"Error in parsing. Uneven # columns detected in parsing"
            f" {set([len(x) for x in blocks[block]['body']])}."
        )

        assert blocks[block]["body"].shape[1] == len(blocks[block]["headers"]), (
            f"Error in parsing. Number of columns {blocks[block]['body'].shape[1]} "
            f"!= number of headers {blocks[block]['headers']}"
        )

        blocks[block] = pd.DataFrame(
            data=blocks[block]["body"], columns=blocks[block]["headers"]
        )

    if "data_" not in blocks:
        raise ValueError(f"Starfile `{starfile}` does not contain a data block!")

    return blocks["data_"], blocks["data_optics"] if "data_optics" in blocks else None


class Stardata:
    """A lightweight class for simple .star file operations."""

    def __init__(
        self,
        sdata: pd.DataFrame,
        data_optics: Optional[pd.DataFrame] = None,
    ) -> None:
        self.df = sdata
        self.data_optics = data_optics

        if self.data_optics is not None:
            if "_rlnOpticsGroup" not in self.data_optics.columns:
                raise ValueError(
                    "Given data optics table does not "
                    "contain a `_rlnOpticsGroup` column!"
                )

            self.data_optics = self.data_optics.set_index("_rlnOpticsGroup", drop=False)

    @staticmethod
    def from_file(filename):
        return Stardata(*parse_star(filename))

    def __len__(self) -> int:
        return len(self.df)

    @staticmethod
    def _write_block(f, data: pd.DataFrame, block_header: str = "data_"):
        f.write(f"{block_header}\n\n")
        f.write("loop_\n")
        f.write("\n".join(data.columns))
        f.write("\n")

        # TODO: Assumes header and df ordering is consistent
        for _, vals in data.iterrows():
            f.write(" ".join([str(val) for val in vals]))
            f.write("\n")

    def write(self, outstar: str):
        with open(outstar, "w") as f:
            f.write("# Created {}\n".format(dt.now()))
            f.write("\n")

            # RELION 3.1
            if self.data_optics is not None:
                self._write_block(f, self.data_optics, block_header="data_optics")
                f.write("\n\n")
                self._write_block(f, self.df, block_header="data_particles")

            # RELION 3.0
            else:
                self._write_block(f, self.df, block_header="data_")

    @property
    def apix(self) -> Union[None, float, np.ndarray]:
        if self.data_optics is not None and "_rlnImagePixelSize" in self.data_optics:
            if "_rlnOpticsGroup" in self.df.columns:
                apix = np.array(
                    [
                        float(self.data_optics.loc[str(g), "_rlnImagePixelSize"])
                        for g in self.df["_rlnOpticsGroup"].values
                    ]
                )
            else:
                apix = float(self.data_optics["_rlnImagePixelSize"][0])

        elif "_rlnImagePixelSize" in self.df:
            apix = self.df["_rlnImagePixelSize"].values.reshape(-1)
        else:
            apix = None

        return apix

    @property
    def resolution(self) -> Union[None, int, np.ndarray]:
        if self.data_optics is not None and "_rlnImageSize" in self.data_optics:
            if "_rlnOpticsGroup" in self.df.columns:
                res = np.array(
                    [
                        int(float(self.data_optics.loc[str(g), "_rlnImageSize"]))
                        for g in self.df["_rlnOpticsGroup"].values
                    ]
                )
            else:
                res = float(self.data_optics["_rlnImageSize"][0])

        elif "_rlnImageSize" in self.df:
            res = self.df["_rlnImageSize"].values.reshape(-1)
        else:
            res = None

        return res


def prefix_paths(mrcs: List, datadir: str):
    mrcs1 = ["{}/{}".format(datadir, os.path.basename(x)) for x in mrcs]
    mrcs2 = ["{}/{}".format(datadir, x) for x in mrcs]
    try:
        for path in set(mrcs1):
            assert os.path.exists(path)
        mrcs = mrcs1
    except AssertionError:
        for path in set(mrcs2):
            assert os.path.exists(path), f"{path} not found"
        mrcs = mrcs2

    return mrcs
