"""Utilities for reading/loading data from .star files."""

import numpy as np
import pandas as pd
from datetime import datetime as dt
from typing import Tuple, Union, Optional, TextIO


def parse_star(starfile: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if not starfile.endswith(".star"):
        raise ValueError(f"{starfile} is not a .star file!")

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


def write_star(
    starfile: str, data: pd.DataFrame, data_optics: Optional[pd.DataFrame] = None
) -> None:
    """Save star data to file, using RELION3.1 or RELION3.0 format if w/o optics."""
    with open(starfile, "w") as f:
        f.write("# Created {}\n".format(dt.now()))
        f.write("\n")

        # RELION3.1
        if data_optics is not None:
            _write_star_block(f, data_optics, block_header="data_optics")
            f.write("\n\n")
            _write_star_block(f, data, block_header="data_particles")

        # RELION3.0
        else:
            _write_star_block(f, data, block_header="data_")


def _write_star_block(
    f: TextIO, data: pd.DataFrame, block_header: str = "data_"
) -> None:
    f.write(f"{block_header}\n\n")
    f.write("loop_\n")
    f.write("\n".join(data.columns))
    f.write("\n")

    # TODO: Assumes header and df ordering is consistent
    for _, vals in data.iterrows():
        f.write(" ".join([str(val) for val in vals]))
        f.write("\n")


class Starfile:
    """A class representing data stored in .star files.

    Attributes
    ----------
    df (pd.DataFrame):  The primary data table of the .star file.
    data_optics (pd.DataFrame):  If RELION3.1, the optics data table in the .star file.

    """

    def __init__(
        self,
        starfile: Optional[str] = None,
        *,
        data: Optional[pd.DataFrame] = None,
        data_optics: Optional[pd.DataFrame] = None,
    ):
        if (starfile is None) == (data is None):
            raise ValueError(
                f"Starfile must be instantiated with "
                f"exactly one of {starfile=} or {data=}!"
            )
        if starfile is not None:
            data, data_optics = parse_star(starfile)

        self.df, self.data_optics = data, data_optics

        if self.relion31:
            if "_rlnOpticsGroup" not in self.data_optics.columns:
                raise ValueError(
                    "Given data optics table does not "
                    "contain a `_rlnOpticsGroup` column!"
                )
            self.data_optics = self.data_optics.set_index("_rlnOpticsGroup", drop=False)

    def write(self, outstar: str) -> None:
        write_star(outstar, data=self.df, data_optics=self.data_optics)

    @property
    def relion31(self) -> bool:
        return self.data_optics is not None

    def __len__(self) -> int:
        return self.df.shape[0]

    def optics_values(
        self, fieldname: str, dtype: Optional[np.dtype] = None
    ) -> Union[None, np.ndarray]:
        """Get sample optics values for a given field using optics table if present."""
        if self.relion31 and fieldname in self.data_optics:
            if "_rlnOpticsGroup" in self.df.columns:
                vals = np.array(
                    [
                        self.data_optics.loc[g, fieldname]
                        for g in self.df["_rlnOpticsGroup"].values
                    ]
                )
            else:
                vals = np.array(
                    [self.data_optics[fieldname][0] for _ in range(self.df.shape[0])]
                )

        elif fieldname in self.df:
            vals = self.df[fieldname].values.reshape(-1)
        else:
            vals = None

        if vals is not None and dtype is not None:
            vals = vals.astype(dtype)

        return vals

    @property
    def apix(self) -> Union[None, np.ndarray]:
        return self.optics_values(fieldname="_rlnImagePixelSize", dtype=np.float32)

    @property
    def resolution(self) -> Union[None, np.ndarray]:
        vals = self.optics_values(fieldname="_rlnImageSize", dtype=np.float32)
        if vals is not None:
            vals = vals.astype(np.int64)

        return vals
