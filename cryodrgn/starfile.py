"""Utilities for reading/loading data from .star files."""

import numpy as np
import pandas as pd
from datetime import datetime as dt
from typing import Tuple, Union, Optional, TextIO, Iterable
from typing_extensions import Self


def parse_star(starfile: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Read the data table in a .star file, and the data optics table if present."""
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
    """Append a DataFrame to a file using the .star data block format."""
    f.write(f"{block_header}\n\n")
    f.write("loop_\n")

    # write the header
    f.write("\n".join(data.columns))
    f.write("\n")

    # write the values in the same order as the DataFrame columns used in the header
    for _, vals in data.iterrows():
        f.write(" ".join([str(val) for val in vals]))
        f.write("\n")


class Starfile:
    """A class representing data stored in .star files.

    Attributes
    ----------
    df (pd.DataFrame):  The primary data table of the .star file.
    data_optics (pd.DataFrame):  If RELION3.1, the optics data table in the .star file.

    Example usage
    -------------
    # If using a file, it can be passed as the lone argument
    > starfile = Starfile("mydata_folder/.particles.star")

    # If using data tables, must use keywords
    > starfile = Starfile(data=stardf, data_optics=optics_df)

    # Can also override data optics table found in file (but not `data` as well!)
    > starfile = Starfile("mydata_folder/.particles.star", data_optics=optics_df)

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

    @classmethod
    def load(cls, starfile: str) -> Self:
        """Convenience instantiation method for backwards compatibility."""
        return cls(starfile=starfile)

    def write(self, outstar: str) -> None:
        """Save these data tables to file using the .star format."""
        write_star(outstar, data=self.df, data_optics=self.data_optics)

    @property
    def relion31(self) -> bool:
        """Whether this file is RELION3.1 format, with a data optics table present."""
        return self.data_optics is not None

    def __len__(self) -> int:
        """The number of particle images described by this file."""
        return self.df.shape[0]

    def get_optics_values(
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

        # If can't find this field in the optics table, look in the primary data table
        elif fieldname in self.df:
            vals = self.df[fieldname].values.reshape(-1)
        else:
            vals = None

        if vals is not None and dtype is not None:
            vals = vals.astype(dtype)

        return vals

    def set_optics_values(self, fieldname: str, vals: Union[float, Iterable]) -> None:
        """Set optics values for a given field, updating optics table if present."""
        if not isinstance(vals, Iterable):
            vals = [vals]
        else:
            vals = list(vals)

        possible_sizes = {1, len(self)}
        if self.relion31:
            possible_sizes |= {self.data_optics.shape[0]}
        if len(vals) not in possible_sizes:
            raise ValueError(
                f"Given optics values have length `{len(vals)}` "
                f"not in {possible_sizes}!"
            )

        if self.relion31 and fieldname in self.data_optics:
            if "_rlnOpticsGroup" in self.df.columns:
                if len(vals) in {1, self.data_optics.shape[0]}:
                    self.data_optics.loc[:, fieldname] = vals
                else:
                    self.df.loc = vals
                    self.data_optics.drop(fieldname, axis=1, inplace=True)
            else:
                if len(vals) != 1:
                    raise ValueError(
                        f"No optics mapping for this .star file, and thus new optics "
                        f"values have to be of length one, given {len(vals)=}!"
                    )
                self.data_optics.loc[:, fieldname] = vals

        # If can't find this field in the optics table, look in the primary data table
        elif fieldname in self.df:
            if len(vals) == self.data_optics.shape[0]:
                self.df.loc[:, fieldname] = np.array(
                    [
                        vals[self.data_optics["_rlnOpticsGroup"].index(g)]
                        for g in self.df["_rlnOpticsGroup"].values
                    ]
                )
            else:
                self.df.loc[:, fieldname] = vals
        else:
            raise ValueError(f"Cannot find {fieldname=} in this .star file!")

    @property
    def apix(self) -> Union[None, np.ndarray]:
        """The A/px of each image in this file."""
        return self.get_optics_values(fieldname="_rlnImagePixelSize", dtype=np.float32)

    @property
    def resolution(self) -> Union[None, np.ndarray]:
        """The resolution of each image in this file."""
        vals = self.get_optics_values(fieldname="_rlnImageSize", dtype=np.float32)
        if vals is not None:
            vals = vals.astype(np.int64)

        return vals

    def to_relion30(self) -> pd.DataFrame:
        """Converts this data into a single data table for use with RELION3.0."""
        r30_df = self.df.copy()
        for field in set(self.data_optics.columns) - set(self.df.columns):
            if "OpticsGroup" not in field:
                r30_df[field] = self.get_optics_values(fieldname=field)

        return r30_df
