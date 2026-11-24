"""Utilities for reading/loading data from .star files.

Example usage
-------------
# Read in the tables of a .star file and look at the first row of each
> from cryodrgn.starfile import parse_star, write_star
> stardata, optics_data = parse_star("output_directory/particles.star")
> print(stardata.df.iloc[0])
> if optics_data is not None:
>     print(optics_data.df.iloc[0])

# The `Starfile` class is useful for easier access to optics values
> from cryodrgn.starfile import Starfile
> starfile = Starfile("output_directory/particles.star")

# These will always be 1D arrays of length same as number of main data table entries
> starfile.apix, starfile.resolution
# Can also access any data fields present in either table
> D = starfile.get_optics_values("_rlnImageSize")
> aberr = starfile.get_optics_values("_rlnSphericalAberration")

"""
import numpy as np
import pandas as pd
from datetime import datetime as dt
from typing import Tuple, Union, Optional, TextIO, Iterable
from typing_extensions import Self


def _coerce_star_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    """Convert numeric .star columns from string to numeric dtypes when possible."""
    for col in df.columns:
        converted = pd.to_numeric(df[col], errors="coerce")
        if converted.notna().all():
            df[col] = converted

    return df


def parse_star(starfile: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Read the data table in a .star file, and the data optics table if present."""
    if not starfile.endswith(".star"):
        raise ValueError(f"{starfile} is not a .star file!")

    comments = list()
    blocks = dict()
    cur_block = None
    with open(starfile, "r") as f:
        while line := f.readline():
            line = line.strip()

            # Here we read in every data table found in the file, and figure out which
            # table is the main one later on (optics table is always `data_optics`)
            if line.startswith("data_"):
                if line in blocks:
                    raise ValueError(f"Multiple `{line}` blocks detected!")

                cur_block = str(line)
                blocks[cur_block] = {"headers": list(), "body": list()}

            # File level comments come before the first `data_` field
            elif cur_block is None:
                comments.append(line)

            # Underscores indicate field names for the subsequent data table values
            elif line.startswith("_"):
                blocks[cur_block]["headers"].append(line.split()[0])

            # Read in data table; ignore everything else
            elif len(line):
                if not line.startswith("#") and not line.startswith("loop_"):
                    vals = line.split()
                    if len(vals):
                        blocks[cur_block]["body"].append(vals)

    blocks = {lbl: block for lbl, block in blocks.items() if block["body"]}
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
        blocks[block] = _coerce_star_dtypes(
            pd.DataFrame(data=blocks[block]["body"], columns=blocks[block]["headers"])
        )

    data_block = None
    for block_lbl, block_vals in blocks.items():
        if block_lbl.startswith("data_") and block_lbl != "data_optics":
            if data_block is None or block_vals.shape[0] > data_block.shape[0]:
                data_block = block_vals.copy()

    if data_block is None:
        raise ValueError(f"Starfile `{starfile}` does not contain a data block!")

    return data_block, blocks["data_optics"] if "data_optics" in blocks else None


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

    # Write the header consisting of the column names and their indices
    for i, col in enumerate(data.columns):
        f.write(" ".join([col, f"#{i+1}"]))
        f.write("\n")

    # Write the values in the same order as the DataFrame columns used in the header
    for _, vals in data.iterrows():
        f.write(" ".join([str(val) for val in vals]))
        f.write("\n")


class Starfile:
    """A class representing data stored in a .star file.

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

    def __eq__(self, other: Self) -> bool:
        """Two .star files are equal if and only if they contain the same data."""
        eq_stat = self.df.equals(other.df) and (self.relion31 == other.relion31)
        if eq_stat and self.relion31:
            eq_stat &= self.data_optics.equals(other.data_optics)

        return eq_stat

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
                    [
                        self.data_optics[fieldname].iloc[0]
                        for _ in range(self.df.shape[0])
                    ]
                )

        # If can't find this field in the optics table, look in the primary data table
        elif fieldname in self.df:
            vals = np.asarray(self.df[fieldname]).ravel()
        else:
            vals = None

        if vals is not None and dtype is not None:
            vals = vals.astype(dtype)

        return vals

    def set_optics_values(self, fieldname: str, vals: Union[float, Iterable]) -> None:
        """Set an existing field, preferring its optics-table column when present.

        Scalars and one-element iterables broadcast to the destination column.
        Otherwise, values may be given per particle row or per optics-table row.
        Per-particle values take precedence when those row counts are equal,
        matching the order returned by ``get_optics_values``. Updates to an
        optics-table column must agree within each particle's optics group.
        """
        if not isinstance(vals, Iterable) or (
            isinstance(vals, np.ndarray) and vals.ndim == 0
        ):
            values = np.asarray([vals])
        else:
            values = np.asarray(list(vals))
        if values.ndim != 1:
            raise ValueError("Optics values must be scalar or one-dimensional!")

        in_optics = self.relion31 and fieldname in self.data_optics
        if not in_optics and fieldname not in self.df:
            raise ValueError(f"Cannot find {fieldname=} in this .star file!")

        n_particles = self.df.shape[0]
        possible_sizes = {1, n_particles}
        if self.relion31:
            possible_sizes |= {self.data_optics.shape[0]}
        if len(values) not in possible_sizes:
            raise ValueError(
                f"Given optics values have length `{len(values)}` "
                f"not in {possible_sizes}!"
            )

        target = self.data_optics if in_optics else self.df
        if len(values) == 1:
            target[fieldname] = values[0]
            return

        per_particle = len(values) == n_particles
        if not in_optics and per_particle:
            self.df[fieldname] = values
            return
        if in_optics and not per_particle:
            # Per-group values follow optics-table row order, not group IDs.
            self.data_optics[fieldname] = values
            return

        if "_rlnOpticsGroup" not in self.df:
            if in_optics and self.data_optics.shape[0] == 1:
                if len(values) == 0:
                    return
                if pd.Series(values).nunique(dropna=False) == 1:
                    self.data_optics[fieldname] = values[0]
                    return
            raise ValueError(
                "Cannot map optics values to particles without _rlnOpticsGroup!"
            )
        if not self.data_optics.index.is_unique:
            raise ValueError("Optics table has duplicate _rlnOpticsGroup IDs!")
        groups = self.df["_rlnOpticsGroup"]
        unknown_groups = groups[~groups.isin(self.data_optics.index)].unique()
        if len(unknown_groups):
            raise ValueError(
                f"Particle optics group(s) {unknown_groups.tolist()} "
                "are not present in the optics table!"
            )

        if in_optics:
            # Preserve unused groups when only the observed particles are updated.
            group_values = self.data_optics[fieldname].to_dict()
            updates = pd.Series(values, index=groups.to_numpy())
            for group, group_updates in updates.groupby(level=0, sort=False):
                if group_updates.nunique(dropna=False) != 1:
                    raise ValueError(
                        f"Conflicting {fieldname} values for optics group {group}!"
                    )
                group_values[group] = group_updates.iloc[0]
            self.data_optics[fieldname] = [
                group_values[group] for group in self.data_optics.index
            ]
        else:
            group_values = dict(zip(self.data_optics.index, values))
            self.df[fieldname] = groups.map(group_values).to_numpy()

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
