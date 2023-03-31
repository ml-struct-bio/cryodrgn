"""
Lightweight parser for starfiles
"""
import os
from datetime import datetime as dt
import numpy as np
import pandas as pd
from typing import Optional, List
import cryodrgn.types as types
from cryodrgn import mrc
from cryodrgn.mrc import LazyImage


class Starfile:
    def __init__(
        self,
        headers: Optional[List],
        df: pd.DataFrame,
        data_optics=None,
        relion31=False,
    ):
        if headers:
            assert headers == list(df.columns), f"{headers} != {df.columns}"
        self.headers = headers or list(df.columns)
        self.df = df
        self.data_optics = data_optics
        self.relion31 = relion31

    def __len__(self):
        return len(self.df)

    @classmethod
    def load(cls, starfile: str):
        # detect star file type
        f = open(starfile, "r")
        BLOCK = "data_"
        line = ""
        while 1:
            for line in f:
                if line.startswith(BLOCK):
                    break
            break
        if line.startswith("data_optics"):
            return cls._parse_relion31(starfile)
        else:
            return cls._parse_block(starfile, block_header="data_")

    @classmethod
    def _parse_relion31(cls, starfile: str):
        data_optics = cls._parse_block(starfile, block_header="data_optics")
        s = cls._parse_block(starfile, block_header="data_particles")
        s.data_optics = data_optics
        s.relion31 = True
        return s

    @classmethod
    def _parse_block(cls, starfile: str, block_header: str = "data_"):
        headers = []
        line = ""
        f = open(starfile, "r")
        # get to data block
        while True:
            for line in f:
                if line.startswith(block_header):
                    break
            break
        # get to header loop
        while True:
            for line in f:
                if line.startswith("loop_"):
                    break
            break
        # get list of column headers
        while True:
            headers = []
            for line in f:
                if line.startswith("_"):
                    headers.append(line)
                else:
                    break
            break
        # assume all subsequent lines until empty line is the body
        headers = [h.strip().split()[0] for h in headers]
        body = [line]
        for line in f:
            if line.strip() == "":
                break
            body.append(line)
        # put data into an array and instantiate as dataframe
        words = [line.strip().split() for line in body]
        words = np.array(words)
        assert (
            words.ndim == 2
        ), f"Error in parsing. Uneven # columns detected in parsing {set([len(x) for x in words])}."
        assert words.shape[1] == len(
            headers
        ), f"Error in parsing. Number of columns {words.shape[1]} != number of headers {len(headers)}"
        data = {h: words[:, i] for i, h in enumerate(headers)}
        df = pd.DataFrame(data=data)
        return cls(headers, df)

    def _write_block(self, f, headers, df, block_header="data_"):
        f.write(f"{block_header}\n\n")
        f.write("loop_\n")
        f.write("\n".join(headers))
        f.write("\n")
        for i in df.index:
            # TODO: Assumes header and df ordering is consistent
            f.write(" ".join([str(v) for v in df.loc[i]]))
            f.write("\n")

    def write(self, outstar: str):
        f = open(outstar, "w")
        f.write("# Created {}\n".format(dt.now()))
        f.write("\n")

        if self.relion31:
            assert self.data_optics is not None
            self._write_block(
                f,
                self.data_optics.headers,
                self.data_optics.df,
                block_header="data_optics",
            )
            f.write("\n\n")
            self._write_block(f, self.headers, self.df, block_header="data_particles")
        else:
            self._write_block(f, self.headers, self.df, block_header="data_")

    def get_particles(self, datadir: Optional[str] = None, lazy: bool = True):
        """
        Return particles of the starfile

        Input:
            datadir (str): Overwrite base directories of particle .mrcs
                Tries both substituting the base path and prepending to the path
            If lazy=True, returns list of LazyImage instances, else np.array
        """
        particles = self.df["_rlnImageName"]

        # format is index@path_to_mrc
        particles = [x.split("@") for x in particles]
        ind = [int(x[0]) - 1 for x in particles]  # convert to 0-based indexing
        mrcs = [x[1] for x in particles]
        if datadir is not None:
            mrcs = prefix_paths(mrcs, datadir)
        for path in set(mrcs):
            assert os.path.exists(path), f"{path} not found"
        header = mrc.parse_header(mrcs[0])
        D = header.D  # image size along one dimension in pixels
        dtype = header.dtype
        stride = dtype.itemsize * D * D
        dataset = [
            LazyImage(f, (D, D), dtype, 1024 + ii * stride) for ii, f in zip(ind, mrcs)
        ]
        if not lazy:
            dataset = np.array([x.get() for x in dataset])
        return dataset


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


def csparc_get_particles(
    csfile: str, datadir: Optional[str] = None, lazy: bool = True
) -> types.ImageArray:
    metadata = np.load(csfile)
    ind = metadata["blob/idx"]  # 0-based indexing
    mrcs = metadata["blob/path"].astype(str).tolist()
    if mrcs[0].startswith(">"):  # Remove '>' prefix from paths
        mrcs = [x[1:] for x in mrcs]
    if datadir is not None:
        mrcs = prefix_paths(mrcs, datadir)
    for path in set(mrcs):
        assert os.path.exists(path), f"{path} not found"
    D = metadata[0]["blob/shape"][0]
    dtype = np.float32
    stride = np.float32().itemsize * D * D
    dataset = [
        LazyImage(f, (D, D), dtype, 1024 + ii * stride) for ii, f in zip(ind, mrcs)
    ]
    if not lazy:
        dataset = np.array([x.get() for x in dataset])
    return dataset
