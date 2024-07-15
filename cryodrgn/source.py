"""This module provides an `ImageSource` class that makes it easy to work with Image data.

An `ImageSource` can be instantiated with a path to a .star/.mrcs/.txt/.cs file, in either lazy or eager mode.
An `images` method is used at runtime to retrieve 3D Tensors for image data at specified indices.
Chunked access is possible using the `chunks()` method.

Typical usage:
  src = source.ImageSource("hand.mrcs", lazy=True)
  im = src.images(np.arange(1000, 2000))
  assert im.shape == (1000, 64, 64)
  ...
  for chunk in src.chunks(chunksize=20):
    assert chunk.shape == (20, 64, 64)
    ...
"""
import os.path
import sys
from datetime import datetime as dt
import struct
from collections import OrderedDict
from collections.abc import Iterable
from concurrent import futures
import numpy as np
import pandas as pd
from typing import List, Tuple, Optional, Union, Any
from typing_extensions import Self
import logging
import torch


logger = logging.getLogger(__name__)


class ImageSource:
    """An ImageSource is a class that returns a copy of the underlying 3D image data, from a .mrcs/.txt/.star file.

    The images(<indices>) method is used to read images at specified indices as torch Tensors.
    <indices> can be a scalar, a slice, a numpy array, or an iterable of indices where we want to query the data.
    Only square images are supported, of side length D pixels.
    The dimensions of the returned Tensor is (<n_images>, D, D).

    The underlying file can be loaded in lazy (default) mode, which is quick, but defers the actual reading of file(s)
    till images(<indices>) is called. In non-lazy mode, the file(s) are read immediately.

    The `images()` call always returns a copy of the data, whether the `ImageSource` is lazy or not.

    Attributes:
        D: An integer indicating the side length (pixels) of the square images in this `ImageSource`.
        n: An integer indicting the total number of images in this `ImageSource`.
        shape: The shape of the underlying data - (n, D, D).

    """

    @staticmethod
    def from_file(
        filepath: str,
        lazy: bool = True,
        indices: Optional[np.ndarray] = None,
        datadir: str = "",
        max_threads: int = 1,
    ):
        ext = os.path.splitext(filepath)[-1][1:]
        if ext == "star":
            return StarfileSource(
                filepath,
                lazy=lazy,
                datadir=datadir,
                indices=indices,
                max_threads=max_threads,
            )
        elif ext in ("mrc", "mrcs"):
            return MRCFileSource(filepath, lazy=lazy, indices=indices)
        elif ext == "txt":
            return TxtFileSource(
                filepath,
                lazy=lazy,
                indices=indices,
                max_threads=max_threads,
            )
        elif ext == "cs":
            return CsSource(
                filepath,
                lazy=lazy,
                datadir=datadir,
                indices=indices,
                max_threads=max_threads,
            )
        else:
            raise RuntimeError(f"Unrecognized file extension {ext}")

    def __init__(
        self,
        D: int,
        n: int,
        filenames: Union[List[str], str, None] = None,
        max_threads: int = 1,
        dtype: str = "float32",
        lazy: bool = True,
        indices: Optional[np.ndarray] = None,
    ):
        self.n = n
        if indices is None:
            self.indices = np.arange(self.n)
        else:
            self.indices = indices
            # If indices is provided, it overrides self.n
            self.n = len(indices)

        self.D = D
        self.shape = self.n, self.D, self.D

        # Some client calls need to access the original filename(s) associated with a source
        # These are traditionally available as the 'fname' attribute of the LazyImage class, hence only used by
        # calling code when lazy=True
        if filenames is None:
            filenames = [""] * self.n
        elif isinstance(filenames, str):
            filenames = [filenames] * self.n
        else:
            assert len(filenames) == self.n, f"{len(filenames)} != {self.n}"
            filenames = list(filenames)
        self.filenames = np.array(filenames)

        self.lazy = lazy
        self.max_threads = max_threads
        self.dtype = dtype

        self.data = None
        if not self.lazy:
            self.data = ArraySource(self._images(self.indices))

    def __len__(self) -> int:
        return self.n

    def __getitem__(self, item) -> torch.Tensor:
        return self.images(item)

    def _convert_to_ndarray(
        self, indices: Optional[Union[np.ndarray, int, slice, Iterable]] = None
    ) -> np.ndarray:
        if isinstance(indices, np.ndarray):
            pass
        elif indices is None:
            indices = np.arange(self.n)
        elif np.isscalar(indices):
            indices = np.array([indices])
        elif isinstance(indices, Iterable):
            indices = np.array(np.fromiter(indices, int))
        elif isinstance(indices, slice):
            start, stop, step = indices.start, indices.stop, indices.step
            start = start or 0
            stop = min(stop, self.n) if stop is not None else self.n
            step = step or 1
            assert (
                start >= 0 and stop >= 0 and step >= 0
            ), "Only positive start/stop/step supported"
            indices = np.arange(start, stop, step)
        else:
            raise TypeError("Unsupported Type for indices")

        assert isinstance(indices, np.ndarray)
        if np.any(indices >= self.n):
            raise ValueError(f"indices should be < {self.n}")

        return indices

    def images(
        self,
        indices: Optional[Union[np.ndarray, int, slice, Iterable]] = None,
        require_contiguous: bool = False,
    ) -> torch.Tensor:
        indices = self._convert_to_ndarray(indices)
        if self.data:  # cached data
            images = self.data._images(indices, require_contiguous=require_contiguous)
        else:
            # Convert incoming caller indices to indices that this ImageSource will use
            if self.indices is not None:
                indices = np.array(self.indices[indices])
            images = self._images(indices, require_contiguous=require_contiguous)

        return torch.from_numpy(images.astype(self.dtype))

    def _images(
        self, indices: np.ndarray, require_contiguous: bool = False
    ) -> np.ndarray:
        """
        Return images at specified indices.
        Args:
            indices: An ndarray of indices
            require_contiguous: Boolean on whether the method should throw an error if image retrieval
            will entail non-contiguous disk access. Callers can employ this if they insist on efficient
            loading and choose to throw an error instead of falling back on inefficient slower loading.
        Returns:
            Images at specified indices.

        """
        raise NotImplementedError("Subclasses must implement this")

    def chunks(self, chunksize: int = 1000):
        """A generator that returns images in chunks of size `chunksize`.

        Returns:
            A 2-tuple of (<indices>, <torch.Tensor>).

        """
        for i in range(0, self.n, chunksize):
            indices = np.arange(i, min(self.n, i + chunksize))
            yield indices, self.images(indices)


class MRCHeader:
    """MRC header class
    # See ref:
    # MRC2014: Extensions to the MRC format header for electron cryo-microscopy and tomography
    # And:
    # https://www.ccpem.ac.uk/mrc_format/mrc2014.php

    """

    ENDIANNESS = "="
    FIELDS = [
        "nx",
        "ny",
        "nz",  # int
        "mode",  # int
        "nxstart",
        "nystart",
        "nzstart",  # int
        "mx",
        "my",
        "mz",  # int
        "xlen",
        "ylen",
        "zlen",  # float
        "alpha",
        "beta",
        "gamma",  # float
        "mapc",
        "mapr",
        "maps",  # int
        "amin",
        "amax",
        "amean",  # float
        "ispg",
        "next",
        "creatid",  # int, int, short, [pad 10]
        "nversion",  # int, [pad 20]
        "nint",
        "nreal",  # short, [pad 20]
        "imodStamp",
        "imodFlags",  # int
        "idtype",
        "lens",
        "nd1",
        "nd2",
        "vd1",
        "vd2",  # short
        "tilt_ox",
        "tilt_oy",
        "tilt_oz",  # float
        "tilt_cx",
        "tilt_cy",
        "tilt_cz",  # float
        "xorg",
        "yorg",
        "zorg",  # float
        "cmap",
        "stamp",
        "rms",  # char[4], float
        "nlabl",
        "labels",
    ]  # int, char[10][80]
    FSTR = "3ii3i3i3f3f3i3f2ih10xi16x2h20x2i6h6f3f4s4sfi800s"

    DTYPE_FOR_MODE = {
        0: np.int8,
        1: np.int16,
        2: np.float32,
        3: "2h",  # complex number from 2 shorts
        4: np.complex64,
        6: np.uint16,
        12: np.float16,
        16: "3B",
    }  # RBG values
    MODE_FOR_DTYPE = {vv: kk for kk, vv in DTYPE_FOR_MODE.items()}

    MACHST_OFFSET = 213
    MACHST_FOR_ENDIANNESS = {"<": b"\x44\x44\x00\x00", ">": b"\x11\x11\x00\x00"}
    ENDIANNESS_FOR_MACHST = {v: k for k, v in MACHST_FOR_ENDIANNESS.items()}

    def __init__(self, header_values, extended_header=b""):
        self.fields = OrderedDict(zip(self.FIELDS, header_values))
        self.extended_header = extended_header
        self.D = self.fields["nx"]
        self.dtype = self.DTYPE_FOR_MODE[self.fields["mode"]]

    def __str__(self):
        return f"Header: {self.fields}\nExtended header: {self.extended_header}"

    @classmethod
    def parse(cls, fname):
        with open(fname, "rb") as f:
            f.seek(cls.MACHST_OFFSET)
            cls.ENDIANNESS = cls.ENDIANNESS_FOR_MACHST.get(f.read(2), "=")

            f.seek(0)
            # prepend endianness specifier to python struct specification
            STRUCT = struct.Struct(cls.ENDIANNESS + cls.FSTR)
            header = cls(STRUCT.unpack(f.read(1024)))

            extbytes = header.fields["next"]
            extended_header = f.read(extbytes)
            header.extended_header = extended_header

        return header

    @classmethod
    def make_default_header(
        cls,
        *,
        nz=None,
        ny=None,
        nx=None,
        data=None,
        is_vol=True,
        Apix=1.0,
        xorg=0.0,
        yorg=0.0,
        zorg=0.0,
    ) -> Self:
        if data is not None:
            nz, ny, nx = data.shape
        assert nz is not None
        assert ny is not None
        assert nx is not None
        ispg = 1 if is_vol else 0

        if is_vol:
            if data is None:
                raise ValueError("If is_vol=True, data array must be specified")

            if isinstance(data, (np.ndarray, torch.Tensor)):
                dmin, dmax, dmean, rms = data.min(), data.max(), data.mean(), data.std()
            elif isinstance(data, ImageSource):
                imgdata = data.images()
                dmin = imgdata.min().item()
                dmax = imgdata.max().item()
                dmean = imgdata.mean().item()
                rms = imgdata.std().item()
            else:
                raise TypeError(f"Unrecognized type of data: `{type(data).__name__}`!")

        else:  # use undefined values for image stacks
            dmin, dmax, dmean, rms = -1, -2, -3, -1

        vals = [
            nx,
            ny,
            nz,
            2,  # mode = 2 for 32-bit float
            0,
            0,
            0,  # nxstart, nystart, nzstart
            nx,
            ny,
            nz,  # mx, my, mz
            Apix * nx,
            Apix * ny,
            Apix * nz,  # cella
            90.0,
            90.0,
            90.0,  # cellb
            1,
            2,
            3,  # mapc, mapr, maps
            dmin,
            dmax,
            dmean,
            ispg,
            0,  # exthd_size
            0,  # creatid
            20140,  # nversion
            0,
            0,  # nint, nreal
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            xorg,
            yorg,
            zorg,
            b"MAP ",
            cls.MACHST_FOR_ENDIANNESS["<" if sys.byteorder == "little" else ">"],
            rms,  # rms
            0,  # nlabl
            b"\x00" * 800,  # labels
        ]

        return cls(vals)

    def write(self, fh):
        STRUCT = struct.Struct(self.FSTR)
        buf = STRUCT.pack(*list(self.fields.values()))
        fh.write(buf)
        fh.write(self.extended_header)

    def get_apix(self):
        return self.fields["xlen"] / self.fields["nx"]

    def update_apix(self, Apix):
        self.fields["xlen"] = self.fields["nx"] * Apix
        self.fields["ylen"] = self.fields["ny"] * Apix
        self.fields["zlen"] = self.fields["nz"] * Apix

    def get_origin(self):
        return self.fields["xorg"], self.fields["yorg"], self.fields["zorg"]

    def update_origin(self, xorg, yorg, zorg):
        self.fields["xorg"] = xorg
        self.fields["yorg"] = yorg
        self.fields["zorg"] = zorg


class MRCFileSource(ImageSource):
    """An ImageSource that reads an .mrc/.mrcs particle stack."""

    def __init__(
        self, filepath: str, lazy: bool = True, indices: Optional[np.ndarray] = None
    ):
        header = MRCHeader.parse(filepath)
        self.header = header
        self.mrcfile_path = filepath
        self.dtype = header.dtype
        self.start = 1024 + header.fields["next"]  # start of image data
        self.nz, self.ny, self.nx = (
            header.fields["nz"],
            header.fields["ny"],
            header.fields["nx"],
        )
        assert self.ny == self.nx, "Only square images supported"
        self.size = self.ny * self.nx
        self.stride = self.dtype().itemsize * self.size

        super().__init__(
            D=self.ny,
            n=self.nz,
            filenames=filepath,
            max_threads=1,
            dtype=self.dtype,
            lazy=lazy,
            indices=indices,
        )

    def _images(
        self,
        indices: np.ndarray,
        data: Optional[np.ndarray] = None,
        tgt_indices: Optional[np.ndarray] = None,
        require_contiguous: bool = False,
    ) -> np.ndarray:
        with open(self.mrcfile_path) as f:
            if data is None:
                data = np.zeros((len(indices), self.D, self.D), dtype=self.dtype)
                assert (
                    tgt_indices is None
                ), "Target indices can only be specified when passing in a pre-allocated array"
                tgt_indices = np.arange(len(indices))
            else:
                if tgt_indices is not None:
                    assert len(tgt_indices) == len(
                        indices
                    ), "indices/tgt_indices length mismatch"
                else:
                    tgt_indices = np.arange(len(indices))

            assert isinstance(tgt_indices, np.ndarray)
            is_contiguous = np.all(indices == indices[0] + np.arange(len(indices)))
            if require_contiguous:
                assert is_contiguous, "MRC indices are not contiguous."

            if is_contiguous:
                f.seek(self.start)
                offset = indices[0] * self.stride
                # 'offset' in the call below is w.r.t the current position of f
                _data = np.fromfile(
                    f, dtype=self.dtype, count=self.size * len(indices), offset=offset
                ).reshape(-1, self.ny, self.nx)
                data[tgt_indices, ...] = _data

            else:
                for index, tgt_index in zip(indices, tgt_indices):
                    f.seek(self.start)
                    offset = index * self.stride
                    # 'offset' in the call below is w.r.t the current position of f
                    _data = np.fromfile(
                        f, dtype=self.dtype, count=self.size, offset=offset
                    ).reshape(self.ny, self.nx)
                    data[tgt_index, ...] = _data

            return data


def parse_mrc(fname: str) -> Tuple[Any, MRCHeader]:  # type: ignore
    # parse the header
    header = MRCHeader.parse(fname)

    # get the number of bytes in extended header
    extbytes = header.fields["next"]
    start = 1024 + extbytes  # start of image data

    dtype = header.dtype
    nz, ny, nx = header.fields["nz"], header.fields["ny"], header.fields["nx"]

    with open(fname, "rb") as fh:
        fh.read(start)  # skip the header + extended header
        array = np.fromfile(fh, dtype=dtype).reshape((nz, ny, nx))

    return array, header


def write_mrc(
    filename: str,
    array: Union[np.ndarray, torch.Tensor, ImageSource],
    header: Optional[MRCHeader] = None,
    Apix: float = 1.0,
    xorg: float = 0.0,
    yorg: float = 0.0,
    zorg: float = 0.0,
    is_vol: Optional[bool] = None,
    transform_fn=None,
    chunksize: int = 1000,
):
    if header is None:
        if is_vol is None:
            is_vol = (
                len(set(array.shape)) == 1
            )  # Guess whether data is vol or image stack
        header = MRCHeader.make_default_header(
            nz=None,
            ny=None,
            nx=None,
            data=array,
            is_vol=is_vol,
            Apix=Apix,
            xorg=xorg,
            yorg=yorg,
            zorg=zorg,
        )
    else:
        # Older versions of MRCHeader had incorrect cmap and stamp fields.
        # Fix these before writing to disk.
        header.fields["cmap"] = b"MAP "
        if header.ENDIANNESS == "=":
            endianness = {"little": "<", "big": ">"}[sys.byteorder]
        else:
            endianness = header.ENDIANNESS
        header.fields["stamp"] = header.MACHST_FOR_ENDIANNESS[endianness]

    if transform_fn is None:
        transform_fn = lambda chunk, indices: chunk  # noqa: E731

    new_dtype = np.dtype(header.dtype).newbyteorder(header.ENDIANNESS)  # type: ignore
    with open(filename, "wb") as f:
        header.write(f)

        if isinstance(array, ImageSource):
            for i, (indices, chunk) in enumerate(array.chunks(chunksize=chunksize)):
                logger.debug(f"Processing chunk {i}")
                chunk = transform_fn(chunk, indices)
                if isinstance(chunk, torch.Tensor):
                    chunk = np.array(chunk.cpu()).astype(new_dtype)
                f.write(chunk.tobytes())
        else:
            indices = np.arange(array.shape[0])
            array = transform_fn(array, indices)
            if isinstance(array, torch.Tensor):
                array = np.array(array.cpu()).astype(new_dtype)

            assert isinstance(array, np.ndarray)
            f.write(array.tobytes())


class ArraySource(ImageSource):
    """A source that is consults an in-memory Numpy array for data.

    An ArraySource is initialized with an ndarray, and indexes into it to return images at specified indices.
    Note that the `indices` argument to `images()` is still a Numpy array, which means that fancy indexing is
    used to get a fresh copy of the requested data. Callers should be mindful of memory usage by passing in a
    reasonable number of indices, or use `chunks()` to iterate through the source.
    """

    def __init__(self, array: np.ndarray):
        if array.ndim == 2:
            array = array[np.newaxis, ...]
        nz, ny, nx = array.shape
        assert ny == nx, "Only square arrays supported"
        self.array = array

        super().__init__(D=ny, n=nz)

    def _images(self, indices: np.ndarray, require_contiguous: bool = False):
        """
        Return ndarray data at specified indices.
        Note that this implementation chooses to ignore `require_contiguous`
        since fancy indexing on a realized ndarray is "fast enough" for all practical purposes.
        """
        return self.array[indices, ...]


class _MRCDataFrameSource(ImageSource):
    def __init__(
        self,
        df: pd.DataFrame,
        datadir: str = "",
        lazy: bool = True,
        indices: Optional[np.ndarray] = None,
        max_threads: int = 1,
    ):
        assert "__mrc_index" in df.columns
        assert "__mrc_filename" in df.columns
        self.df = df

        if datadir:
            self.df["__mrc_filepath"] = self.df["__mrc_filename"].apply(
                lambda filename: os.path.join(datadir, os.path.basename(filename))
            )
        else:
            self.df["__mrc_filepath"] = self.df["__mrc_filename"]

        # Peek into the first mrc file to get image size
        D = MRCFileSource(self.df["__mrc_filepath"][0]).D
        self._sources = {
            filepath: MRCFileSource(filepath)
            for filepath in self.df["__mrc_filepath"].unique()
        }
        super().__init__(
            D=D,
            n=len(self.df),
            max_threads=max_threads,
            lazy=lazy,
            indices=indices,
        )

    def _images(self, indices: np.ndarray, require_contiguous: bool = False):
        def load_single_mrcs(filepath, df):
            src = self._sources[filepath]

            # df.index indicates the positions where the data needs to be inserted -> return for use by caller
            return df.index, src._images(
                df["__mrc_index"].to_numpy(), require_contiguous=require_contiguous
            )

        data = np.zeros((len(indices), self.D, self.D), dtype=self.dtype)

        # Create a DataFrame corresponding to the indices we're interested in
        batch_df = self.df.iloc[indices].reset_index(drop=True)
        groups = batch_df.groupby("__mrc_filepath")
        max_threads = min(self.max_threads, len(groups))
        with futures.ThreadPoolExecutor(max_threads) as executor:
            to_do = []
            for filepath, _df in groups:
                future = executor.submit(load_single_mrcs, filepath, _df)
                to_do.append(future)

            for future in futures.as_completed(to_do):
                data_indices, _data = future.result()
                for idx, d in enumerate(data_indices):
                    data[d] = _data[idx, :, :] if _data.ndim == 3 else _data
        return data


class CsSource(_MRCDataFrameSource):
    def __init__(
        self,
        filepath: str,
        datadir: str = "",
        lazy: bool = True,
        indices: Optional[np.ndarray] = None,
        max_threads: int = 1,
    ):
        metadata = np.load(filepath)
        blob_indices = metadata["blob/idx"]
        blob_paths = metadata["blob/path"].astype(str).tolist()
        n = len(blob_indices)
        assert len(blob_paths) == n

        # Remove leading ">" from paths, if present
        if blob_paths[0].startswith(">"):
            blob_paths = [p[1:] for p in blob_paths]

        df = pd.DataFrame({"__mrc_index": blob_indices, "__mrc_filename": blob_paths})

        if datadir:
            if not os.path.isabs(datadir):
                datadir = os.path.join(os.path.dirname(filepath), datadir)
        else:
            datadir = os.path.dirname(filepath)

        super().__init__(
            df=df, datadir=datadir, lazy=lazy, indices=indices, max_threads=max_threads
        )


class TxtFileSource(_MRCDataFrameSource):
    def __init__(
        self,
        filepath: str,
        lazy: bool = True,
        indices: Optional[np.ndarray] = None,
        max_threads: int = 1,
    ):
        _paths = []
        filepath_dir = os.path.dirname(filepath)
        for line in open(filepath).readlines():
            path = line.strip()
            if not os.path.isabs(path):
                _paths.append(os.path.join(filepath_dir, path))
            else:
                _paths.append(path)

        _source_lengths = [MRCFileSource(path, lazy=True).n for path in _paths]
        mrc_filename, mrc_index = [], []
        for path, length in zip(_paths, _source_lengths):
            mrc_filename.extend([path] * length)
            mrc_index.append(np.arange(length))
        mrc_index = np.concatenate(mrc_index)
        df = pd.DataFrame(
            data={"__mrc_filename": mrc_filename, "__mrc_index": mrc_index}
        )
        super().__init__(df=df, lazy=lazy, indices=indices, max_threads=max_threads)


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


class StarfileSource(_MRCDataFrameSource):
    def __init__(
        self,
        filepath: str,
        datadir: str = "",
        lazy: bool = True,
        indices: Optional[np.ndarray] = None,
        max_threads: int = 1,
    ):
        sdata, self.data_optics = parse_star(filepath)
        sdata[["__mrc_index", "__mrc_filename"]] = sdata["_rlnImageName"].str.split(
            "@", n=1, expand=True
        )
        sdata["__mrc_index"] = pd.to_numeric(sdata["__mrc_index"]) - 1

        if not datadir:
            datadir = os.path.dirname(filepath)

        super().__init__(
            df=sdata,
            datadir=os.path.abspath(datadir),
            lazy=lazy,
            indices=indices,
            max_threads=max_threads,
        )

    def __len__(self):
        return len(self.df)

    @staticmethod
    def _write_block(f, data: pd.DataFrame, block_header: str = "data_"):
        f.write(f"{block_header}\n\n")
        f.write("loop_\n")
        f.write("\n".join(data.columns))
        f.write("\n")

        # TODO: Assumes header and df ordering is consistent
        for i, vals in data.iterrows():
            f.write(" ".join([str(val) for val in vals]))
            f.write("\n")

    def write(self, outstar: str):
        with open(outstar, "w") as f:
            f.write("# Created {}\n".format(dt.now()))
            f.write("\n")

            # RELION 3.1
            if self.data_optics:
                self._write_block(f, self.data_optics, block_header="data_optics")
                f.write("\n\n")
                self._write_block(f, self.df, block_header="data_particles")

            # RELION 3.0
            else:
                self._write_block(f, self.df, block_header="data_")

    def get_apix(self):
        pass


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
