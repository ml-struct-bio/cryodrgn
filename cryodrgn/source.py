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
from collections.abc import Iterable
from concurrent import futures
import numpy as np
import pandas as pd
from typing import List, Optional, Union
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


class MRCFileSource(ImageSource):
    """An ImageSource that reads an .mrc/.mrcs particle stack."""

    def __init__(
        self, filepath: str, lazy: bool = True, indices: Optional[np.ndarray] = None
    ):
        from cryodrgn.mrc import MRCHeader

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


class StarfileSource(_MRCDataFrameSource):
    def __init__(
        self,
        filepath: str,
        datadir: str = "",
        lazy: bool = True,
        indices: Optional[np.ndarray] = None,
        max_threads: int = 1,
    ):
        from cryodrgn.starfile import Starfile

        df = Starfile.load(filepath).df
        df[["__mrc_index", "__mrc_filename"]] = df["_rlnImageName"].str.split(
            "@", n=1, expand=True
        )
        df["__mrc_index"] = pd.to_numeric(df["__mrc_index"]) - 1

        if not datadir:
            datadir = os.path.dirname(filepath)

        super().__init__(
            df=df,
            datadir=os.path.abspath(datadir),
            lazy=lazy,
            indices=indices,
            max_threads=max_threads,
        )


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
