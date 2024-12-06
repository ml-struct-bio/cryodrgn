"""Classes for reading and using particle image data from various file formats.

This module contains the class hierarchy used by cryoDRGN for loading image stacks.
These stacks can be saved to file in a variety of formats such as
.mrcs, .star, and .txt, each of which is handled by its own class.

The base class of the hierarchy, `ImageSource`, contains the `.from_file()` method,
which accepts a path to a .star/.mrcs/.txt/.cs file, detects the file format, and
instantiates the appropriate child class.
An `images` method is used at runtime to retrieve 3D Tensors for image data
at specified indices. Chunked access is possible using the `chunks()` method.

See also
--------
`cryodrgn.dataset.ImageDataset` — Using `ImageSource` data in torch.data workflows
`cryodrgn.mrcfile`, `cryodrgn.starfile` — Utilities for specific file formats

Example usage
-------------
> from cryodrgn.source import ImageSource
> # load testing image dataset found in this repository; .mrcs is detected automatically
> src = ImageSource.from_file("tests/data/hand.mrcs", lazy=True)

> # get a 10x64x64 torch Tensor with the images at indices 10...19
> im = src.images(range(10, 20))

> # iterate through all images in the dataset, twenty images at a time
> # this is useful for avoiding memory issues when loading large datasets
> for indices, chunk in src.chunks(chunksize=20):
>    assert chunk.shape == (20, 64, 64)

"""
import os.path
from collections.abc import Iterable
from concurrent import futures
import numpy as np
import pandas as pd
from typing import List, Iterator, Optional, Union, Callable
import logging
import torch
from cryodrgn.mrcfile import MRCHeader, write_mrc, fix_mrc_header
from cryodrgn.starfile import parse_star, Starfile

logger = logging.getLogger(__name__)


class ImageSource:
    """A class that returns the 3D image data in a .mrcs/.txt/.star file stack.

    `images(<indices>)` method is used to read in images as torch Tensors.
    <indices> can be omitted to get all images at once;
    it can also be a scalar, a slice, a numpy array, or an iterable of indices
    where we want to query the data.
    Only square images are supported, of side length D pixels.
    The dimensions of the returned Tensor is (<n_images>, D, D).

    The underlying file can be loaded using the `from_file` method in lazy (default)
    mode, which is quick, but defers the actual reading of file(s) until
    images(<indices>) is called.
    In non-lazy mode, the file(s) are read immediately.

    The `images()` call always returns a copy of the data,
    whether the `ImageSource` is lazy or not.

    Arguments
    ---------
    D           Side length (pixels) of the square images in this stack.
    n           Total number of images in this stack.
    filenames   The file(s) containing the images in this stack.
    lazy        Whether to load the images in this stack immediately or on demand.
    indices     Filter the images using these indices.

    Attributes
    ----------
    orig_n (int): Total number of images originally in this stack.
    n (int): The number of images in this stack after `indices` was applied.
    shape (tuple): The shape of the underlying image data tensor - `(n, D, D)`.
    data (np.array): The image stack data loaded as a matrix.
                     Will be `None` if using lazy loading mode.
    """

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
        self.orig_n = n
        self.D = D

        # If indices is provided, it overrides self.n (but not self.orig_n)
        self.indices = np.arange(n) if indices is None else indices
        self.n = len(self.indices)
        self.shape = self.n, self.D, self.D

        # Some client calls need to access the original filename(s) associated with a
        # source; these are traditionally available as the 'fname' attribute of the
        # LazyImage class, hence only used by calling code when lazy=True
        if filenames is None:
            filenames = ["" for _ in range(self.n)]
        elif isinstance(filenames, str):
            filenames = [filenames for _ in range(self.n)]
        else:
            if len(filenames) != self.n:
                raise ValueError(f"`{len(filenames)=}` != `{n=}`")
            filenames = list(filenames)

        self.filenames = np.array(filenames)
        self.max_threads = max_threads
        self.dtype = dtype

        if lazy:
            self.data = None
        else:
            assert D is not None
            array = self._images(self.indices)
            if array.ndim == 2:
                array = array[np.newaxis, ...]

            nz, ny, nx = array.shape
            assert ny == nx, "Only square arrays supported"
            self.data = array

    @staticmethod
    def from_file(
        filepath: str,
        lazy: bool = True,
        indices: Optional[np.ndarray] = None,
        datadir: Optional[str] = None,
        max_threads: int = 1,
    ):
        """Utility for creating an `ImageSource` object from any given file type.

        This method instantiates the appropriate child class of the `ImageSource` class
        according to the extension of the given file pointing to a stack of particle
        images.

        Arguments
        ---------
        filepath:   A path to a stack of particle images.
        lazy:       Whether to load the images in this stack immediately or on demand.
        indices:    Filter the images using these indices.
        datadir:    Path prefix to particle stack files if loading relative stack
                    paths from a .star or .cs file.
        max_threads:    Maximum number of CPU cores to use for data loading.

        """
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
            raise ValueError(
                f"Unrecognized ImageSource file extension `{ext}` not in .star, "
                ".mrc/.mrcs, .txt, or .cs!"
            )

    def __len__(self) -> int:
        return self.n

    @property
    def lazy(self) -> bool:
        """Whether the images were loaded lazily or directly as a numpy array."""
        return self.data is None

    def __getitem__(self, item) -> torch.Tensor:
        return self.images(item)

    def __eq__(self, other):
        return np.allclose(self.images(), other.images())

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
            raise TypeError(f"Unsupported type `{type(indices)}` for indices")

        assert isinstance(indices, np.ndarray)
        if np.any(indices >= self.n):
            raise ValueError(f"indices should be < {self.n}")

        return indices

    def images(
        self,
        indices: Optional[Union[np.ndarray, int, slice, Iterable]] = None,
        require_contiguous: bool = False,
        as_numpy: bool = False,
    ) -> Union[np.ndarray, torch.Tensor]:
        """Retrieve the image data contained in this stack, or a subset thereof.

        Arguments
        ---------
        indices:     The subset of images we want to return.
        require_contiguous:
            Whether the method should throw an error if image retrieval will entail
            non-contiguous disk access. Callers can employ this if they insist on
            efficient loading and choose to throw an error instead of falling back
            on inefficient slower loading.
        as_numpy:    Cast image data to a numpy array (instead of PyTorch tensor).

        Returns
        -------
        images: Image data at specified indices.

        """
        indices = self._convert_to_ndarray(indices)

        if self.lazy:
            # Convert incoming caller indices to indices that this ImageSource will use
            indices = np.array(self.indices[indices])
            images = self._images(indices, require_contiguous=require_contiguous)
        else:
            images = self.data[indices, ...]  # cached data when not using lazy mode

        assert images.dtype == self.dtype, (
            f"Class `{self.__class__.__name__}` has implemented an `_images` method "
            f"that does not return arrays of numpy dtype `{self.dtype}` !"
        )
        if not as_numpy:
            images = torch.from_numpy(images)

        return images

    def _images(
        self, indices: np.ndarray, require_contiguous: bool = False
    ) -> np.ndarray:
        """Retrieve images at specified indices from this stack's files."""
        raise NotImplementedError("Subclasses of `ImageSource` must implement this!")

    def chunks(
        self,
        chunksize: int = 1000,
        indices: Optional[np.ndarray] = None,
    ) -> Iterable[tuple[np.ndarray, torch.Tensor]]:
        """A generator that returns images in chunks of size `chunksize`.

        Returns:
            A 2-tuple of (<indices>, <torch.Tensor>).

        """
        use_indices = np.arange(self.n) if indices is None else np.array(indices)
        for i in range(0, len(use_indices), chunksize):
            indx = use_indices[np.arange(i, min(len(use_indices), i + chunksize))]

            yield indx, self.images(indx)

    @property
    def apix(self) -> Union[None, float, np.ndarray]:
        """The angstroms per pixels for the images in this source."""
        return None

    def write_mrc(
        self,
        output_file: str,
        header: Optional[MRCHeader] = None,
        indices: Optional[np.ndarray] = None,
        transform_fn: Optional[Callable] = None,
        chunksize: Optional[int] = None,
    ) -> None:
        """Save this source's data to a .mrc file, using chunking if necessary."""
        if header is None and hasattr(self, "header"):
            header = self.header
        if indices is None:
            indices = np.arange(self.n)

        if chunksize is None:
            header_args = {"Apix": self.apix or 1.0} if header is None else dict()
            write_mrc(
                output_file,
                self.images(indices),
                header=header,
                transform_fn=transform_fn,
                **header_args,
            )

        else:
            if header is None:
                header = MRCHeader.make_default_header(
                    nz=len(indices),
                    ny=self.D,
                    nx=self.D,
                    dtype=self.dtype,
                    is_vol=False,
                    Apix=self.apix if self.apix is not None else 1.0,
                )
            else:
                header = fix_mrc_header(header=header)

            new_dtype = np.dtype(header.dtype).newbyteorder(header.ENDIANNESS)  # type: ignore
            with open(output_file, "wb") as f:
                header.write(f)

                for i, (indx, chunk) in enumerate(self.chunks(chunksize, indices)):
                    logger.debug(f"Processing chunk {i}")

                    if transform_fn is not None:
                        chunk = transform_fn(chunk, indx)
                    if isinstance(chunk, torch.Tensor):
                        chunk = np.array(chunk.cpu()).astype(new_dtype)

                    f.write(chunk.tobytes())

    def get_default_mrc_header(self) -> MRCHeader:
        return MRCHeader.make_default_header(data=self.images().numpy(), is_vol=False)


class MRCFileSource(ImageSource):
    """An image stack saved as a single .mrc/.mrcs file."""

    def __init__(
        self, filepath: str, lazy: bool = True, indices: Optional[np.ndarray] = None
    ):
        self.header = MRCHeader.parse(filepath)
        self.mrcfile_path = filepath
        self.dtype = self.header.dtype
        self.start = 1024 + self.header.fields["next"]  # start of image data
        orig_n, self.ny, self.nx = (
            self.header.fields["nz"],
            self.header.fields["ny"],
            self.header.fields["nx"],
        )
        assert self.ny == self.nx, "Only square images supported"
        self.size = self.ny * self.nx
        self.stride = self.dtype().itemsize * self.size

        # Adjust the header for the index filter we are applying
        if indices is not None:
            self.header.fields["nz"] = len(indices)

        super().__init__(
            D=self.ny,
            n=orig_n,
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
        if data is None:
            data = np.zeros((len(indices), self.D, self.D), dtype=self.dtype)
            if tgt_indices is not None:
                raise ValueError(
                    "Target indices can only be specified when passing "
                    "in a pre-allocated array using `data`!"
                )
            tgt_indices = np.arange(len(indices))
        else:
            if tgt_indices is not None:
                if len(tgt_indices) != len(indices):
                    raise ValueError(
                        f"`indices` ({len(indices)}) and `tgt_indices` "
                        f"({len(tgt_indices)}) length mismatch!"
                    )
            else:
                tgt_indices = np.arange(len(indices))

        assert isinstance(tgt_indices, np.ndarray)
        is_contiguous = np.all(indices == indices[0] + np.arange(len(indices)))
        if require_contiguous and not is_contiguous:
            raise ValueError("MRC indices are not contiguous!")

        with open(self.mrcfile_path) as f:
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

    def write(
        self,
        output_file: str,
        header: Optional[MRCHeader] = None,
        transform_fn: Optional[Callable] = None,
        chunksize: Optional[int] = None,
    ):
        self.write_mrc(output_file, header, transform_fn, chunksize)

    @property
    def apix(self) -> float:
        return self.header.apix


class _MRCDataFrameSource(ImageSource):
    """Base class for image stacks saved across a collection of .mrc/.mrcs files.

    These stacks use a single file as a reference to a collection of .mrc/.mrcs files
    storing the actual data. Examples include .star, .cs, and .txt files.

    Attributes
    ----------
    df (pd.DataFrame):  The table listing the constituent parts of this stack.
    datadir (str):  Optional path used by .cs and .star files to prepend to file names.
    _sources (dict[str, MRCFileSource])
        Index of the .mrc/.mrcs files in this collection; keys are the file paths
        and values are the data in each loaded lazily.
    """

    def __init__(
        self,
        df: pd.DataFrame,
        datadir: Optional[str] = None,
        lazy: bool = True,
        indices: Optional[np.ndarray] = None,
        max_threads: int = 1,
    ) -> None:
        assert "__mrc_index" in df.columns
        assert "__mrc_filename" in df.columns
        self.df = df
        self.datadir = datadir
        self.df["__mrc_filepath"] = self.df["__mrc_filename"].apply(self.parse_filename)

        self._sources = {
            filepath: MRCFileSource(filepath) if os.path.exists(filepath) else None
            for filepath in self.df["__mrc_filepath"].unique()
        }

        # Peek into the first mrc file to get image size
        D = None
        for filepath, src in self._sources.items():
            if isinstance(src, MRCFileSource):
                D = src.D
                break

        super().__init__(
            D=D,
            n=len(self.df),
            max_threads=max_threads,
            lazy=lazy,
            indices=indices,
        )

    def _images(
        self, indices: np.ndarray, require_contiguous: bool = False
    ) -> np.ndarray:
        def load_single_mrcs(filepath, df):
            src = self._sources[filepath]

            # `df.index` indicates the positions where the data needs to be inserted
            # and returned for use by caller
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

    @property
    def sources(self) -> Iterator[tuple[str, MRCFileSource]]:
        return iter(self._sources.items())

    def parse_filename(self, filename: str) -> str:
        """Get the complete path to an image stack using `self.datadir` if necessary.

        This function is used for operations such as getting the `__mrc_filepath` field
        from `__mrc_filename`. We first try appending the stack's full file path to the
        datadir; if that path does not yield an extant file we try again by appending
        just the stack's file name to the datadir before throwing an error. This second
        approach is for cases where the image stack path contains a directory structure
        leftover from processing.

        """
        newname = (
            os.path.abspath(filename) if os.path.isabs(filename) else str(filename)
        )
        if self.datadir is not None:
            if os.path.exists(newname):
                pass
            elif os.path.exists(os.path.join(self.datadir, newname)):
                newname = os.path.join(self.datadir, newname)
            elif os.path.exists(os.path.join(self.datadir, os.path.basename(newname))):
                newname = os.path.join(self.datadir, os.path.basename(newname))
            else:
                raise ValueError(
                    f"Cannot find file `{newname}` under `{self.datadir=}`!"
                )

        return newname


class CsSource(_MRCDataFrameSource):
    """Image stacks created using the cryoSPARC format, saved to a .cs file."""

    def __init__(
        self,
        filepath: str,
        datadir: str = "",
        lazy: bool = True,
        indices: Optional[np.ndarray] = None,
        max_threads: int = 1,
    ):
        csdata = np.load(filepath)
        blob_indices = csdata["blob/idx"]
        blob_paths = csdata["blob/path"].astype(str).tolist()
        assert len(blob_paths) == len(blob_indices)

        # If --datadir is not given, we assume that image stacks are in the
        # same location as this .cs file
        if not datadir:
            datadir = os.path.dirname(filepath)

        # Remove leading ">" from paths, if present
        blob_paths = [p[1:] if p.startswith(">") else p for p in blob_paths]
        df = pd.DataFrame({"__mrc_index": blob_indices, "__mrc_filename": blob_paths})

        super().__init__(
            df=df, datadir=datadir, lazy=lazy, indices=indices, max_threads=max_threads
        )


class TxtFileSource(_MRCDataFrameSource):
    """Image stacks indexed using a .txt file listing a .mrcs stack on each line.

    Note that .txt files differ from .cs and .star files in that the filenames contained
    therein are always assumed to be stated relative to the directory
    the .txt file is in; thus we don't need a --datadir.
    """

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

        _source_lengths = [MRCHeader.parse(path).N for path in _paths]
        mrc_filename, mrc_index = [], []
        for path, length in zip(_paths, _source_lengths):
            mrc_filename.extend([path] * length)
            mrc_index.append(np.arange(length))

        mrc_index = np.concatenate(mrc_index)
        df = pd.DataFrame({"__mrc_filename": mrc_filename, "__mrc_index": mrc_index})

        super().__init__(df=df, lazy=lazy, indices=indices, max_threads=max_threads)

    def write(self, output_file: str):
        """Save the list of stacks referenced in this source as its .txt version."""
        with open(output_file, "w") as f:
            f.write("\n".join(self.df["__mrc_filename"].unique()))


class StarfileSource(_MRCDataFrameSource, Starfile):
    """Image stacks indexed using a .star file in RELION3.0 or RELION3.1 format.

    In RELION3.1 format, these files will have an optics table that lists parameters
    for images grouped by optics parameters.
    See `Starfile.optics_values()` for how these optics values, such as A/px and
    resolution, are retrieved on an image-wise basis.

    Attributes
    ----------
    df (pd.DataFrame):  The primary data table in the .star file.
    data_optics (pd.Dataframe): `None` if RELION3.1
    """

    def __init__(
        self,
        filename: str,
        datadir: Optional[str] = None,
        lazy: bool = True,
        indices: Optional[np.ndarray] = None,
        max_threads: int = 1,
    ) -> None:
        sdata, data_optics = parse_star(filename)
        Starfile.__init__(self, data=sdata, data_optics=data_optics)
        self.df = None

        # If --datadir is not given, we assume that image stacks are in the
        # same location as this .star file
        if not datadir:
            datadir = os.path.dirname(filename)

        sdata[["__mrc_index", "__mrc_filename"]] = sdata["_rlnImageName"].str.split(
            "@", n=1, expand=True
        )
        sdata["__mrc_index"] = pd.to_numeric(sdata["__mrc_index"]) - 1

        _MRCDataFrameSource.__init__(
            self,
            df=sdata,
            datadir=os.path.abspath(datadir) if datadir else None,
            lazy=lazy,
            indices=indices,
            max_threads=max_threads,
        )
