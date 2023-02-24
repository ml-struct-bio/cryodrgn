import os.path
from collections.abc import Iterable
from concurrent import futures
import numpy as np
import pandas as pd
from typing import Union, List, Optional
import logging
import torch

logger = logging.getLogger(__name__)


class ImageSource:
    @staticmethod
    def from_file(filepath: str, *args, **kwargs):
        ext = os.path.splitext(filepath)[-1][1:]
        assert ext in (
            "star",
            "mrcs",
            "mrc",
            "txt",
            "cs",
        ), f"Unknown file extension {ext}"
        source = getattr(ImageSource, f"from_{ext}")(filepath, *args, **kwargs)
        return source

    @staticmethod
    def from_star(filepath: str, *args, **kwargs):
        return StarfileSource(filepath, *args, **kwargs)

    @staticmethod
    def from_mrcs(filepath: str, *args, **kwargs):
        return MRCFileSource(filepath, *args, **kwargs)

    @staticmethod
    def from_mrc(filepath: str, *args, **kwargs):
        return MRCFileSource(filepath, *args, **kwargs)

    @staticmethod
    def from_txt(filepath: str, *args, **kwargs):
        return TxtFileSource(filepath, *args, **kwargs)

    @staticmethod
    def from_cs(filepath: str, *args, **kwargs):
        return CsSource(filepath, *args, **kwargs)

    def __init__(
        self,
        D: int,
        n: int,
        filenames: Union[List[str], str, None] = None,
        n_workers: int = 1,
        dtype: str = "float32",
        lazy: bool = True,
        indices: Optional[np.ndarray] = None,
        **kwargs,
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
            assert len(filenames) == self.n
            filenames = list(filenames)
        self.filenames = np.array(filenames)

        self.lazy = lazy
        self.n_workers = n_workers
        self.dtype = dtype

        self.cache = None
        if not self.lazy:
            self.cache = ArraySource(self._images(self.indices))

    def __len__(self):
        return self.n

    def __getitem__(self, item):
        return self.images(item)

    def _convert_to_ndarray(self, indices) -> np.ndarray:
        if indices is None:
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

        if np.any(indices >= self.n):
            raise ValueError(f"indices should be < {self.n}")

        assert isinstance(indices, np.ndarray)
        return indices

    def images(
        self, indices: Optional[Union[int, np.ndarray, slice]] = None
    ) -> torch.Tensor:

        if self.cache:
            return self.cache._images(indices)
        else:
            indices = self._convert_to_ndarray(indices)
            # Convert incoming caller indices to indices that this ImageSource will use
            if self.indices is not None:
                indices = np.array(self.indices[indices])
            images = self._images(indices).astype(self.dtype)

        return torch.tensor(images)

    def _images(self, indices: np.ndarray) -> np.ndarray:
        raise NotImplementedError("Subclasses must implement this")

    def chunks(self, chunksize=1000):
        for i in range(0, self.n, chunksize):
            _slice = slice(i, i + chunksize)
            yield np.arange(_slice.start, _slice.stop), self.images(_slice)


class ArraySource(ImageSource):
    def __init__(self, data: np.ndarray):
        if data.ndim == 2:
            data = data[np.newaxis, ...]
        nz, ny, nx = data.shape
        assert ny == nx, "Only square arrays supported"
        self.data = torch.tensor(data)

        super().__init__(D=ny, n=nz)

    def _images(self, indices: Optional[Union[int, np.ndarray, slice]] = None):
        if isinstance(indices, np.ndarray):
            logger.warning(
                "Indexing inside an eager ImageSource with an ndarray creates a copy; consider using slice or scalar indexing instead"
            )
        if indices is None:
            indices = slice(0, self.n)
        return self.data[indices, ...]


class MRCFileSource(ImageSource):
    def __init__(self, filepath: str, *args, **kwargs):
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

        super().__init__(D=self.ny, n=self.nz, filenames=filepath, *args, **kwargs)

    def _images(
        self,
        indices: np.ndarray,
        data: Optional[np.ndarray] = None,
        tgt_indices: Optional[np.ndarray] = None,
    ):

        with open(self.mrcfile_path) as f:
            if data is None:
                data = np.empty((len(indices), self.D, self.D), dtype=self.dtype)
                assert (
                    tgt_indices is None
                ), "Target indices can only be specified when passing in a preallocated array"
                tgt_indices = np.arange(len(indices))
            else:
                if tgt_indices is not None:
                    assert len(tgt_indices) == len(
                        indices
                    ), "indices/tgt_indices length mismatch"
                else:
                    tgt_indices = np.arange(len(indices))

            assert isinstance(tgt_indices, np.ndarray)

            for (index, tgt_index) in zip(indices, tgt_indices):
                f.seek(self.start)
                offset = index * self.stride
                # 'offset' in the call below is w.r.t the current position of f
                _data = np.fromfile(
                    f, dtype=self.dtype, count=self.size, offset=offset
                ).reshape(self.ny, self.nx)
                data[tgt_index, : self.ny, : self.nx] = _data

            return data


class TxtFileSource(ImageSource):
    def __init__(self, filepath: str, n_workers: int = 1, *args, **kwargs):

        _paths = []
        filepath_dir = os.path.dirname(filepath)
        for line in open(filepath).readlines():
            path = line.strip()
            if not os.path.isabs(path):
                _paths.append(os.path.join(filepath_dir, path))
            else:
                _paths.append(path)
        self.sources = [MRCFileSource(path, *args, **kwargs) for path in _paths]

        # We'll only look at the header from the first .mrcs file, and assume that all headers are compatible
        header = self.sources[0].header
        self.ny, self.nx = (
            header.fields["ny"],
            header.fields["nx"],
        )
        assert self.ny == self.nx, "Only square images supported"

        # Calculate the total length of this source by consulting the individual constituent sources
        _source_lengths = np.array([s.n for s in self.sources])
        self.nz = sum(_source_lengths)

        # Maintain a list of (<start>, <end>) tuples indicating the positions of each source
        _s_i = np.add.accumulate(_source_lengths)
        self.source_intervals = [(0, _s_i[0])] + [
            (_s_i[i], _s_i[i + 1]) for i in range(0, len(_s_i) - 1)
        ]
        n_workers = min(n_workers, len(self.source_intervals))

        super().__init__(D=self.ny, n=self.nz, n_workers=n_workers, *args, **kwargs)

    def _images(self, indices: np.ndarray):
        def load_single_mrcs(
            data: np.ndarray,
            src: MRCFileSource,
            src_indices: np.ndarray,
            tgt_indices: np.ndarray,
        ):
            src._images(indices=src_indices, data=data, tgt_indices=tgt_indices)

        data = np.empty((len(indices), self.D, self.D), dtype=self.dtype)

        with futures.ThreadPoolExecutor(self.n_workers) as executor:
            to_do = []
            for source_i, (source_start_index, source_end_index) in enumerate(
                self.source_intervals
            ):
                tgt_indices = np.nonzero(
                    (source_start_index <= indices) & (indices < source_end_index)
                )[0]
                src_indices = indices[tgt_indices] - source_start_index
                src_images = src_indices.size

                if src_images > 0:
                    future = executor.submit(
                        load_single_mrcs,
                        data=data,
                        src=self.sources[source_i],
                        src_indices=src_indices,
                        tgt_indices=tgt_indices,
                    )
                    to_do.append(future)

            for future in futures.as_completed(to_do):
                res, exc = future.result(), future.exception()
                if exc is not None:
                    raise exc

        return data


class _MRCDataFrameSource(ImageSource):
    def __init__(
        self, df: pd.DataFrame, datadir: str = "", n_workers: int = 1, *args, **kwargs
    ):
        assert "__mrc_index" in df.columns
        assert "__mrc_filename" in df.columns
        self.df = df

        self.df["__mrc_filepath"] = self.df["__mrc_filename"].apply(
            lambda filename: os.path.join(datadir, filename)
        )

        # Peek into the first mrc file to get image size
        D = MRCFileSource(self.df["__mrc_filepath"][0]).D
        super().__init__(
            D=D,
            n=len(self.df),
            filenames=df["__mrc_filename"],
            n_workers=n_workers,
            *args,
            **kwargs,
        )

    def _images(self, indices: np.ndarray):
        def load_single_mrcs(filepath, df):
            src = MRCFileSource(filepath)
            # df.index indicates the positions where the data needs to be inserted -> return for use by caller
            return df.index, src._images(df["__mrc_index"])

        data = np.empty((len(indices), self.D, self.D), dtype=self.dtype)

        # Create a DataFrame corresponding to the indices we're interested in
        batch_df = self.df.iloc[indices].reset_index(drop=True)
        groups = batch_df.groupby("__mrc_filepath")
        n_workers = min(self.n_workers, len(groups))

        with futures.ThreadPoolExecutor(n_workers) as executor:
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
    def __init__(self, filepath: str, datadir: str = "", *args, **kwargs):
        from cryodrgn.starfile import Starfile

        df = Starfile.load(filepath).df
        df[["__mrc_index", "__mrc_filename"]] = df["_rlnImageName"].str.split(
            "@", 1, expand=True
        )
        df["__mrc_index"] = pd.to_numeric(df["__mrc_index"]) - 1

        if datadir:
            if not os.path.isabs(datadir):
                datadir = os.path.join(os.path.dirname(filepath), datadir)
        else:
            datadir = os.path.dirname(filepath)

        super().__init__(df=df, datadir=datadir, *args, **kwargs)


class CsSource(_MRCDataFrameSource):
    def __init__(self, filepath: str, datadir: str = "", *args, **kwargs):
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

        super().__init__(df=df, datadir=datadir, *args, **kwargs)
