import os.path
from collections.abc import Iterable
from concurrent import futures
import numpy as np
import pandas as pd
from typing import Union, List, Optional

import torch.utils.data
from cryodrgn.starfile import Starfile
from cryodrgn.mrc import MRCHeader


class LazyImage:
    def __init__(self, src, item):
        self.src = src
        self.item = item
        self.shape = src.shape[0], src.shape[1]

    def get(self):
        return self.src.images(self.item)


class ImageSource:
    @staticmethod
    def from_file(filepath: str, *args, **kwargs):
        ext = os.path.splitext(filepath)[-1][1:]
        assert ext in ("star", "mrcs", "txt", "cs"), f"Unknown file extension {ext}"
        source = getattr(ImageSource, f"from_{ext}")(filepath, *args, **kwargs)

        lazy = kwargs.get("lazy", False)
        if lazy:
            return source
        else:
            return source.images(slice(0, None))

    @staticmethod
    def from_star(filepath: str, *args, **kwargs):
        return StarfileSource(filepath, *args, **kwargs)

    @staticmethod
    def from_mrcs(filepath: str, *args, **kwargs):
        return MRCFileSource(filepath, *args, **kwargs)

    @staticmethod
    def from_txt(filepath: str, *args, **kwargs):
        return TxtFileSource(filepath, *args, **kwargs)

    @staticmethod
    def from_cs(filepath: str, *args, **kwargs):
        return CsSource(filepath, *args, **kwargs)

    def __init__(
        self,
        L: int,
        n: int,
        filenames: Union[List[str], str, None] = None,
        n_workers: int = 1,
        dtype: str = "float32",
        lazy: bool = True,  # remove after old API is gone - everything is lazy
        preallocated: bool = False,
        **kwargs
    ):
        self.L = L + int(preallocated)
        self.n = n
        self.shape = self.L, self.L, self.n
        self.lazy = lazy

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

        self.n_workers = n_workers
        self.dtype = dtype

    def __len__(self):
        return self.n

    def __getitem__(self, item):
        return LazyImage(src=self, item=item)

    def images(self, indices: Union[int, np.ndarray, slice]) -> np.ndarray:
        if np.isscalar(indices):
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
            indices = np.array(np.arange(start, stop, step))
        else:
            raise TypeError("Unsupported Type for indices")

        if np.any(indices >= self.n):
            raise ValueError(f"indices should be < {self.n}")

        images = self._images(indices).astype(self.dtype)
        if images.ndim == 3 and images.shape[0] == 1:
            images = images.squeeze(axis=0)

        return images

    def _images(self, indices: np.ndarray):
        raise NotImplementedError("Subclasses must implement this")


class MRCFileSource(ImageSource):
    def __init__(self, filepath: str, *args, **kwargs):
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
        self.stride = self.dtype().itemsize * self.ny * self.nx

        super().__init__(L=self.ny, n=self.nz, filenames=filepath, *args, **kwargs)

    def _images(self, indices: np.ndarray, data: Optional[np.ndarray]) -> np.ndarray:
        # TODO: Load array in contiguous chunks instead of using an iterative loop
        with open(self.mrcfile_path) as f:
            if data is None:
                data = np.empty((len(indices), self.ny, self.nx), dtype=self.dtype)
            for i, index in enumerate(indices):
                # TODO: Do this concurrently in several threads
                f.seek(
                    self.start
                )  # TODO: Do not seek in every iteration - find a better way!
                count = self.ny * self.nx
                offset = index * self.stride
                _data = np.fromfile(
                    f, dtype=self.dtype, count=count, offset=offset
                ).reshape(self.ny, self.nx)
                data[i, : self.ny, : self.nx] = _data
            return data


class TxtFileSource(ImageSource):
    def __init__(self, filepath: str, n_workers: int = 1, *args, **kwargs):

        _paths = []
        for line in open(filepath).readlines():
            line = line.strip()
            if not os.path.isabs(line):
                _paths.append(os.path.join(os.path.dirname(filepath), line))
            else:
                _paths.append(line)
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

        super().__init__(L=self.ny, n=self.nz, n_workers=n_workers, *args, **kwargs)

    def _images(self, indices: np.ndarray) -> np.ndarray:
        def load_single_mrcs(
            data: np.ndarray,
            source: MRCFileSource,
            source_start_index: int,
            indices_for_source: np.ndarray,
        ):
            # return source._images(indices[indices_for_source] - source_start_index)
            # data[indices_for_source] = source._images(indices[indices_for_source] - source_start_index, data=data)
            source._images(indices[indices_for_source] - source_start_index, data=data)

        data = np.empty((len(indices), self.L, self.L), dtype=self.dtype)

        with futures.ThreadPoolExecutor(self.n_workers) as executor:
            to_do = {}
            for source_i, (source_start_index, source_end_index) in enumerate(
                self.source_intervals
            ):
                indices_for_source = np.where(
                    (source_start_index <= indices) & (indices < source_end_index)
                )[0]
                if indices_for_source.size > 0:
                    future = executor.submit(
                        load_single_mrcs,
                        data,
                        self.sources[source_i],
                        source_start_index,
                        indices_for_source,
                    )
                    to_do[future] = indices_for_source

            for future in futures.as_completed(to_do):
                # data[to_do[future]] = future.result()
                pass

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
        L = MRCFileSource(self.df["__mrc_filepath"][0]).L
        super().__init__(
            L=L,
            n=len(self.df),
            filenames=df["__mrc_filename"],
            n_workers=n_workers,
            *args,
            **kwargs,
        )

    def _images(self, indices: np.ndarray) -> np.ndarray:
        def load_single_mrcs(filepath, df):
            arr = MRCFileSource(filepath)[df["__mrc_index"]]
            # df.index indicates the positions where the data needs to be inserted -> return for use by caller
            return df.index, arr

        data = np.empty((len(indices), self.L, self.L), dtype=self.dtype)

        # Create a DataFrame corresponding to the indices we're interested in
        df = self.df.iloc[indices].reset_index(drop=True)
        groups = df.groupby("__mrc_filepath")
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
