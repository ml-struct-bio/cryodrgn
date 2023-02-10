import numpy as np
import multiprocessing as mp
import logging
import torch
from torch.utils import data
from cryodrgn import fft
from cryodrgn.source import ImageSource
from cryodrgn.utils import window_mask

logger = logging.getLogger(__name__)


class MyMRCData(data.Dataset):
    def __init__(
        self,
        mrcfile,
        tilt_mrcfile=None,
        lazy=True,
        norm=None,
        keepreal=False,
        invert_data=False,
        ind=None,  # TODO
        window=True,
        datadir=None,
        window_r=0.85,
        max_threads=16,
    ):
        assert ind is None, "ind not supported yet"
        assert not keepreal, "Not implemented yet"
        self.src = ImageSource.from_file(
            mrcfile, lazy=lazy, datadir=datadir, preallocated=True
        )
        if tilt_mrcfile is None:
            self.tilt_src = None
        else:
            self.tilt_src = ImageSource.from_file(
                tilt_mrcfile, lazy=lazy, datadir=datadir, preallocated=True
            )

        ny, nx = self.src.L, self.src.L
        assert (ny - 1) % 2 == 0, "Image size must be odd."

        self.N = self.src.n
        self.D = ny
        self.invert_data = invert_data
        self.window = window_mask(ny - 1, window_r, 0.99) if window else None
        self.max_threads = min(max_threads, mp.cpu_count())
        self.norm = norm or self.estimate_normalization(n=1000)

    def estimate_normalization(self, n=None):
        n = min(n, self.N) if n is not None else self.N
        indices = range(0, self.N, self.N // n)
        imgs = self.src.images(indices)

        if self.window is not None:
            imgs[..., : self.D - 1, : self.D - 1] *= self.window

        fft.ht2_center(
            imgs[..., :-1, :-1],
            inplace=True,
            chunksize=1000,
            n_workers=self.max_threads,
        )

        if self.invert_data:
            imgs *= -1

        imgs = fft.symmetrize_ht(imgs, preallocated=True)
        norm = (0, torch.std(imgs))
        logger.info("Normalizing HT by {} +/- {}".format(*norm))
        return norm

    def _process(self, data):
        if data is None:
            return None
        if data.ndim == 2:
            data = data[np.newaxis, ...]
        if self.window is not None:
            data[..., : self.D - 1, : self.D - 1] *= self.window
        fft.ht2_center(data[..., :-1, :-1], inplace=True)
        if self.invert_data:
            data *= -1
        fft.symmetrize_ht(data, preallocated=True)
        data = (data - self.norm[0]) / self.norm[1]
        return data

    def __len__(self):
        return self.N

    def __getitem__(self, index):
        if isinstance(index, list):
            index = torch.Tensor(index).to(torch.long)

        particles = self._process(self.src.images(index))
        if self.tilt_src is None:
            # If no tilt data is present because a tilt_mrcfile was not specified,
            # we simply return a reference to the particle data to avoid consuming
            # any more memory while conforming to torch.Dataset's type/shape expectations,
            # and rely on the caller to properly interpret it.
            # TODO: Find a more robust way to do this.
            tilt = particles
        else:
            tilt = self._process(self.tilt_src.images(index))

        return particles, tilt, index
