import numpy as np
import multiprocessing as mp
import logging
import torch
from torch.utils import data
from cryodrgn import fft
from cryodrgn.source import ImageSource
from cryodrgn.utils import window_mask

logger = logging.getLogger(__name__)


class ImageDataset(data.Dataset):
    def __init__(
        self,
        mrcfile,
        tilt_mrcfile=None,
        lazy=True,
        norm=None,
        keepreal=False,
        invert_data=False,
        ind=None,
        window=True,
        datadir=None,
        window_r=0.85,
        max_threads=16,
        device: str = "cpu",
    ):
        assert ind is None, "ind not supported yet"
        assert not keepreal, "Not implemented yet"
        datadir = datadir or ""
        self.src = ImageSource.from_file(
            mrcfile, lazy=lazy, datadir=datadir, indices=ind
        )
        if tilt_mrcfile is None:
            self.tilt_src = None
        else:
            self.tilt_src = ImageSource.from_file(
                tilt_mrcfile, lazy=lazy, datadir=datadir
            )

        ny = self.src.D
        assert ny % 2 == 0, "Image size must be even."

        self.N = self.src.n
        self.D = ny + 1  # after symmetrization
        self.invert_data = invert_data
        self.window = window_mask(ny, window_r, 0.99) if window else None
        self.max_threads = min(max_threads, mp.cpu_count())
        self.norm = norm or self.estimate_normalization()
        self.device = device

    def estimate_normalization(self, n=None):
        n = min(n, self.N) if n is not None else self.N
        indices = range(0, self.N, self.N // n)
        imgs = self.src.images(indices)

        if self.window is not None:
            imgs *= self.window

        if self.max_threads > 1:
            with mp.Pool(self.max_threads) as p:
                imgs = torch.asarray(p.map(fft.ht2_center, imgs), dtype=torch.float32)
        else:
            particleslist = []
            for img in imgs:
                particleslist.append(fft.ht2_center(img))
            imgs = torch.stack(particleslist)

        if self.invert_data:
            imgs *= -1

        imgs = fft.symmetrize_ht(imgs)
        norm = (0, torch.std(imgs))
        logger.info("Normalizing HT by {} +/- {}".format(*norm))
        return norm

    def _process(self, data):
        if data.ndim == 2:
            data = data[np.newaxis, ...]
        if self.window is not None:
            data *= self.window
        data = fft.ht2_center(data)
        if self.invert_data:
            data *= -1
        data = fft.symmetrize_ht(data)
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

        if isinstance(index, int):
            logger.debug(f"ImageDataset returning images at index ({index})")
        else:
            logger.debug(
                f"ImageDataset returning images for {len(index)} indices ({index[0]}..{index[-1]})"
            )

        return particles.to(self.device), tilt.to(self.device), index
