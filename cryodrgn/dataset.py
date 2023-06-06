import numpy as np
import logging
import torch
from torch.utils import data
from typing import Tuple, Union
from cryodrgn import fft
from cryodrgn.source import ImageSource
from cryodrgn.utils import window_mask

from torch.utils.data import DataLoader
from torch.utils.data.sampler import BatchSampler, RandomSampler

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
        device: Union[str, torch.device] = "cpu",
    ):
        assert not keepreal, "Not implemented yet"
        datadir = datadir or ""
        self.ind = ind
        self.src = ImageSource.from_file(
            mrcfile,
            lazy=lazy,
            datadir=datadir,
            indices=ind,
            max_threads=max_threads,
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
        self.window = window_mask(ny, window_r, 0.99).to(device) if window else None
        norm = norm or self.estimate_normalization()
        self.norm = [float(x) for x in norm]
        self.device = device
        self.lazy = lazy

    def estimate_normalization(self, n=1000):
        n = min(n, self.N) if n is not None else self.N
        indices = range(0, self.N, self.N // n)  # FIXME: what if the data is not IID??
        imgs = self.src.images(indices)

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

        particles = self._process(self.src.images(index).to(self.device))
        if self.tilt_src is None:
            # If no tilt data is present because a tilt_mrcfile was not specified,
            # we simply return a reference to the particle data to avoid consuming
            # any more memory while conforming to torch.Dataset's type/shape expectations,
            # and rely on the caller to properly interpret it.
            # TODO: Find a more robust way to do this.
            tilt = particles
        else:
            tilt = self._process(self.tilt_src.images(index).to(self.device))

        if isinstance(index, int):
            logger.debug(f"ImageDataset returning images at index ({index})")
        else:
            logger.debug(
                f"ImageDataset returning images for {len(index)} indices ({index[0]}..{index[-1]})"
            )

        return particles, tilt, index

    def get_slice(self, start: int, stop: int) -> torch.Tensor:
        assert self.tilt_src is None
        return self.src.images(slice(start, stop), require_contiguous=True).numpy()


class DataShuffler:
    def __init__(
        self, dataset: ImageDataset, batch_size, buffer_size, dtype=np.float32
    ):
        if not all(dataset.src.indices == np.arange(len(dataset))):
            raise NotImplementedError(
                "Sorry dude, --ind is not supported for the data shuffler. "
                "The purpose of the shuffler is to load chunks contiguously during lazy loading on huge datasets, which doesn't work with --ind. "
                "If you really need this, maybe you should probably use `--ind` during preprocessing (e.g. cryodrgn downsample)."
            )
        self.dataset = dataset
        self.batch_size = batch_size
        self.buffer_size = buffer_size
        self.dtype = dtype
        assert self.buffer_size % self.batch_size == 0, (
            self.buffer_size,
            self.batch_size,
        )  # FIXME
        self.batch_capacity = self.buffer_size // self.batch_size
        assert self.buffer_size <= len(self.dataset), (
            self.buffer_size,
            len(self.dataset),
        )

    def __iter__(self):
        return _DataShufflerIterator(self)


class _DataShufflerIterator:
    def __init__(self, shuffler: DataShuffler):
        self.dataset = shuffler.dataset
        self.buffer_size = shuffler.buffer_size
        self.batch_size = shuffler.batch_size
        self.batch_capacity = shuffler.batch_capacity
        self.dtype = shuffler.dtype

        self.buffer = np.empty(
            (self.buffer_size, self.dataset.D - 1, self.dataset.D - 1), dtype=self.dtype
        )
        self.index_buffer = np.full((self.buffer_size,), -1, dtype=np.int64)
        self.num_batches = (
            len(self.dataset) // self.batch_size
        )  # FIXME off-by-one? Nah, lets leave the last batch behind
        self.chunk_order = torch.randperm(self.num_batches)
        self.count = 0
        self.flush_remaining = -1  # at the end of the epoch, got to flush the buffer
        # pre-fill
        logger.info("Pre-filling data shuffler buffer...")
        for i in range(self.batch_capacity):
            chunk, chunk_indices = self._get_next_chunk()
            self.buffer[i * self.batch_size : (i + 1) * self.batch_size] = chunk
            self.index_buffer[
                i * self.batch_size : (i + 1) * self.batch_size
            ] = chunk_indices
        logger.info(
            f"Filled buffer with {self.buffer_size} images ({self.batch_capacity} contiguous chunks)."
        )

    def _get_next_chunk(self) -> Tuple[torch.Tensor, np.ndarray]:
        chunk_idx = int(self.chunk_order[self.count])
        self.count += 1
        particles = self.dataset.get_slice(
            chunk_idx * self.batch_size, (chunk_idx + 1) * self.batch_size
        )
        particle_indices = np.arange(
            chunk_idx * self.batch_size, (chunk_idx + 1) * self.batch_size
        )
        return particles, particle_indices

    def __iter__(self):
        return self

    def __next__(self) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Returns a batch of images, and the indices of those images in the dataset.

        The buffer starts filled with `batch_capacity` random contiguous chunks.
        Each time a batch is requested, `batch_size` random images are selected from the buffer,
        and refilled with the next random contiguous chunk from disk.

        Once all the chunks have been fetched from disk, the buffer is randomly permuted and then
        flushed sequentially.
        """
        if self.count == self.num_batches and self.flush_remaining == -1:
            logger.info(
                "Finished fetching chunks. Flushing buffer for remaining batches..."
            )
            # since we're going to flush the buffer sequentially, we need to shuffle it first
            perm = np.random.permutation(self.buffer_size)
            self.buffer = self.buffer[perm]
            self.index_buffer = self.index_buffer[perm]
            self.flush_remaining = self.buffer_size

        if self.flush_remaining != -1:
            # we're in flush mode, just return chunks out of the buffer
            assert self.flush_remaining % self.batch_size == 0
            if self.flush_remaining == 0:
                raise StopIteration()
            particles = self.buffer[
                self.flush_remaining - self.batch_size : self.flush_remaining
            ]
            particle_indices = self.index_buffer[
                self.flush_remaining - self.batch_size : self.flush_remaining
            ]
            self.flush_remaining -= self.batch_size
        else:
            indices = np.random.choice(
                self.buffer_size, size=self.batch_size, replace=False
            )
            particles = self.buffer[indices]
            particle_indices = self.index_buffer[indices]
            chunk, chunk_indices = self._get_next_chunk()
            self.buffer[indices] = chunk
            self.index_buffer[indices] = chunk_indices

        particles = torch.from_numpy(particles)
        particle_indices = torch.from_numpy(particle_indices)
        particles = self.dataset._process(particles.to(self.dataset.device))
        return particles, particles, particle_indices


def make_dataloader(
    data: ImageDataset, *, batch_size: int, num_workers: int = 0, shuffler_size: int = 0
):
    if shuffler_size > 0:
        assert data.lazy, "Only enable a data shuffler for lazy loading"
        return DataShuffler(data, batch_size=batch_size, buffer_size=shuffler_size)
    else:
        return DataLoader(
            data,
            num_workers=num_workers,
            sampler=BatchSampler(
                RandomSampler(data), batch_size=batch_size, drop_last=False
            ),
            batch_size=None,
            multiprocessing_context="spawn" if num_workers > 0 else None,
        )
