"""Classes for using particle image datasets in PyTorch learning methods.

This module contains classes that implement various preprocessing and data access
methods acting on the image data stored in a cryodrgn.source.ImageSource class.
These methods are used by learning methods such as those used in volume reconstruction
algorithms; the classes are thus implemented as children of torch.utils.data.Dataset
to allow them to inherit behaviour such as batch training.

For example, during initialization, ImageDataset initializes an ImageSource class and
then also estimates normalization parameters, a non-trivial computational step. When
image data is retrieved during model training using __getitem__, the data is whitened
using these parameters.

"""
import numpy as np
from collections import Counter, OrderedDict

import logging
import torch
from typing import Optional, Tuple, Union
from scipy.spatial.transform import Rotation
from cryodrgn import fft
from cryodrgn.source import ImageSource, StarfileSource, parse_star
from cryodrgn.masking import spherical_window_mask

from torch.utils.data import DataLoader
from torch.utils.data.sampler import BatchSampler, RandomSampler, SequentialSampler

logger = logging.getLogger(__name__)


class ImageDataset(torch.utils.data.Dataset):
    def __init__(
        self,
        mrcfile,
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
        self.keepreal = keepreal
        datadir = datadir or ""
        self.ind = ind
        self.src = ImageSource.from_file(
            mrcfile,
            lazy=lazy,
            datadir=datadir,
            indices=ind,
            max_threads=max_threads,
        )
        ny = self.src.D
        assert ny % 2 == 0, "Image size must be even."
        self.N = self.src.n
        self.D = ny + 1  # after symmetrization
        self.invert_data = invert_data

        if window:
            self.window = spherical_window_mask(D=ny, in_rad=window_r, out_rad=0.99).to(
                device
            )
        else:
            self.window = None

        norm = norm or self.estimate_normalization()
        norm_real = self.estimate_normalization_real()
        self.norm_real = [float(x) for x in norm_real]
        self.norm = [float(x) for x in norm]
        self.device = device
        self.lazy = lazy

        if np.issubdtype(self.src.dtype, np.integer):
            self.window = self.window.int()

    def estimate_normalization(self, n=1000):
        n = min(n, self.N) if n is not None else self.N
        indices = range(0, self.N, self.N // n)  # FIXME: what if the data is not IID??

        imgs = torch.stack([fft.ht2_center(img) for img in self.src.images(indices)])
        if self.invert_data:
            imgs *= -1

        imgs = fft.symmetrize_ht(imgs)
        norm = (0, torch.std(imgs))
        logger.info("Normalizing HT by {} +/- {}".format(*norm))

        return norm

    def estimate_normalization_real(self, n=1000):
        n = min(n, self.N) if n is not None else self.N
        indices = range(0, self.N, self.N // n)  # FIXME: what if the data is not IID??
        imgs = self.src.images(indices)
        norm = (torch.mean(imgs), torch.std(imgs))
        logger.info("Normalized real space images by {} +/- {}".format(*norm))

        return norm

    def _process(self, data):
        if data.ndim == 2:
            data = data[np.newaxis, ...]
        if self.window is not None:
            data *= self.window
        if self.invert_data:
            data *= -1

        f_data = fft.ht2_center(data)
        f_data = fft.symmetrize_ht(f_data)
        f_data = (f_data - self.norm[0]) / self.norm[1]

        if self.keepreal:
            r_data = (data - self.norm_real[0]) / self.norm_real[1]
        else:
            r_data = None

        return r_data, f_data

    def __len__(self):
        return self.N

    def __getitem__(self, index):
        if isinstance(index, list):
            index = torch.Tensor(index).to(torch.long)

        r_particles, f_particles = self._process(self.src.images(index).to(self.device))

        # this is why it is tricky for index to be allowed to be a list!
        if r_particles is not None and len(r_particles.shape) == 2:
            r_particles = r_particles[np.newaxis, ...]
        if f_particles is not None and len(f_particles.shape) == 2:
            f_particles = f_particles[np.newaxis, ...]

        if isinstance(index, (int, np.integer)):
            logger.debug(f"ImageDataset returning images at index ({index})")
        else:
            logger.debug(
                f"ImageDataset returning images for {len(index)} indices:"
                f" ({index[0]}..{index[-1]})"
            )

        return {
            "y_real": r_particles,
            "y": f_particles,
            "r_tilt": None,
            "tilt": None,
            "index": index,
        }

    def get_slice(
        self, start: int, stop: int
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        return (
            self.src.images(slice(start, stop), require_contiguous=True).numpy(),
            None,
        )


class TiltSeriesData(ImageDataset):
    """
    Class representing tilt series
    """

    def __init__(
        self,
        tiltstar,
        ntilts=None,
        random_tilts=False,
        ind=None,
        voltage=None,
        expected_res=None,
        dose_per_tilt=None,
        angle_per_tilt=None,
        tilt_axis_angle=0.0,
        stack_tilts=False,
        **kwargs,
    ):
        # Note: ind is the indices of the *tilts*, not the particles
        super().__init__(tiltstar, ind=ind, **kwargs)

        # Parse unique particles from _rlnGroupName
        star_df, _ = parse_star(tiltstar)
        assert isinstance(self.src, StarfileSource)
        # star_df = self.src.df
        if ind is not None:
            star_df = star_df.loc[ind]

        if "_rlnGroupName" in star_df.columns:
            group_name = list(star_df["_rlnGroupName"])
        elif "_rlnGroupNumber" in star_df.columns:
            group_name = list(star_df["_rlnGroupNumber"])
        else:
            raise ValueError(
                "No tilt-series group name or number column found in star file!"
            )

        particles = OrderedDict()
        for ii, gn in enumerate(group_name):
            if gn not in particles:
                particles[gn] = []
            particles[gn].append(ii)
        self.particles = [np.asarray(pp, dtype=int) for pp in particles.values()]
        self.Np = len(particles)
        self.ctfscalefactor = np.asarray(
            star_df["_rlnCtfScalefactor"], dtype=np.float32
        )
        rank_proxy, rank_column = self._tilt_rank_proxy(star_df)
        self.tilt_numbers = np.zeros(self.N)
        for i, ind in enumerate(self.particles):
            # Rank 0 is the first image of the series (lowest dose, when known).
            # The kept stack is that order, not the order of rows in the star file.
            sort_idxs = rank_proxy[ind].argsort()
            ranks = np.empty(len(ind), dtype=int)
            ranks[sort_idxs[::-1]] = np.arange(len(ind))
            self.tilt_numbers[ind] = ranks
            self.particles[i] = ind[np.argsort(ranks)]

        self.tilt_numbers = torch.tensor(self.tilt_numbers).to(self.device)
        logger.info(f"Loaded {self.N} tilts for {self.Np} particles")
        logger.info(f"Tilt order within each particle follows {rank_column}")
        counts = Counter(group_name)
        unique_counts = set(counts.values())
        logger.info(f"{unique_counts} tilts per particle")

        self.counts = counts
        self.ntilts = ntilts or min(unique_counts)
        assert self.ntilts <= min(unique_counts)
        self.random_tilts = random_tilts
        self.voltage = voltage
        self.dose_per_tilt = dose_per_tilt
        self.stack_tilts = stack_tilts
        self.subtomogram_averaging = bool(stack_tilts)

        # Geometric tilt of dose-rank i. Positive steps come in pairs before the
        # matching negative pair: 0, +a, +2a, -a, -2a, +3a, +4a, ...
        # That is the schedule whose magnitudes match _rlnAngleTilt on the
        # purified-yeast tilt series. The Hagen alternation (0, +a, -a, +2a, -2a)
        # assigns those images to the wrong projection direction.
        self.tilt_angles = None
        self.tilt_rots = None
        self.tilt_scheme_angles = None
        if angle_per_tilt is not None:
            rank_np = self.tilt_numbers.detach().cpu().numpy().astype(int)
            full_scheme = self._paired_tilt_scheme(
                int(rank_np.max()) + 1, angle_per_tilt
            )
            self.tilt_angles = torch.tensor(
                np.abs(full_scheme)[rank_np], dtype=torch.float32, device=self.device
            )
            tilt_scheme = self._paired_tilt_scheme(self.ntilts, angle_per_tilt)
            logger.info(
                "Tilt scheme (deg, dose order): %s",
                np.array2string(np.asarray(tilt_scheme), precision=2),
            )
            tilt_rots = [
                self.tilt_rotation_matrix(tilt_axis_angle, t) for t in tilt_scheme
            ]
            self.tilt_rots = torch.tensor(np.stack(tilt_rots)).float()
            self.tilt_scheme_angles = torch.tensor(tilt_scheme).float()

    @staticmethod
    def tilt_rotation_matrix(tilt_axis_angle_deg: float, tilt_deg: float) -> np.ndarray:
        """Extrinsic ZYZ rotation: tilt axis, then stage angle, then zero.

        This is the convention used by the March 2026 drgnai trainer and by
        RELION. Intrinsic ``zyz`` agrees with it only at zero stage tilt. With
        a tilt axis near -100 degrees the two differ by about 1.5 degrees per
        degree of stage tilt, which sends pose search to a different orientation.
        """
        return Rotation.from_euler(
            "ZYZ",
            [
                tilt_axis_angle_deg * np.pi / 180.0,
                float(tilt_deg) * np.pi / 180.0,
                0.0,
            ],
        ).as_matrix()

    @staticmethod
    def _tilt_rank_proxy(star_df):
        """Value whose descending order is dose order. Lowest dose gets rank 0."""
        if "_rlnMicrographPreExposure" in star_df.columns:
            # Negate so the smallest accumulated dose sorts last in argsort and
            # therefore receives rank 0.
            proxy = -np.asarray(star_df["_rlnMicrographPreExposure"], dtype=np.float32)
            return proxy, "_rlnMicrographPreExposure (ascending dose)"
        if "_rlnCtfBfactor" in star_df.columns:
            proxy = np.asarray(star_df["_rlnCtfBfactor"], dtype=np.float32)
            return proxy, "_rlnCtfBfactor"
        if "_rlnCtfScalefactor" in star_df.columns:
            proxy = np.asarray(star_df["_rlnCtfScalefactor"], dtype=np.float32)
            return proxy, "_rlnCtfScalefactor"
        raise ValueError(
            "Cannot order tilts: the star file has none of "
            "_rlnMicrographPreExposure, _rlnCtfBfactor, or _rlnCtfScalefactor."
        )

    @staticmethod
    def _paired_tilt_scheme(n_tilts: int, angle_per_tilt: float) -> list[float]:
        """0, +a, +2a, -a, -2a, +3a, +4a, -3a, -4a, ... truncated to ``n_tilts``."""
        if n_tilts <= 0:
            return []
        tilt = [0.0]
        step = float(angle_per_tilt)
        k = step
        while len(tilt) < n_tilts:
            tilt.extend([k, k + step, -k, -(k + step)])
            k += 2.0 * step
        return tilt[:n_tilts]

    def __len__(self):
        return self.Np

    def __getitem__(self, index) -> dict[str, torch.Tensor]:
        if isinstance(index, list):
            index = torch.Tensor(index).to(torch.long)
        tilt_indices = []

        for ii in index:
            if self.random_tilts:
                tilt_index = np.random.choice(
                    self.particles[ii], self.ntilts, replace=False
                )
            else:
                # take the first ntilts
                tilt_index = self.particles[ii][0 : self.ntilts]
            tilt_indices.append(tilt_index)

        tilt_indices = np.concatenate(tilt_indices)
        r_images, f_images = self._process(
            self.src.images(tilt_indices).to(self.device)
        )
        if self.stack_tilts:
            f_images = f_images.reshape(-1, self.ntilts, *f_images.shape[-2:])
            r_images = r_images.reshape(-1, self.ntilts, *r_images.shape[-2:])

        return {
            "y": f_images,
            "y_real": r_images,
            "tilt_index": torch.as_tensor(tilt_indices, dtype=torch.long),
            "index": index,
        }

    def get_tilting_func(self):
        """Expand a particle rotation across the dose-symmetric tilt scheme."""
        if self.tilt_rots is None:
            raise ValueError(
                "angle_per_tilt is required to build the subtomogram tilt scheme."
            )

        def tilting_func(rots):
            tilts = self.tilt_rots.to(rots.device)
            return torch.sum(tilts[..., None] * rots[..., None, None, :, :], -2)

        return tilting_func

    @classmethod
    def parse_particle_tilt(
        cls, tiltstar: str
    ) -> tuple[list[np.ndarray], dict[np.int64, int]]:
        star_df, _ = parse_star(tiltstar)

        if "_rlnGroupName" in star_df.columns:
            group_name = list(star_df["_rlnGroupName"])
        elif "_rlnGroupNumber" in star_df.columns:
            group_name = list(star_df["_rlnGroupNumber"])
        else:
            raise ValueError(
                "No tilt-series group name or number column found in star file!"
            )

        particles = OrderedDict()
        for ii, gn in enumerate(group_name):
            if gn not in particles:
                particles[gn] = []
            particles[gn].append(ii)

        particles = [np.asarray(pp, dtype=int) for pp in particles.values()]
        particles_to_tilts = particles
        tilts_to_particles = {}

        for i, j in enumerate(particles):
            for jj in j:
                tilts_to_particles[jj] = i

        return particles_to_tilts, tilts_to_particles

    @classmethod
    def particles_to_tilts(
        cls, particles_to_tilts: list[np.ndarray], particles: np.ndarray
    ) -> np.ndarray:
        tilts = [particles_to_tilts[int(i)] for i in particles]
        tilts = np.concatenate(tilts)

        return tilts

    @classmethod
    def tilts_to_particles(cls, tilts_to_particles, tilts):
        particles = [tilts_to_particles[i] for i in tilts]
        particles = np.array(sorted(set(particles)))
        return particles

    def get_tilt(self, index):
        return super().__getitem__(index)

    def get_tilt_particle(self, index) -> int:
        """Get the particle index for a given tilt index."""
        for p_i, p_tilts in enumerate(self.particles):
            if index in p_tilts:
                return p_i

        return None

    def get_slice(self, start: int, stop: int) -> Tuple[np.ndarray, np.ndarray]:
        # we have to fetch all the tilts to stay contiguous, and then subset
        tilt_indices = [self.particles[index] for index in range(start, stop)]
        cat_tilt_indices = np.concatenate(tilt_indices)
        images = self.src.images(cat_tilt_indices, require_contiguous=True)

        tilt_masks = []
        for tilt_idx in tilt_indices:
            tilt_mask = np.zeros(len(tilt_idx), dtype=bool)
            if self.random_tilts:
                tilt_mask_idx = np.random.choice(
                    len(tilt_idx), self.ntilts, replace=False
                )
                tilt_mask[tilt_mask_idx] = True
            else:
                i = (len(tilt_idx) - self.ntilts) // 2
                tilt_mask[i : i + self.ntilts] = True
            tilt_masks.append(tilt_mask)
        tilt_masks = np.concatenate(tilt_masks)
        selected_images = images[tilt_masks]
        selected_tilt_indices = cat_tilt_indices[tilt_masks]

        return selected_images.numpy(), selected_tilt_indices

    def critical_exposure(self, freq):
        assert (
            self.voltage is not None
        ), "Critical exposure calculation requires voltage"

        assert (
            self.voltage == 300 or self.voltage == 200
        ), "Critical exposure calculation requires 200kV or 300kV imaging"

        # From Grant and Grigorieff, 2015
        scale_factor = 1
        if self.voltage == 200:
            scale_factor = 0.75
        critical_exp = torch.pow(freq, -1.665)
        critical_exp = torch.mul(critical_exp, scale_factor * 0.245)
        return torch.add(critical_exp, 2.81)

    def get_dose_filters(self, tilt_index, lattice, Apix, *, apply_tilt_cosine=True):
        """Grant–Grigorieff exposure filter, one weight per Fourier pixel.

        ``apply_tilt_cosine`` multiplies the filter by ``cos(alpha)``. The tilt
        VAE uses that extra scale. Subtomogram ``abinit`` does not: the March
        2026 trainer's filter is the exposure term alone, and the pose-search
        score has to use the same weights as the training residual.
        """
        D = lattice.D

        N = len(tilt_index)
        freqs = lattice.freqs2d / Apix  # D/A
        x = freqs[..., 0]
        y = freqs[..., 1]
        s2 = x**2 + y**2
        s = torch.sqrt(s2)

        cumulative_dose = self.tilt_numbers[tilt_index] * self.dose_per_tilt
        cd_tile = torch.repeat_interleave(cumulative_dose, D * D).view(N, -1)

        ce = self.critical_exposure(s).to(self.device)
        ce_tile = ce.repeat(N, 1)

        oe_tile = ce_tile * 2.51284  # Optimal exposure
        oe_mask = (cd_tile < oe_tile).long()

        freq_correction = torch.exp(-0.5 * cd_tile / ce_tile)
        freq_correction = torch.mul(freq_correction, oe_mask)
        if apply_tilt_cosine:
            angle_correction = torch.cos(self.tilt_angles[tilt_index] * np.pi / 180)
            ac_tile = torch.repeat_interleave(angle_correction, D * D).view(N, -1)
            freq_correction = torch.mul(freq_correction, ac_tile)

        return freq_correction.float()

    def optimal_exposure(self, freq):
        return 2.51284 * self.critical_exposure(freq)


class DataShuffler:
    def __init__(
        self, dataset: ImageDataset, batch_size, buffer_size, dtype=np.float32
    ):
        if not all(dataset.src.indices == np.arange(dataset.N)):
            raise NotImplementedError(
                "--ind is not supported for the data shuffler. "
                "The purpose of the shuffler is to load chunks contiguously during "
                "lazy loading on huge datasets, which doesn't work with --ind subsets. "
                "We recommend instead using --ind during preprocessing (e.g. with "
                "`cryodrgn downsample`) if you aim to use the shuffler or simply "
                "pass --lazy for on-the-fly data loading (potentially slower)."
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
        self.ntilts = getattr(dataset, "ntilts", 1)  # FIXME

    def __iter__(self):
        return _DataShufflerIterator(self)


class _DataShufflerIterator:
    def __init__(self, shuffler: DataShuffler):
        self.dataset = shuffler.dataset
        self.buffer_size = shuffler.buffer_size
        self.batch_size = shuffler.batch_size
        self.batch_capacity = shuffler.batch_capacity
        self.dtype = shuffler.dtype
        self.ntilts = shuffler.ntilts

        self.buffer = np.empty(
            (self.buffer_size, self.ntilts, self.dataset.D - 1, self.dataset.D - 1),
            dtype=self.dtype,
        )
        self.index_buffer = np.full((self.buffer_size,), -1, dtype=np.int64)
        self.tilt_index_buffer = np.full(
            (self.buffer_size, self.ntilts), -1, dtype=np.int64
        )
        self.num_batches = (
            len(self.dataset) // self.batch_size
        )  # FIXME off-by-one? Nah, lets leave the last batch behind
        self.chunk_order = torch.randperm(self.num_batches)
        self.count = 0
        self.flush_remaining = -1  # at the end of the epoch, got to flush the buffer
        # pre-fill
        logger.info("Pre-filling data shuffler buffer...")
        for i in range(self.batch_capacity):
            chunk, maybe_tilt_indices, chunk_indices = self._get_next_chunk()
            self.buffer[i * self.batch_size : (i + 1) * self.batch_size] = chunk
            self.index_buffer[
                i * self.batch_size : (i + 1) * self.batch_size
            ] = chunk_indices
            if maybe_tilt_indices is not None:
                self.tilt_index_buffer[
                    i * self.batch_size : (i + 1) * self.batch_size
                ] = maybe_tilt_indices
        logger.info(
            f"Filled buffer with {self.buffer_size} images ({self.batch_capacity} contiguous chunks)."
        )

    def _get_next_chunk(self) -> Tuple[np.ndarray, Optional[np.ndarray], np.ndarray]:
        chunk_idx = int(self.chunk_order[self.count])
        self.count += 1
        particles, maybe_tilt_indices = self.dataset.get_slice(
            chunk_idx * self.batch_size, (chunk_idx + 1) * self.batch_size
        )
        particle_indices = np.arange(
            chunk_idx * self.batch_size, (chunk_idx + 1) * self.batch_size
        )
        particles = particles.reshape(
            self.batch_size, self.ntilts, *particles.shape[1:]
        )
        if maybe_tilt_indices is not None:
            maybe_tilt_indices = maybe_tilt_indices.reshape(
                self.batch_size, self.ntilts
            )
        return particles, maybe_tilt_indices, particle_indices

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
            tilt_indices = self.tilt_index_buffer[
                self.flush_remaining - self.batch_size : self.flush_remaining
            ]
            self.flush_remaining -= self.batch_size
        else:
            indices = np.random.choice(
                self.buffer_size, size=self.batch_size, replace=False
            )
            particles = self.buffer[indices]
            particle_indices = self.index_buffer[indices]
            tilt_indices = self.tilt_index_buffer[indices]

            chunk, maybe_tilt_indices, chunk_indices = self._get_next_chunk()
            self.buffer[indices] = chunk
            self.index_buffer[indices] = chunk_indices
            if maybe_tilt_indices is not None:
                self.tilt_index_buffer[indices] = maybe_tilt_indices

        particles = torch.from_numpy(particles)
        particle_indices = torch.from_numpy(particle_indices)
        tilt_indices = torch.from_numpy(tilt_indices)

        # merge the batch and tilt dimension
        particles = particles.view(-1, *particles.shape[2:])
        tilt_indices = tilt_indices.view(-1, *tilt_indices.shape[2:])

        r_particles, f_particles = self.dataset._process(
            particles.to(self.dataset.device)
        )
        # print('ZZZ', particles.shape, tilt_indices.shape, particle_indices.shape)
        return {
            "y": f_particles,
            "y_real": r_particles,
            "tilt_index": tilt_indices,
            "index": particle_indices,
        }


def make_dataloader(
    data: ImageDataset,
    *,
    batch_size: int,
    num_workers: int = 0,
    shuffler_size: int = 0,
    shuffle: bool = True,
    seed: Optional[int] = None,
):
    if shuffler_size > 0 and shuffle:
        assert data.lazy, "Only enable a data shuffler for lazy loading"
        return DataShuffler(data, batch_size=batch_size, buffer_size=shuffler_size)
    else:
        # see https://github.com/zhonge/cryodrgn/pull/221#discussion_r1120711123
        # for discussion of why we use BatchSampler, etc.
        if shuffle:
            generator = None if seed is None else torch.Generator().manual_seed(seed)
            sampler = RandomSampler(data, generator=generator)
        else:
            sampler = SequentialSampler(data)

        return DataLoader(
            data,
            num_workers=num_workers,
            sampler=BatchSampler(sampler, batch_size=batch_size, drop_last=False),
            batch_size=None,
            multiprocessing_context="spawn" if num_workers > 0 else None,
        )
