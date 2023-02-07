# type: ignore

import numpy as np
import multiprocessing as mp
import os
import time
from multiprocessing import Pool
import logging
import torch
from torch.utils import data
from cryodrgn import fft, mrc, starfile, utils

logger = logging.getLogger(__name__)


def load_particles(mrcs_txt_star, lazy=False, datadir=None, preallocated=False):
    """
    Load particle stack from either a .mrcs file, a .star file, a .txt file containing paths to .mrcs files, or a
    cryosparc particles.cs file.

    lazy (bool): Return numpy array if True, or return list of LazyImages
    datadir (str or None): Base directory overwrite for .star or .cs file parsing
    """

    from cryodrgn.source import ImageSource

    src = ImageSource.from_file(
        mrcs_txt_star, lazy=lazy, datadir=datadir, preallocated=preallocated
    )
    if lazy:
        return src
    else:
        return src[:]


class LazyMRCData(data.Dataset):
    """
    Class representing an .mrcs stack file -- images loaded on the fly
    """

    def __init__(
        self,
        mrcfile,
        norm=None,
        keepreal=False,
        invert_data=False,
        ind=None,
        window=True,
        datadir=None,
        window_r=0.85,
        preallocated=False,
    ):
        assert not keepreal, "Not implemented error"
        particles = load_particles(mrcfile, True, datadir=datadir, preallocated=preallocated)

        if ind is not None:
            particles = [particles[x] for x in ind]
        N = len(particles)
        ny, nx = particles[0].shape
        assert ny == nx, "Images must be square"
        if preallocated:
            assert (
                (ny - 1) % 2 == 0
            ), "Image size must be even. Is this a preprocessed dataset? Use the --preprocessed flag if so."
        else:
            assert (
                ny % 2 == 0
            ), "Image size must be even. Is this a preprocessed dataset? Use the --preprocessed flag if so."
        logger.info("Loaded {} {}x{} images".format(N, ny, nx))
        self.particles = particles
        self.N = N
        self.preallocated = preallocated
        self.D = ny if preallocated else (ny + 1)  # after symmetrizing HT
        self.invert_data = invert_data
        if norm is None:
            norm = self.estimate_normalization()
        self.norm = norm
        self.window = window_mask(ny - int(preallocated), window_r, 0.99) if window else None

    def estimate_normalization(self, n=1000):
        n = min(n, self.N)
        imgs = torch.stack(
            [
                fft.ht2_center(self.particles[i][:-int(self.preallocated), :-int(self.preallocated)])
                for i in range(0, self.N, self.N // n)
            ]
        )

        if self.invert_data:
            imgs *= -1

        imgs = fft.symmetrize_ht(imgs, preallocated=False)
        norm = [torch.mean(imgs), torch.std(imgs)]
        norm[0] = 0
        logger.info("Normalizing HT by {} +/- {}".format(*norm))
        return norm

    def get(self, i):
        img = self.particles[i]
        if img.ndim == 2:
            img = img[np.newaxis, ...]
        if self.window is not None:
            img[..., :self.D-1, :self.D-1] *= self.window
        fft.ht2_center(img[..., :-int(self.preallocated), :-int(self.preallocated)], inplace=True)
        if self.invert_data:
            img *= -1
        fft.symmetrize_ht(img, preallocated=self.preallocated)
        img = (img - self.norm[0]) / self.norm[1]
        return img

    def __len__(self):
        return self.N

    def __getitem__(self, index):
        if type(index) == list:
            index = np.array(index)
        return self.get(index), index


def window_mask(D, in_rad, out_rad):
    assert D % 2 == 0
    x0, x1 = torch.meshgrid(
        torch.linspace(-1, 1, D+1, dtype=torch.float32)[:-1],
        torch.linspace(-1, 1, D+1, dtype=torch.float32)[:-1],
    )
    r = (x0**2 + x1**2) ** 0.5
    mask = torch.minimum(torch.tensor(1.0), torch.maximum(torch.tensor(0.0), 1 - (r - in_rad) / (out_rad - in_rad)))
    return mask


class MRCData(data.Dataset):
    """
    Class representing an .mrcs stack file
    """

    def __init__(
        self,
        mrcfile,
        norm=None,
        keepreal=False,
        invert_data=False,
        ind=None,
        window=True,
        datadir=None,
        max_threads=16,
        window_r=0.85,
        preallocated=False,
    ):
        if keepreal:
            raise NotImplementedError

        if ind is not None:
            particles = load_particles(mrcfile, True, datadir=datadir, preallocated=preallocated)
            particles = torch.tensor([particles[i].get() for i in ind])
        else:
            particles = load_particles(mrcfile, False, datadir=datadir, preallocated=preallocated)

        N, ny, nx = particles.shape
        assert ny == nx, "Images must be square"
        if preallocated:
            assert (
                (ny - 1) % 2 == 0
            ), "Image size must be even. Is this a preprocessed dataset? Use the --preprocessed flag if so."
        else:
            assert (
                ny % 2 == 0
            ), "Image size must be even. Is this a preprocessed dataset? Use the --preprocessed flag if so."
        logger.info("Loaded {} {}x{} images".format(N, ny, nx))

        self.N = N
        self.preallocated = preallocated
        self.D = ny if preallocated else (ny + 1)

        # Real space window
        if window:
            logger.info(f"Windowing images with radius {window_r}")
            particles[..., :self.D-1, :self.D-1] *= window_mask(ny - int(preallocated), window_r, 0.99)

        # compute HT
        logger.info("Computing FFT")
        max_threads = min(max_threads, mp.cpu_count())
        _start = time.perf_counter()
        fft.ht2_center(particles[..., :-int(preallocated), :-int(preallocated)], inplace=True, chunksize=1000, n_workers=max_threads)
        _end = time.perf_counter()
        logger.info(f"Converted to FFT in {_end-_start}")

        if invert_data:
            particles *= -1

        # symmetrize HT
        logger.info("Symmetrizing image data")
        particles = fft.symmetrize_ht(particles, preallocated=preallocated)

        # normalize
        if norm is None:
            norm = [0, None]

        logger.info("Normalizing image data")
        # particles = (particles - norm[0]) / norm[1]
        fft.normalize(
            particles,
            mean=norm[0],
            std=norm[1],
            inplace=True,
            chunksize=1000,
            n_workers=1,
        )

        self.particles = particles
        self.norm = norm
        self.keepreal = keepreal
        if keepreal:
            self.particles_real = particles_real  # noqa: F821
            logger.info(
                "Normalized real space images by {}".format(particles_real.std())
            )  # noqa: F821
            self.particles_real /= particles_real.std()  # noqa: F821

    def __len__(self):
        return self.N

    def __getitem__(self, index):
        if type(index) == list:
            index = np.array(index)
        return self.particles[index], index

    def get(self, index):
        return self.particles[index]


class PreprocessedMRCData(data.Dataset):
    """ """

    def __init__(self, mrcfile, norm=None, ind=None, lazy=False):
        particles = load_particles(mrcfile, lazy=lazy)
        self.lazy = lazy
        if ind is not None:
            particles = particles[ind]

        self.particles = particles
        self.N = len(particles)
        if self.lazy:
            self.D = particles[0].get().shape[0]  # ny + 1 after symmetrizing HT
        else:
            self.D = particles.shape[1]  # ny + 1 after symmetrizing HT

        logger.info(f"Loaded {len(particles)} {self.D}x{self.D} images")
        if norm is None:
            norm = list(self.calc_statistic())
            norm[0] = 0

        if not lazy:
            self.particles = (self.particles - norm[0]) / norm[1]

        logger.info("Normalized HT by {} +/- {}".format(*norm))
        self.norm = norm

    def calc_statistic(self):
        if self.lazy:
            max_size = min(10000, self.N)
            sample_index = torch.sort(
                torch.random.choice(
                    torch.arange(self.N), max(int(0.1 * self.N), max_size), replace=False
                )
            )
            print("--lazy mode, sample 10% of samples to calculate standard error...")
            data = []
            for d in sample_index:
                data.append(self.particles[d].get())
            data = torch.stack(data, 0)
            mean, std = torch.mean(data), torch.std(data)
        else:
            mean, std = torch.mean(self.particles), torch.std(self.particles)
        # print(f"std={std}, mean={mean}")
        return mean, std

    def __len__(self):
        return self.N

    def __getitem__(self, index):
        if self.lazy:
            return (self.particles[index].get() - self.norm[0]) / self.norm[1], index
        else:
            return self.particles[index], index

    def get(self, index):
        if self.lazy:
            return (self.particles[index].get() - self.norm[0]) / self.norm[1]
        else:
            return self.particles[index]


class TiltMRCData(data.Dataset):
    """
    Class representing an .mrcs tilt series pair
    """

    def __init__(
        self,
        mrcfile,
        mrcfile_tilt,
        norm=None,
        keepreal=False,
        invert_data=False,
        ind=None,
        window=True,
        datadir=None,
        window_r=0.85
    ):
        if ind is not None:
            particles_real = load_particles(mrcfile, True, datadir)
            particles_tilt_real = load_particles(mrcfile_tilt, True, datadir)
            particles_real = torch.tensor(
                [particles_real[i].get() for i in ind], dtype=torch.float32
            )
            particles_tilt_real = torch.array(
                [particles_tilt_real[i].get() for i in ind], dtype=torch.float32
            )
        else:
            particles_real = load_particles(mrcfile, False, datadir)
            particles_tilt_real = load_particles(mrcfile_tilt, False, datadir)

        N, ny, nx = particles_real.shape
        assert ny == nx, "Images must be square"
        assert (
            ny % 2 == 0
        ), "Image size must be even. Is this a preprocessed dataset? Use the --preprocessed flag if so."
        logger.info("Loaded {} {}x{} images".format(N, ny, nx))
        assert particles_tilt_real.shape == (
            N,
            ny,
            nx,
        ), "Tilt series pair must have same dimensions as untilted particles"
        logger.info("Loaded {} {}x{} tilt pair images".format(N, ny, nx))

        # Real space window
        if window:
            m = window_mask(ny, window_r, 0.99)
            particles_real *= m
            particles_tilt_real *= m

        # compute HT
        particles = torch.tensor([fft.ht2_center(img) for img in particles_real]).astype(
            torch.float32
        )
        particles_tilt = torch.tensor(
            [fft.ht2_center(img) for img in particles_tilt_real]
        ).astype(torch.float32)
        if invert_data:
            particles *= -1
            particles_tilt *= -1

        # symmetrize HT
        particles = fft.symmetrize_ht(particles)
        particles_tilt = fft.symmetrize_ht(particles_tilt)

        # normalize
        if norm is None:
            norm = [torch.mean(particles), torch.std(particles)]
            norm[0] = 0
        particles = (particles - norm[0]) / norm[1]
        particles_tilt = (particles_tilt - norm[0]) / norm[1]
        logger.info("Normalized HT by {} +/- {}".format(*norm))

        self.particles = particles
        self.particles_tilt = particles_tilt
        self.norm = norm
        self.N = N
        self.D = particles.shape[1]
        self.keepreal = keepreal
        if keepreal:
            self.particles_real = particles_real
            self.particles_tilt_real = particles_tilt_real

    def __len__(self):
        return self.N

    def __getitem__(self, index):
        return self.particles[index], self.particles_tilt[index], index

    def get(self, index):
        return self.particles[index], self.particles_tilt[index]


# TODO: LazyTilt
