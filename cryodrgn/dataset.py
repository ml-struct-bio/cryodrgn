# type: ignore

import numpy as np

try:
    import cupy as cp
except ImportError:
    cp = None

import multiprocessing as mp
import os
import multiprocessing
import logging
from torch.utils import data

from cryodrgn import fft, mrc, starfile


logger = logging.getLogger(__name__)


def load_particles(mrcs_txt_star, lazy=False, datadir=None):
    """
    Load particle stack from either a .mrcs file, a .star file, a .txt file containing paths to .mrcs files, or a
    cryosparc particles.cs file.

    lazy (bool): Return numpy array if True, or return list of LazyImages
    datadir (str or None): Base directory overwrite for .star or .cs file parsing
    """
    if mrcs_txt_star.endswith(".txt"):
        particles = mrc.parse_mrc_list(mrcs_txt_star, lazy=lazy)
    elif mrcs_txt_star.endswith(".star"):
        # not exactly sure what the default behavior should be for the data paths if parsing a starfile
        try:
            particles = starfile.Starfile.load(mrcs_txt_star).get_particles(
                datadir=datadir, lazy=lazy
            )
        except Exception as e:
            if datadir is None:
                datadir = os.path.dirname(
                    mrcs_txt_star
                )  # assume .mrcs files are in the same director as the starfile
                particles = starfile.Starfile.load(mrcs_txt_star).get_particles(
                    datadir=datadir, lazy=lazy
                )
            else:
                raise RuntimeError(e)
    elif mrcs_txt_star.endswith(".cs"):
        particles = starfile.csparc_get_particles(mrcs_txt_star, datadir, lazy)
    else:
        particles, _ = mrc.parse_mrc(mrcs_txt_star, lazy=lazy)

    return particles


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
        use_cupy=False,
    ):
        assert not keepreal, "Not implemented error"
        particles = load_particles(mrcfile, True, datadir=datadir)
        if ind is not None:
            particles = [particles[x] for x in ind]
        N = len(particles)
        ny, nx = particles[0].get().shape
        assert ny == nx, "Images must be square"
        assert (
            ny % 2 == 0
        ), "Image size must be even. Is this a preprocessed dataset? Use the --preprocessed flag if so."
        logger.info("Loaded {} {}x{} images".format(N, ny, nx))
        self.particles = particles
        self.N = N
        self.D = ny + 1  # after symmetrizing HT
        self.invert_data = invert_data
        self.use_cupy = use_cupy  # estimate_normalization may need access to self.use_cupy, so save it first
        if norm is None:
            norm = self.estimate_normalization()
        self.norm = [float(x) for x in norm]
        self.window = window_mask(ny, window_r, 0.99) if window else None

    def estimate_normalization(self, n=1000):
        pp = cp if (self.use_cupy and cp is not None) else np

        n = min(n, self.N)
        imgs = pp.asarray(
            [
                fft.ht2_center(self.particles[i].get())
                for i in range(0, self.N, self.N // n)
            ]
        )
        if self.invert_data:
            imgs *= -1
        imgs = fft.symmetrize_ht(imgs)
        norm = [pp.mean(imgs), pp.std(imgs)]
        norm[0] = 0
        logger.info("Normalizing HT by {} +/- {}".format(*norm))
        return norm

    def get(self, i):
        pp = cp if (self.use_cupy and cp is not None) else np

        img = self.particles[i].get()
        if self.window is not None:
            img *= self.window
        img = fft.ht2_center(img).astype(pp.float32)
        if self.invert_data:
            img *= -1
        img = fft.symmetrize_ht(img)
        img = (img - self.norm[0]) / self.norm[1]
        return img

    def __len__(self):
        return self.N

    def __getitem__(self, index):
        return self.get(index), index


def window_mask(D, in_rad, out_rad, use_cupy=False):
    pp = cp if (use_cupy and cp is not None) else np

    assert D % 2 == 0
    x0, x1 = pp.meshgrid(
        pp.linspace(-1, 1, D, endpoint=False, dtype=pp.float32),
        pp.linspace(-1, 1, D, endpoint=False, dtype=pp.float32),
    )
    r = (x0**2 + x1**2) ** 0.5
    mask = pp.minimum(1.0, pp.maximum(0.0, 1 - (r - in_rad) / (out_rad - in_rad)))
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
        use_cupy=False,
    ):
        pp = cp if (use_cupy and cp is not None) else np

        if keepreal:
            raise NotImplementedError
        if ind is not None:
            particles = load_particles(mrcfile, True, datadir=datadir)
            particles = pp.array([particles[i].get() for i in ind])
        else:
            particles = load_particles(mrcfile, False, datadir=datadir)
        N, ny, nx = particles.shape
        assert ny == nx, "Images must be square"
        assert (
            ny % 2 == 0
        ), "Image size must be even. Is this a preprocessed dataset? Use the --preprocessed flag if so."
        logger.info("Loaded {} {}x{} images".format(N, ny, nx))

        # Real space window
        if window:
            logger.info(f"Windowing images with radius {window_r}")
            particles *= window_mask(ny, window_r, 0.99)

        # compute HT
        logger.info("Computing FFT")
        max_threads = min(max_threads, mp.cpu_count())
        if max_threads > 1:
            logger.info(f"Spawning {max_threads} processes")
            context = multiprocessing.get_context("spawn")
            with context.Pool(max_threads) as p:
                particles = pp.asarray(
                    p.map(fft.ht2_center, particles), dtype=pp.float32
                )
        else:
            particles = pp.asarray(
                [fft.ht2_center(img) for img in particles], dtype=pp.float32
            )
            logger.info("Converted to FFT")

        if invert_data:
            particles *= -1

        # symmetrize HT
        logger.info("Symmetrizing image data")
        particles = fft.symmetrize_ht(particles)

        # normalize
        if norm is None:
            norm = [pp.mean(particles), pp.std(particles)]
            norm[0] = 0
        particles = (particles - norm[0]) / norm[1]
        logger.info("Normalized HT by {} +/- {}".format(*norm))

        self.particles = particles
        self.N = N
        self.D = particles.shape[1]  # ny + 1 after symmetrizing HT
        self.norm = [float(x) for x in norm]
        self.keepreal = keepreal
        self.use_cupy = use_cupy
        if keepreal:
            self.particles_real = particles_real  # noqa: F821
            logger.info(
                "Normalized real space images by {}".format(
                    particles_real.std()  # noqa: F821
                )
            )
            self.particles_real /= particles_real.std()  # noqa: F821

    def __len__(self):
        return self.N

    def __getitem__(self, index):
        return self.particles[index], index

    def get(self, index):
        return self.particles[index]


class PreprocessedMRCData(data.Dataset):
    """ """

    def __init__(self, mrcfile, norm=None, ind=None, lazy=False, use_cupy=False):
        self.use_cupy = use_cupy
        if ind is not None:
            # First lazy load to avoid loading the whole dataset
            particles = load_particles(mrcfile, True)
            if not lazy:
                # Then, load the desired particles specified by ind
                particles = np.array([particles[i].get() for i in ind])
            else:
                particles = particles[ind]
        else:
            particles = load_particles(mrcfile, lazy=lazy)
        self.lazy = lazy
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
        self.norm = [float(x) for x in norm]

    def calc_statistic(self):
        pp = cp if (self.use_cupy and cp is not None) else np

        if self.lazy:
            max_size = min(10000, self.N)
            sample_index = pp.sort(
                pp.random.choice(
                    pp.arange(self.N), max(int(0.1 * self.N), max_size), replace=False
                )
            )
            print("--lazy mode, sample 10% of samples to calculate standard error...")
            data = []
            for d in sample_index:
                data.append(self.particles[d].get())
            data = pp.stack(data, 0)
            mean, std = pp.mean(data), pp.std(data)
        else:
            mean, std = pp.mean(self.particles), pp.std(self.particles)
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
        window_r=0.85,
        use_cupy=False,
    ):
        pp = cp if (use_cupy and cp is not None) else np

        if ind is not None:
            particles_real = load_particles(mrcfile, True, datadir)
            particles_tilt_real = load_particles(mrcfile_tilt, True, datadir)
            particles_real = pp.array(
                [particles_real[i].get() for i in ind], dtype=pp.float32
            )
            particles_tilt_real = pp.array(
                [particles_tilt_real[i].get() for i in ind], dtype=pp.float32
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
        particles = pp.asarray([fft.ht2_center(img) for img in particles_real]).astype(
            pp.float32
        )
        particles_tilt = pp.asarray(
            [fft.ht2_center(img) for img in particles_tilt_real]
        ).astype(pp.float32)
        if invert_data:
            particles *= -1
            particles_tilt *= -1

        # symmetrize HT
        particles = fft.symmetrize_ht(particles)
        particles_tilt = fft.symmetrize_ht(particles_tilt)

        # normalize
        if norm is None:
            norm = [pp.mean(particles), pp.std(particles)]
            norm[0] = 0
        particles = (particles - norm[0]) / norm[1]
        particles_tilt = (particles_tilt - norm[0]) / norm[1]
        logger.info("Normalized HT by {} +/- {}".format(*norm))

        self.particles = particles
        self.particles_tilt = particles_tilt
        self.norm = [float(x) for x in norm]
        self.N = N
        self.D = particles.shape[1]
        self.keepreal = keepreal
        self.use_cupy = use_cupy
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
