import numpy as np
import torch
from torch.utils import data
import os

from . import fft
from . import mrc
from . import utils
from . import starfile

log = utils.log

def load_particles(mrcs_txt_star, lazy=False, datadir=None, relion31=False):
    '''
    Load particle stack from either a .mrcs file, a .star file, a .txt file containing paths to .mrcs files, or a cryosparc particles.cs file

    lazy (bool): Return numpy array if True, or return list of LazyImages
    datadir (str or None): Base directory overwrite for .star or .cs file parsing
    '''
    if mrcs_txt_star.endswith('.txt'):
        particles = mrc.parse_mrc_list(mrcs_txt_star, lazy=lazy)
    elif mrcs_txt_star.endswith('.star'):
        # not exactly sure what the default behavior should be for the data paths if parsing a starfile
        try:
            particles = starfile.Starfile.load(mrcs_txt_star, relion31=relion31).get_particles(datadir=datadir, lazy=lazy)
        except Exception as e:
            if datadir is None:
                datadir = os.path.dirname(mrcs_txt_star) # assume .mrcs files are in the same director as the starfile
                particles = starfile.Starfile.load(mrcs_txt_star, relion31=relion31).get_particles(datadir=datadir, lazy=lazy)
            else: raise RuntimeError(e)
    elif mrcs_txt_star.endswith('.cs'):
        particles = starfile.csparc_get_particles(mrcs_txt_star, datadir, lazy)
    else:
        particles, _ = mrc.parse_mrc(mrcs_txt_star, lazy=lazy)
    return particles


class LazyMRCData(data.Dataset):
    '''
    Class representing an .mrcs stack file -- images loaded on the fly
    '''
    def __init__(self, mrcfile, norm=None, keepreal=False, invert_data=False, ind=None, window=True, datadir=None, relion31=False):
        assert not keepreal, 'Not implemented error'
        particles = load_particles(mrcfile, True, datadir, relion31=False)
        if ind is not None:
            particles = [particles[x] for x in ind]
        N = len(particles)
        ny, nx = particles[0].get().shape
        assert ny == nx, "Images must be square"
        assert ny % 2 == 0, "Image size must be even"
        log('Loaded {} {}x{} images'.format(N, ny, nx))
        self.particles = particles
        self.N = N
        self.D = ny + 1 # after symmetrizing HT
        self.invert_data = invert_data
        if norm is None:
            norm = self.estimate_normalization()
        self.norm = norm
        self.window = window_mask(ny, .85, .99) if window else None

    def estimate_normalization(self, n=1000):
        n = min(n,self.N)
        imgs = np.asarray([fft.ht2_center(self.particles[i].get()) for i in range(0,self.N, self.N//n)])
        if self.invert_data: imgs *= -1
        imgs = fft.symmetrize_ht(imgs)
        norm = [np.mean(imgs), np.std(imgs)]
        norm[0] = 0
        log('Normalizing HT by {} +/- {}'.format(*norm))
        return norm

    def get(self, i):
        img = self.particles[i].get()
        if self.window is not None:
            img *= self.window
        img = fft.ht2_center(img).astype(np.float32)
        if self.invert_data: img *= -1
        img = fft.symmetrize_ht(img)
        img = (img - self.norm[0])/self.norm[1]
        return img

    def __len__(self):
        return self.N

    def __getitem__(self, index):
        return self.get(index), index

def window_mask(D, in_rad, out_rad):
    assert D % 2 == 0
    x0, x1 = np.meshgrid(np.linspace(-1, 1, D, endpoint=False, dtype=np.float32), 
                         np.linspace(-1, 1, D, endpoint=False, dtype=np.float32))
    r = (x0**2 + x1**2)**.5
    mask = np.minimum(1.0, np.maximum(0.0, 1 - (r-in_rad)/(out_rad-in_rad)))
    return mask

class MRCData(data.Dataset):
    '''
    Class representing an .mrcs stack file
    '''
    def __init__(self, mrcfile, norm=None, keepreal=False, invert_data=False, ind=None, window=True, datadir=None, relion31=False):
        particles_real = load_particles(mrcfile, False, datadir, relion31=False)
        if ind is not None:
            particles_real = particles_real[ind]
        N, ny, nx = particles_real.shape
        assert ny == nx, "Images must be square"
        assert ny % 2 == 0, "Image size must be even"
        log('Loaded {} {}x{} images'.format(N, ny, nx))

        # Real space window
        if window:
            particles_real *= window_mask(ny, .85, .99)

        # compute HT
        particles = np.asarray([fft.ht2_center(img) for img in particles_real])
        particles = particles.astype(np.float32)
        if invert_data: particles *= -1

        # symmetrize HT
        particles = fft.symmetrize_ht(particles)

        # normalize
        if norm is None:
            norm  = [np.mean(particles), np.std(particles)]
            norm[0] = 0
        particles = (particles - norm[0])/norm[1]
        log('Normalized HT by {} +/- {}'.format(*norm))

        self.particles = particles
        self.N = N
        self.D = particles.shape[1] # ny + 1 after symmetrizing HT
        self.norm = norm
        self.keepreal = keepreal
        if keepreal:
            self.particles_real = particles_real
            log('Normalized real space images by {}'.format(particles_real.std()))
            self.particles_real /= particles_real.std()

    def __len__(self):
        return self.N

    def __getitem__(self, index):
        return self.particles[index], index

    def get(self, index):
        return self.particles[index]

class TiltMRCData(data.Dataset):
    '''
    Class representing an .mrcs tilt series pair
    '''
    def __init__(self, mrcfile, mrcfile_tilt, norm=None, keepreal=False, invert_data=False, ind=None, window=True, datadir=None):
        # load untilted
        particles_real = load_particles(mrcfile, False, datadir)

        # load tilt series
        particles_tilt_real = load_particles(mrcfile_tilt, False, datadir)

        # filter dataset
        if ind is not None:
            particles_real = particles_real[ind]
            particles_tilt_real = particles_tilt_real[ind]
        N, ny, nx = particles_real.shape
        assert ny == nx, "Images must be square"
        assert ny % 2 == 0, "Image size must be even"
        log('Loaded {} {}x{} images'.format(N, ny, nx))
        assert particles_tilt_real.shape == (N, ny, nx), "Tilt series pair must have same dimensions as untilted particles"
        log('Loaded {} {}x{} tilt pair images'.format(N, ny, nx))

        # Real space window
        if window:
            m = window_mask(ny, .85, .99)
            particles_real *= m
            particles_tilt_real *= m 

        # compute HT
        particles = np.asarray([fft.ht2_center(img) for img in particles_real]).astype(np.float32)
        particles_tilt = np.asarray([fft.ht2_center(img) for img in particles_tilt_real]).astype(np.float32)
        if invert_data: 
            particles *= -1
            particles_tilt *= -1

        # symmetrize HT
        particles = fft.symmetrize_ht(particles)
        particles_tilt = fft.symmetrize_ht(particles_tilt)

        # normalize
        if norm is None:
            norm  = [np.mean(particles), np.std(particles)]
            norm[0] = 0
        particles = (particles - norm[0])/norm[1]
        particles_tilt = (particles_tilt - norm[0])/norm[1]
        log('Normalized HT by {} +/- {}'.format(*norm))

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
