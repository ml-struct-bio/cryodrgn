import numpy as np
import torch
from torch.utils import data

import fft
import mrc

class LazyMRCData(data.Dataset):
    def __init__(self, *args):
        raise NotImplementedError

class MRCData(data.Dataset):
    '''
    Class representing an .mrcs stack file
    '''
    def __init__(self, mrcfile, norm=None, keepreal=False):
        particles_real, _, _ = mrc.parse_mrc(mrcfile)
        N, ny, nx = particles_real.shape
        assert ny == nx, "Images must be square"
        assert ny % 2 == 0, "Image size must be even"
        # compute FT
        particles = np.asarray([fft.fft2_center(img) for img in particles_real])
        particles = np.stack([particles.real.astype(np.float32), particles.imag.astype(np.float32)], -1)
        # normalize
        if norm is None:
            norm  = [np.mean(particles), np.std(particles)]
            norm[0] = 0
        particles = (particles - norm[0])/norm[1]
        self.particles = particles
        self.N = N
        self.D = ny
        self.norm = norm
        self.keepreal = keepreal
        if keepreal:
            self.particles_real = particles_real

    def __len__(self):
        return self.N

    def __getitem__(self, index):
        return self.particles[index], index

class TiltMRCData(data.Dataset):
    '''
    Class representing an .mrcs tilt series pair
    '''

    def __init__(self, mrcfile, mrcfile_tilt, norm=None, keepreal=False):
        particles_real, _, _ = mrc.parse_mrc(mrcfile)
        N, ny, nx = particles_real.shape
        assert ny == nx, "Images must be square"
        assert ny % 2 == 0, "Image size must be even"
        particles_tilt, _, _ = mrc.parse_mrc(mrcfile_tilt)
        assert particles_tilt.shape == (N, ny, nx), "Tilt series pair must have same dimensions as untilted particles"
        # compute FT
        particles = np.asarray([fft.fft2_center(img) for img in particles_real])
        particles = np.stack([particles.real.astype(np.float32), particles.imag.astype(np.float32)], -1)
        particles_tilt = np.asarray([fft.fft2_center(img) for img in particles_tilt])
        particles_tilt = np.stack([particles_tilt.real.astype(np.float32), particles_tilt.imag.astype(np.float32)], -1)
        # normalize
        if norm is None:
            norm  = [np.mean(particles), np.std(particles)]
            norm[0] = 0
        particles = (particles - norm[0])/norm[1]
        particles_tilt = (particles_tilt - norm[0])/norm[1]
        self.particles = particles
        self.particles_tilt = particles_tilt
        self.norm = norm
        self.N = N
        self.D = ny
        self.keepreal = keepreal
        if keepreal:
            self.particles_real = particles_real

    def __len__(self):
        return self.N

    def __getitem__(self, index):
        return self.particles[index], self.particles_tilt[index], index



   

