import numpy as np
import torch
from torch.utils import data

import fft
import mrc
import utils

log = utils.log

class LazyMRCData(data.Dataset):
    '''
    Class representing an .mrcs stack file -- images loaded on the fly
    '''
    def __init__(self, mrcfile, norm=None, keepreal=False, invert_data=False):
        assert not keepreal, 'Not implemented error'
        if mrcfile.endswith('.txt'):
            particles = mrc.parse_mrc_list(mrcfile, lazy=True)
        else:
            particles, _, _ = mrc.parse_mrc(mrcfile, lazy=True)
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
        log('Normalizing HT by {} +/- {}'.format(*norm))
        self.norm = norm

    def estimate_normalization(self, n=1000):
        n = min(n,self.N)
        imgs = np.asarray([fft.ht2_center(self.particles[i].get()) for i in range(0,self.N, self.N//n)])
        if self.invert_data: imgs *= -1
        imgs = fft.symmetrize_ht(imgs)
        norm = [0, np.max(np.abs(imgs))/10]
        return norm

    def get(self, i):
        img = self.particles[i].get()
        img = fft.ht2_center(img).astype(np.float32)
        if self.invert_data: img *= -1
        img = fft.symmetrize_ht(img)
        img = (img - self.norm[0])/self.norm[1]
        return img

    def __len__(self):
        return self.N

    def __getitem__(self, index):
        return self.get(index), index


class MRCData(data.Dataset):
    '''
    Class representing an .mrcs stack file
    '''
    def __init__(self, mrcfile, norm=None, keepreal=False, invert_data=False):
        if mrcfile.endswith('.txt'):
            particles_real = mrc.parse_mrc_list(mrcfile)
        else:
            particles_real, _, _ = mrc.parse_mrc(mrcfile)
        N, ny, nx = particles_real.shape
        assert ny == nx, "Images must be square"
        assert ny % 2 == 0, "Image size must be even"
        log('Loaded {} {}x{} images'.format(N, ny, nx))
        # compute HT
        particles = np.asarray([fft.ht2_center(img) for img in particles_real])
        particles = particles.astype(np.float32)
        if invert_data: particles *= -1

        # symmetrize HT
        particles = fft.symmetrize_ht(particles)

        # normalize
        if norm is None:
            norm = [0, np.max(np.abs(particles))/10]
        particles = (particles - norm[0])/norm[1]
        log('Normalized HT by {} +/- {}'.format(*norm))

        self.particles = particles
        self.N = N
        self.D = particles.shape[1] # ny + 1 after symmetrizing HT
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

    def __init__(self, mrcfile, mrcfile_tilt, norm=None, keepreal=False, invert_data=False):
        if mrcfile.endswith('.txt'):
            particles_real = mrc.parse_mrc_list(mrcfile)
        else:
            particles_real, _, _ = mrc.parse_mrc(mrcfile)
        if mrcfile_tilt.endswith('.txt'):
            particles_tilt = mrc.parse_mrc_list(mrcfile_tilt)
        else:
            particles_tilt, _, _ = mrc.parse_mrc(mrcfile_tilt)
        N, ny, nx = particles_real.shape
        assert ny == nx, "Images must be square"
        assert ny % 2 == 0, "Image size must be even"
        log('Loaded {} {}x{} images'.format(N, ny, nx))
        assert particles_tilt.shape == (N, ny, nx), "Tilt series pair must have same dimensions as untilted particles"
        log('Loaded {} {}x{} tilt pair images'.format(N, ny, nx))

        # compute HT
        particles = np.asarray([fft.ht2_center(img) for img in particles_real]).astype(np.float32)
        particles_tilt = np.asarray([fft.ht2_center(img) for img in particles_tilt]).astype(np.float32)
        if invert_data: 
            particles *= -1
            particles_tilt *= -1

        # symmetrize HT
        particles = fft.symmetrize_ht(particles)
        particles_tilt = fft.symmetrize_ht(particles_tilt)

        # normalize
        if norm is None:
            norm = [0, np.max(np.abs(particles))/10]
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

    def __len__(self):
        return self.N

    def __getitem__(self, index):
        return self.particles[index], self.particles_tilt[index], index
