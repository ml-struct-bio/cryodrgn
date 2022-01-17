'''Lattice object'''

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from . import utils

log = utils.log

class Lattice:
    def __init__(self, D, extent=0.5, ignore_DC=True, device=None):
        assert D % 2 == 1, "Lattice size must be odd"
        x0, x1 = np.meshgrid(np.linspace(-extent, extent, D, endpoint=True), 
                             np.linspace(-extent, extent, D, endpoint=True))
        coords = np.stack([x0.ravel(),x1.ravel(),np.zeros(D**2)],1).astype(np.float32)
        self.coords = torch.tensor(coords, device=device)
        self.extent = extent
        self.D = D
        self.D2 = int(D/2)

        # todo: center should now just be 0,0; check Lattice.rotate...
        # c = 2/(D-1)*(D/2) -1 
        # self.center = torch.tensor([c,c]) # pixel coordinate for img[D/2,D/2]
        self.center = torch.tensor([0.,0.], device=device)

        self.square_mask = {}
        self.circle_mask = {}
        
        self.freqs2d = self.coords[:,0:2]/extent/2
        
        self.ignore_DC = ignore_DC
        self.device = device

    def get_downsample_coords(self, d):
        assert d % 2 == 1
        extent = self.extent * (d-1) / (self.D-1)
        x0, x1 = np.meshgrid(np.linspace(-extent, extent, d, endpoint=True), 
                             np.linspace(-extent, extent, d, endpoint=True))
        coords = np.stack([x0.ravel(),x1.ravel(),np.zeros(d**2)],1).astype(np.float32)
        return torch.tensor(coords, device=self.device)

    def get_square_lattice(self, L):
        b,e = self.D2-L, self.D2+L+1
        center_lattice = self.coords.view(self.D,self.D,3)[b:e,b:e,:].contiguous().view(-1,3)
        return center_lattice

    def get_square_mask(self, L):
        '''Return a binary mask for self.coords which restricts coordinates to a centered square lattice'''
        if L in self.square_mask:
            return self.square_mask[L]
        assert 2*L+1 <= self.D, 'Mask with size {} too large for lattice with size {}'.format(L,self.D)
        log('Using square lattice of size {}x{}'.format(2*L+1,2*L+1))
        b,e = self.D2-L, self.D2+L
        c1 = self.coords.view(self.D,self.D,3)[b,b]
        c2 = self.coords.view(self.D,self.D,3)[e,e]
        m1 = self.coords[:,0] >= c1[0]
        m2 = self.coords[:,0] <= c2[0]
        m3 = self.coords[:,1] >= c1[1]
        m4 = self.coords[:,1] <= c2[1]
        mask = m1*m2*m3*m4
        self.square_mask[L] = mask
        if self.ignore_DC:
            raise NotImplementedError
        return mask

    def get_circular_mask(self, R):
        '''Return a binary mask for self.coords which restricts coordinates to a centered circular lattice'''
        if R in self.circle_mask:
            return self.circle_mask[R]
        assert 2*R+1 <= self.D, 'Mask with radius {} too large for lattice with size {}'.format(R,self.D)
        log('Using circular lattice with radius {}'.format(R))
        r = R/(self.D//2)*self.extent
        mask = self.coords.pow(2).sum(-1) <= r**2
        if self.ignore_DC:
            assert self.coords[self.D**2//2].sum() == 0.0
            mask[self.D**2//2] = 0
        self.circle_mask[R] = mask
        return mask

    def rotate(self, images, theta):
        '''
        images: BxYxX
        theta: Q, in radians
        '''
        images = images.expand(len(theta), *images.shape) # QxBxYxX
        cos = torch.cos(theta)
        sin = torch.sin(theta)
        rot = torch.stack([cos, sin, -sin, cos], 1).view(-1, 2, 2)
        grid = self.coords[:,0:2]/self.extent @ rot # grid between -1 and 1
        grid = grid.view(len(rot), self.D, self.D, 2) # QxYxXx2
        offset = self.center - grid[:,self.D2,self.D2] # Qx2
        grid += offset[:,None,None,:]
        rotated = F.grid_sample(images, grid) # QxBxYxX
        return rotated.transpose(0,1) # BxQxYxX

    def translate_ft(self, img, t, mask=None):
        '''
        Translate an image by phase shifting its Fourier transform
        
        Inputs:
            img: FT of image (B x img_dims x 2)
            t: shift in pixels (B x T x 2)
            mask: Mask for lattice coords (img_dims x 1)

        Returns:
            Shifted images (B x T x img_dims x 2) 

        img_dims can either be 2D or 1D (unraveled image) 
        '''
        # F'(k) = exp(-2*pi*k*x0)*F(k)
        coords = self.freqs2d if mask is None else self.freqs2d[mask]
        img = img.unsqueeze(1) # Bx1xNx2
        t = t.unsqueeze(-1) # BxTx2x1 to be able to do bmm
        tfilt = coords @ t * -2 * np.pi # BxTxNx1
        tfilt = tfilt.squeeze(-1) # BxTxN
        c = torch.cos(tfilt) # BxTxN
        s = torch.sin(tfilt) # BxTxN
        return torch.stack([img[...,0]*c-img[...,1]*s,img[...,0]*s+img[...,1]*c],-1)

    def translate_ht(self, img, t, mask=None):
        '''
        Translate an image by phase shifting its Hartley transform
        
        Inputs:
            img: HT of image (B x img_dims)
            t: shift in pixels (B x T x 2)
            mask: Mask for lattice coords (img_dims x 1)

        Returns:
            Shifted images (B x T x img_dims) 

        img must be 1D unraveled image, symmetric around DC component
        '''
        #H'(k) = cos(2*pi*k*t0)H(k) + sin(2*pi*k*t0)H(-k)
        coords = self.freqs2d if mask is None else self.freqs2d[mask]
        center = int(len(coords)/2)
        img = img.unsqueeze(1) # Bx1xN
        t = t.unsqueeze(-1) # BxTx2x1 to be able to do bmm
        tfilt = coords @ t * 2 * np.pi # BxTxNx1
        tfilt = tfilt.squeeze(-1) # BxTxN
        c = torch.cos(tfilt) # BxTxN
        s = torch.sin(tfilt) # BxTxN
        return c*img + s*img[:,:,np.arange(len(coords)-1,-1,-1)]


class EvenLattice(Lattice):
    '''For a DxD lattice where D is even, we set D/2,D/2 pixel to the center'''
    def __init__(self, D, extent=0.5, ignore_DC=False, device=None):
        # centered and scaled xy plane, values between -1 and 1
        # endpoint=False since FT is not symmetric around origin
        assert D % 2 == 0, "Lattice size must be even"
        if ignore_DC: raise NotImplementedError
        x0, x1 = np.meshgrid(np.linspace(-1, 1, D, endpoint=False), 
                             np.linspace(-1, 1, D, endpoint=False))
        coords = np.stack([x0.ravel(),x1.ravel(),np.zeros(D**2)],1).astype(np.float32)
        self.coords = torch.tensor(coords, device=device)
        self.extent = extent
        self.D = D
        self.D2 = int(D/2)

        c = 2/(D-1)*(D/2) -1 
        self.center = torch.tensor([c,c], device=device) # pixel coordinate for img[D/2,D/2]
        
        self.square_mask = {}
        self.circle_mask = {}

        self.ignore_DC = ignore_DC
        self.device = device

    def get_downsampled_coords(self, d):
        raise NotImplementedError
