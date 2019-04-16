'''Lattice object'''

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

import lie_tools
import so3_grid
import utils

log = utils.log


class Lattice:
    def __init__(self, D):
        # centered and scaled xy plane, values between -1 and 1
        # endpoint=False since FT is not symmetric around origin
        x0, x1 = np.meshgrid(np.linspace(-1, 1, D, endpoint=False), 
                             np.linspace(-1, 1, D, endpoint=False))
        coords = np.stack([x0.ravel(),x1.ravel(),np.zeros(D**2)],1).astype(np.float32)
        self.coords = torch.tensor(coords)
        self.D = D
        self.D2 = int(D/2)
        c = 2/(D-1)*(D/2) -1 
        self.center = torch.tensor([c,c]) # pixel coordinate for img[D/2,D/2]
        
        self.square_mask = {}
        self.circle_mask = {}

    def get_square_lattice(self, L):
        b,e = self.D2-L, self.D2+L+1
        center_lattice = self.coords.view(self.D,self.D,3)[b:e,b:e,:].contiguous().view(-1,3)
        return center_lattice

    def get_square_mask(self, L):
        '''Return a binary mask for self.coords which restricts coordinates to a centered square lattice'''
        if L in self.square_mask:
            return self.square_mask[L]
        assert 2*L+1 < self.D, 'Mask with size {} too large for lattice with size {}'.format(L,D)
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
        return mask

    def get_circular_mask(self, R):
        '''Return a binary mask for self.coords which restricts coordinates to a centered circular lattice'''
        if R in self.circle_mask:
            return self.circle_mask[R]
        assert 2*R+1 < self.D, 'Mask with radius {} too large for lattice with size {}'.format(R,D)
        log('Using circular lattice with radius {}'.format(R))
        r = 2*R/self.D
        mask = self.coords.pow(2).sum(-1) < r**2
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
        grid = self.coords[:,0:2] @ rot
        grid = grid.view(len(rot), self.D, self.D, 2) # QxYxXx2
        offset = self.center - grid[:,self.D2,self.D2] # Qx2
        grid += offset[:,None,None,:]
        rotated = F.grid_sample(images, grid) # QxBxYxX
        return rotated.transpose(0,1) # BxQxYxX



