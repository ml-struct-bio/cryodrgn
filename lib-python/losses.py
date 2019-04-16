"""Equivariance loss for Encoder."""

from math import pi
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from lie_tools import s2s1rodrigues

class EquivarianceLoss(nn.Module):
    """Equivariance Loss for SO(2) subgroup."""

    def __init__(self, model, ny, nx):
        super().__init__()
        self.model = model
        self.ny = ny
        self.nx = nx

    def forward(self, img, encoding):
        n = img.shape[0]
        theta = torch.rand(n) * 2 * pi

        img = torch.unsqueeze(img, 1) 
        img_rot = self.rotate(img, theta)
        img_rot = torch.squeeze(img_rot)
        img_rot_enc = self.model.encode(img_rot)[0]

        diffs = (encoding - img_rot_enc).pow(2).view(n, -1).sum(-1)
        #import pickle
        #local = locals()
        #debug = {kk:local[kk] for kk in ['diffs', 'img_rot_enc', 'img_rot', 'enc_rot', 'g', 's1', 'v', 'theta', 'n', 'encoding', 'img']}
        #pickle.dump(debug, open('debug_equiv.pkl','wb'))
        #sys.exit(1)
        return diffs.mean()

    def rotate(self, img, theta):
        cos = torch.cos(theta)
        sin = torch.sin(theta)
        rotT = torch.stack([cos, sin, -sin, cos], 1).view(-1, 2, 2)
        grid = self.model.lattice[:,0:2] @ rotT # get x and y dimensions of lattice (no z)
        grid = grid.view(-1, self.ny, self.nx, 2)
        return F.grid_sample(img, grid)

def expand_dim(x, n, dim=0):
    if dim < 0:
        dim = x.dim()+dim+1
    return x.unsqueeze(dim).expand(*[-1]*dim, n, *[-1]*(x.dim()-dim))
