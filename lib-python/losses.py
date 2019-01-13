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
        x0, x1 = np.meshgrid(np.linspace(-1, 1, nx, endpoint=False), # FT is not symmetric around origin
                             np.linspace(-1, 1, ny, endpoint=False))
        lattice = np.stack([x0.ravel(),x1.ravel()],1).astype(np.float32)
        self.lattice = torch.from_numpy(lattice)
        self.ny = ny
        self.nx = nx

    def forward(self, img, encoding):
        assert encoding.shape[-2:] == (3, 3), "Rotation matrix input required"
        n = img.shape[0]
        theta = torch.rand(n, device=encoding.device) * 2 * pi
        v = torch.tensor([0, 0, 1], dtype=torch.float32, device=encoding.device)
        s1 = torch.stack((torch.cos(theta), torch.sin(theta)), 1)
        g = s2s1rodrigues(expand_dim(v, n), s1)

        enc_rot = g.bmm(encoding)
        img = torch.unsqueeze(img, 1) 
        img_rot = self.rotate(img, theta)
        img_rot = torch.squeeze(img_rot)
        img_rot_enc = self.model.latent_encoder(self.model.encoder(img_rot))[0]

        diffs = (enc_rot - img_rot_enc).pow(2).view(n, -1).sum(-1)
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
        grid = self.lattice @ rotT
        grid = grid.view(-1, self.ny, self.nx, 2)
        return F.grid_sample(img, grid)

def expand_dim(x, n, dim=0):
    if dim < 0:
        dim = x.dim()+dim+1
    return x.unsqueeze(dim).expand(*[-1]*dim, n, *[-1]*(x.dim()-dim))
