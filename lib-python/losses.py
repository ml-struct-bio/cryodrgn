"""Equivariance loss for Encoder."""
from math import pi
import torch
import torch.nn as nn
import torch.nn.functional as F
from lie_tools import s2s1rodrigues

class EquivarianceLoss(nn.Module):
    """Equivariance Loss for SO(2) subgroup."""

    def __init__(self, model):
        super().__init__()
        self.model = model

    def forward(self, img, encoding):
        assert encoding.shape[-2:] == (3, 3), "Rotation matrix input required"
        n = img.shape[0]
        theta = torch.rand(n, device=encoding.device) * 2 * pi
        v = torch.tensor([1, 0, 0], dtype=torch.float32, device=encoding.device)
        s1 = torch.stack((torch.cos(theta), torch.sin(theta)), 1)
        g = s2s1rodrigues(expand_dim(v, n), s1)

        enc_rot = g.bmm(encoding)
        img = torch.unsqueeze(img, 1) 
        img_rot = self.rotate(img, theta)
        img_rot = torch.squeeze(img_rot)
        img_rot_enc = self.model.latent_encoder(self.model.encoder(img_rot))[0]

        diffs = (enc_rot - img_rot_enc).pow(2).view(n, -1).sum(-1)
        return diffs.mean()

    @staticmethod
    def rotate(img, theta):
        cos = torch.cos(theta)
        sin = torch.sin(theta)
        zero = torch.zeros_like(theta)
        affine = torch.stack([cos, -sin, zero, sin, cos, zero], 1).view(-1, 2, 3)
        grid = F.affine_grid(affine, img.size())
        return F.grid_sample(img, grid)

def expand_dim(x, n, dim=0):
    if dim < 0:
        dim = x.dim()+dim+1
    return x.unsqueeze(dim).expand(*[-1]*dim, n, *[-1]*(x.dim()-dim))
