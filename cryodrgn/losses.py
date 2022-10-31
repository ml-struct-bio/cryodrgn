"""Equivariance loss for Encoder"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


class EquivarianceLoss(nn.Module):
    """Equivariance loss for SO(2) subgroup."""

    def __init__(self, model, D):
        super().__init__()
        self.model = model
        self.D = D

    def forward(self, img, encoding):
        """Rotate img by a random amount and compute difference in latent encoding"""
        n = img.shape[0]
        theta = torch.rand(n) * 2 * np.pi
        img = torch.unsqueeze(img, 1)
        img_rot = self.rotate(img, theta)
        img_rot = torch.squeeze(img_rot)
        img_rot_enc = self.model.encode(img_rot)[0]
        diffs = (encoding - img_rot_enc).pow(2).view(n, -1).sum(-1)
        return diffs.mean()

    def rotate(self, img, theta):
        cos = torch.cos(theta)
        sin = torch.sin(theta)
        rotT = torch.stack([cos, sin, -sin, cos], 1).view(-1, 2, 2)
        grid = self.model.lattice.coords[:, 0:2] @ rotT
        grid = grid.view(-1, self.D, self.D, 2)
        return F.grid_sample(img, grid)
