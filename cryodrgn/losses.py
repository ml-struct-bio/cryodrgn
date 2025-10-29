"""Equivariance loss for Encoder"""

from __future__ import annotations
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

        return F.grid_sample(img, grid, align_corners=False)


# Loss utilities for DRGN-AI training (suffix style) #


def l2_frequency_bias(
    y: torch.Tensor,
    freqs2d: torch.Tensor,
    mask: torch.Tensor,
    resolution: int,
    alpha: int = 1,
) -> torch.Tensor:
    loss = torch.mean(
        torch.sum(
            (resolution * freqs2d[mask].norm(p=2, dim=-1)) ** alpha * y**2, dim=-1
        )
    )
    return loss


def kl_divergence_conf(latent_variables_dict: dict[str, torch.Tensor]) -> torch.Tensor:
    z_mu = latent_variables_dict["z"]
    z_logvar = latent_variables_dict["z_logvar"]
    return torch.mean(
        torch.sum((1.0 + z_logvar - z_mu.pow(2) - z_logvar.exp()) / 2.0, dim=1), dim=0
    )


def l1_regularizer(x: torch.Tensor) -> torch.Tensor:
    return torch.mean(torch.sum(torch.abs(x), dim=-1))
