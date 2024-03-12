"""Utilities for losses used by encoders."""

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


def _sqrt_positive_part(x: torch.Tensor) -> torch.Tensor:
    """
    Returns torch.sqrt(torch.max(0, x))
    but with a zero subgradient where x is 0.
    """
    ret = torch.zeros_like(x)
    positive_mask = x > 0
    ret[positive_mask] = torch.sqrt(x[positive_mask])
    return ret


def matrix_to_quaternion(matrix: torch.Tensor) -> torch.Tensor:
    """
    Convert rotations given as rotation matrices to quaternions.

    Args:
        matrix: Rotation matrices as tensor of shape (..., 3, 3).

    Returns:
        quaternions with real part first, as tensor of shape (..., 4).
    """
    if matrix.size(-1) != 3 or matrix.size(-2) != 3:
        raise ValueError(f"Invalid rotation matrix shape {matrix.shape}.")

    batch_dim = matrix.shape[:-2]
    m00, m01, m02, m10, m11, m12, m20, m21, m22 = torch.unbind(
        matrix.reshape(batch_dim + (9,)), dim=-1
    )

    q_abs = _sqrt_positive_part(
        torch.stack(
            [
                1.0 + m00 + m11 + m22,
                1.0 + m00 - m11 - m22,
                1.0 - m00 + m11 - m22,
                1.0 - m00 - m11 + m22,
            ],
            dim=-1,
        )
    )

    # we produce the desired quaternion multiplied by each of r, i, j, k
    quat_by_rijk = torch.stack(
        [
            # pyre-fixme[58]: `**` is not supported for operand types `Tensor` and
            #  `int`.
            torch.stack([q_abs[..., 0] ** 2, m21 - m12, m02 - m20, m10 - m01], dim=-1),
            # pyre-fixme[58]: `**` is not supported for operand types `Tensor` and
            #  `int`.
            torch.stack([m21 - m12, q_abs[..., 1] ** 2, m10 + m01, m02 + m20], dim=-1),
            # pyre-fixme[58]: `**` is not supported for operand types `Tensor` and
            #  `int`.
            torch.stack([m02 - m20, m10 + m01, q_abs[..., 2] ** 2, m12 + m21], dim=-1),
            # pyre-fixme[58]: `**` is not supported for operand types `Tensor` and
            #  `int`.
            torch.stack([m10 - m01, m20 + m02, m21 + m12, q_abs[..., 3] ** 2], dim=-1),
        ],
        dim=-2,
    )

    # We floor here at 0.1 but the exact level is not important; if q_abs is small,
    # the candidate won't be picked.
    flr = torch.tensor(0.1).to(dtype=q_abs.dtype, device=q_abs.device)
    quat_candidates = quat_by_rijk / (2.0 * q_abs[..., None].max(flr))

    # if not for numerical problems, quat_candidates[i] should be same (up to a sign),
    # forall i; we pick the best-conditioned one (with the largest denominator)

    return quat_candidates[
        F.one_hot(q_abs.argmax(dim=-1), num_classes=4) > 0.5, :
    ].reshape(batch_dim + (4,))


def l2_frequency_bias(y, freqs2d, mask, resolution, alpha=1):
    """
    y: [..., n_pts],
    freqs2d: [resolution ** 2, 2]
    mask: [resolution ** 2]
    resolution: int
    alpha: int
    """
    loss = torch.mean(
        torch.sum((resolution * freqs2d[mask].norm(p=2, dim=-1)) ** alpha * y**2, -1)
    )
    return loss


def kl_divergence_conf(latent_variables_dict):
    z_mu = latent_variables_dict["z"]
    z_logvar = latent_variables_dict["z_logvar"]
    kld = torch.mean(
        torch.sum((1.0 + z_logvar - z_mu.pow(2) - z_logvar.exp()) / 2.0, dim=1), dim=0
    )

    return kld


def l1_regularizer(x):
    """
    x: [batch_size, dim]

    output: [1]
    """
    return torch.mean(torch.sum(torch.abs(x), -1))
