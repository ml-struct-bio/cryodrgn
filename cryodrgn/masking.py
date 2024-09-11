"""Filters applied to lattice coordinates as part of training."""

import numpy as np
import torch
from scipy.ndimage import distance_transform_edt, binary_dilation
import logging
from typing import Optional, Union
from cryodrgn.lattice import Lattice

logger = logging.getLogger(__name__)


def spherical_window_mask(
    vol: Optional[Union[np.ndarray, torch.Tensor]] = None,
    *,
    D: Optional[int] = None,
    in_rad: float = 1.0,
    out_rad: float = 1.0,
) -> torch.Tensor:
    """Create a radial mask centered within a square image with a soft or hard edge.

    Given a volume or a volume's dimension, this function creates a masking array with
    values of 1.0 for points within `in_rad` of the image's center, values of 0.0 for
    points beyond `out_rad` of the center, and linearly-interpolated values between 0.0
    and 1.0 for points located between the two given radii.

    The default radii values create a mask circumscribed against the borders of the
    image with a hard edge.

    Arguments
    ---------
    vol:        A volume array to create a mask for.
    D:          Side length of the (square) image the mask is for.
    in_rad      Inner radius (fractional float between 0 and 1)
                inside which all values are 1.0
    out_rad     Outer radius (fractional float between 0 and 1)
                beyond which all values are 0.0

    Returns
    -------
    mask    A 2D torch.Tensor of shape (D,D) with mask values between
            0 (masked) and 1 (unmasked) inclusive.

    """
    if (vol is None) == (D is None):
        raise ValueError("Either `vol` or `D` must be specified!")
    if vol is not None:
        D = vol.shape[0]

    assert D % 2 == 0
    assert in_rad <= out_rad
    x0, x1 = torch.meshgrid(
        torch.linspace(-1, 1, D + 1, dtype=torch.float32)[:-1],
        torch.linspace(-1, 1, D + 1, dtype=torch.float32)[:-1],
        indexing="ij",
    )
    dists = (x0**2 + x1**2) ** 0.5

    # Create a mask with a hard edge which goes directly from 1.0 to 0.0
    if in_rad == out_rad:
        mask = (dists <= out_rad).float()

    # Create a mask with a soft edge between `in_rad` and `out_rad`
    else:
        mask = torch.minimum(
            torch.tensor(1.0),
            torch.maximum(torch.tensor(0.0), 1 - (dists - in_rad) / (out_rad - in_rad)),
        )

    return mask


def cosine_dilation_mask(
    vol: Union[np.ndarray, torch.Tensor],
    threshold: Optional[float] = None,
    dilation: int = 25,
    edge_dist: int = 15,
    apix: float = 1.0,
) -> np.ndarray:
    threshold = threshold or np.percentile(vol, 99.99) / 2
    logger.info(f"A/px={apix:.5g}; Threshold={threshold:.5g}")
    x = np.array(vol >= threshold).astype(bool)

    dilate_val = int(dilation // apix)
    if dilate_val:
        logger.info(f"Dilating initial vol>={threshold:3g} mask by {dilate_val} px")
        x = binary_dilation(x, iterations=dilate_val).astype(float)
    else:
        logger.info("no mask dilation applied")

    dist_val = edge_dist / apix
    logger.info(f"Width of cosine edge: {dist_val:.2f} px")
    if dist_val:
        y = distance_transform_edt(~x.astype(bool))
        y[y > dist_val] = dist_val
        z = np.cos(np.pi * y / dist_val / 2)
    else:
        z = x.astype(float)

    return z.round(6)


class CircularMask:
    """A circular lattice coordinate filter that is not updated over training."""

    def __init__(self, lattice: Lattice, radius: int) -> None:
        self.lattice = lattice
        self.binary_mask = self.lattice.get_circular_mask(radius)
        self.current_radius = radius

    def update_radius(self, radius: int) -> None:
        self.binary_mask = self.lattice.get_circular_mask(radius)
        self.current_radius = radius

    def update_batch(self, total_images_count: int) -> None:
        pass

    def update_epoch(self, n_frequencies: int) -> None:
        pass

    def get_lf_submask(self) -> torch.Tensor:
        return self.lattice.get_circular_mask(self.current_radius // 2)[
            self.binary_mask
        ]

    def get_hf_submask(self) -> torch.Tensor:
        return ~self.get_lf_submask()


class FrequencyMarchingMask(CircularMask):
    """Circular lattice coordinate filters that are broadened as training proceeds."""

    def __init__(
        self,
        lattice: Lattice,
        radius: int,
        radius_max: int,
        add_one_every: int = 100000,
    ) -> None:
        super().__init__(lattice, radius)
        self.radius_max = radius_max
        self.radius_init = radius
        self.add_one_every = add_one_every

    def update_batch(self, total_images_count) -> None:
        new_radius = int(self.radius_init + total_images_count / self.add_one_every)

        if self.current_radius < new_radius <= self.radius_max:
            self.update_radius(new_radius)
            logger.info(
                f"Frequency marching mask updated, new radius = {self.current_radius}"
            )

    def update_epoch(self, n_frequencies: int) -> None:
        self.update_radius(min(self.current_radius + n_frequencies, self.radius_max))

    def reset(self) -> None:
        self.update_radius(self.radius_init)


class FrequencyMarchingExpMask(FrequencyMarchingMask):
    def __init__(
        self,
        lattice: Lattice,
        radius: int,
        radius_max: int,
        add_one_every: int = 100000,
        exp_factor=0.05,
    ) -> None:
        super().__init__(lattice, radius, radius_max, add_one_every)
        self.exp_factor = exp_factor

    def update_batch(self, total_images_count: int) -> None:
        new_radius = int(
            self.radius_init
            + np.exp((total_images_count / self.add_one_every) * self.exp_factor)
            - (1.0 / self.exp_factor)
        )

        if self.current_radius < new_radius <= self.radius_max:
            self.update_radius(new_radius)
            logger.info(
                f"Exp. Frequency marching mask updated, "
                f"new radius = {self.current_radius}"
            )
