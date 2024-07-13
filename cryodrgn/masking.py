"""Filters applied to lattice coordinates as part of training."""

import numpy as np
import torch
import logging
from cryodrgn.lattice import Lattice

logger = logging.getLogger(__name__)


class CircularMask:
    """A circular lattice coordinate filter."""

    def __init__(self, lattice: Lattice, radius: int) -> None:
        self.lattice = lattice
        self.binary_mask = self.lattice.get_circular_mask(radius)
        self.current_radius = radius

    def update_radius(self, radius: int) -> None:
        self.binary_mask = self.lattice.get_circular_mask(radius)
        self.current_radius = radius

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

    def update(self, total_images_count) -> None:
        new_radius = int(self.radius_init + total_images_count / self.add_one_every)

        if self.current_radius < new_radius <= self.radius_max:
            self.update_radius(new_radius)
            logger.info(
                f"Frequency marching mask updated, new radius = {self.current_radius}"
            )

    def update_epoch(self, n_freqs: int) -> None:
        self.update_radius(min(self.current_radius + n_freqs, self.radius_max))

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

    def update(self, total_images_count: int) -> None:
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
