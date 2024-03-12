"""Masking filters used in learning."""

import numpy as np
import logging

logger = logging.getLogger(__name__)


def get_circular_mask(lattice, radius):
    """
    lattice: Lattice
    radius: float (in pixels)

    output: [resolution**2]
    """
    coords = lattice.coords
    resolution = lattice.D

    assert (
        2 * radius + 1 <= resolution
    ), f"Mask with radius {radius} too large for lattice with size {resolution}"

    r = radius / (resolution // 2) * lattice.extent
    mask = coords.pow(2).sum(-1) <= r**2
    if lattice.ignore_DC:
        mask[resolution**2 // 2] = 0

    return mask.cpu()


class CircularMask:
    def __init__(self, lattice, radius):
        """
        lattice: Lattice
        radius: float

        binary_mask: [resolution**2]
        """
        self.lattice = lattice
        self.binary_mask = get_circular_mask(lattice, radius)
        self.current_radius = radius

    def update_radius(self, radius):
        """
        radius: float
        """
        self.binary_mask = get_circular_mask(self.lattice, radius)
        self.current_radius = radius

    def get_lf_submask(self):
        return get_circular_mask(self.lattice, self.current_radius // 2)[
            self.binary_mask
        ]

    def get_hf_submask(self):
        return ~self.get_lf_submask()


class FrequencyMarchingMask(CircularMask):
    def __init__(self, lattice, radius_max, radius=3, add_one_every=100000):
        """
        lattice: Lattice
        radius: float
        add_one_every: int
        exp_factor: float
        """
        super().__init__(lattice, radius)
        self.radius_max = radius_max
        self.radius_init = radius
        self.add_one_every = add_one_every
        logger.info(f"Frequency marching initialized at r = {radius}")

    def update(self, total_images_count):
        """
        total_images_count: int
        """
        new_radius = int(self.radius_init + total_images_count / self.add_one_every)

        if self.current_radius < new_radius <= self.radius_max:
            self.update_radius(new_radius)
            logger.info(f"Mask updated. New radius = {self.current_radius}")

    def update_epoch(self, n_freqs):
        """
        n_freqs: int
        """
        self.update_radius(min(self.current_radius + n_freqs, self.radius_max))

    def reset(self):
        self.update_radius(self.radius_init)


class FrequencyMarchingExpMask(FrequencyMarchingMask):
    def __init__(
        self, lattice, radius_max, radius=3, add_one_every=100000, exp_factor=0.05
    ):
        """
        lattice: Lattice
        radius: float
        add_one_every: int
        exp_factor: float
        """
        super().__init__(lattice, radius)
        self.radius_max = radius_max
        self.radius_init = radius
        self.add_one_every = add_one_every
        self.exp_factor = exp_factor
        logger.info(f"Frequency marching initialized at r = {radius}")

    def update(self, total_images_count):
        """
        total_images_count: int
        """
        new_radius = int(
            self.radius_init
            + np.exp((total_images_count / self.add_one_every) * self.exp_factor)
            - (1.0 / self.exp_factor)
        )

        if self.current_radius < new_radius <= self.radius_max:
            self.update_radius(new_radius)
            logger.info(f"Mask updated. New radius = {self.current_radius}")
