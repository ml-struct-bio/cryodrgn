"""Regenerate healpy_grid.json, the cached nested HEALPix grids for Nside 2-128.

The shipped file was made with healpy; pix2ang_nest reproduces it to within 1 ulp.
"""

import json

import numpy as np

from cryodrgn.so3_grid import pix2ang_nest

x = {}
for r in range(7):
    Nside = 2 ** (r + 1)
    Npix = 12 * Nside * Nside
    theta, phi = pix2ang_nest(Nside, np.arange(Npix))
    x[Nside] = [theta.tolist(), phi.tolist()]

with open("healpy_grid.json", "w") as hf:
    json.dump(x, hf)
