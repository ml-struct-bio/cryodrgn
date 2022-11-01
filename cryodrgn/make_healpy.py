import json

import healpy
import numpy as np

x = {}
for r in range(7):
    Nside = 2 ** (r + 1)
    Npix = 12 * Nside * Nside
    theta, phi = healpy.pix2ang(Nside, np.arange(Npix), nest=True, lonlat=False)
    x[Nside] = [theta.tolist(), phi.tolist()]

with open("healpy_grid.json", "w") as hf:
    json.dump(x, hf)
