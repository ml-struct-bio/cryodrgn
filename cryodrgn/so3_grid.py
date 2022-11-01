"""
Implementation of Yershova et al. "Generating uniform incremental
grids on SO(3) using the Hopf fribration"
"""

import json
import os

import numpy as np


def grid_s1(resol):
    Npix = 6 * 2**resol
    dt = 2 * np.pi / Npix
    grid = np.arange(Npix) * dt + dt / 2
    return grid


def grid_s2(resol):
    Nside = 2**resol
    Npix = 12 * Nside * Nside
    theta, phi = pix2ang(Nside, np.arange(Npix), nest=True)
    return theta, phi


def hopf_to_quat(theta, phi, psi):
    """
    Hopf coordinates to quaternions
    theta: [0,pi)
    phi: [0, 2pi)
    psi: [0, 2pi)
    """
    ct = np.cos(theta / 2)
    st = np.sin(theta / 2)
    quat = np.array(
        [
            ct * np.cos(psi / 2),
            ct * np.sin(psi / 2),
            st * np.cos(phi + psi / 2),
            st * np.sin(phi + psi / 2),
        ]
    )
    return quat.T.astype(np.float32)


def grid_SO3(resol):
    theta, phi = grid_s2(resol)
    psi = grid_s1(resol)
    quat = hopf_to_quat(
        np.repeat(theta, len(psi)),  # repeats each element by len(psi)
        np.repeat(phi, len(psi)),  # repeats each element by len(psi)
        np.tile(psi, len(theta)),
    )  # tiles the array len(theta) times
    return quat  # hmm convert to rot matrix?


def s2_grid_SO3(resol):
    theta, phi = grid_s2(resol)
    quat = hopf_to_quat(theta, phi, np.zeros((len(phi),)))
    return quat


# Neighbor finding


def get_s1_neighbor(mini, curr_res):
    """
    Return the 2 nearest neighbors on S1 at the next resolution level
    """
    Npix = 6 * 2 ** (curr_res + 1)
    dt = 2 * np.pi / Npix
    # return np.array([2*mini, 2*mini+1])*dt + dt/2
    # the fiber bundle grid on SO3 is weird
    # the next resolution level's nearest neighbors in SO3 are not
    # necessarily the nearest neighbor grid points in S1
    # include the 13 neighbors for now... eventually learn/memoize the mapping
    ind = np.arange(2 * mini - 1, 2 * mini + 3)
    if ind[0] < 0:
        ind[0] += Npix
    return ind * dt + dt / 2, ind


def get_s2_neighbor(mini, cur_res):
    """
    Return the 4 nearest neighbors on S2 at the next resolution level
    """
    Nside = 2 ** (cur_res + 1)
    ind = np.arange(4) + 4 * mini
    return pix2ang(Nside, ind, nest=True), ind


def get_base_ind(ind, base):
    """
    Return the corresponding S2 and S1 grid index for an index on the base SO3 grid
    """
    Np = 6 * 2**base
    psii = ind % Np
    thetai = ind // Np
    return np.stack((thetai, psii), axis=1)


def get_neighbor(quat, s2i, s1i, cur_res):
    """
    Return the 8 nearest neighbors on SO3 at the next resolution level
    """
    (theta, phi), s2_nexti = get_s2_neighbor(s2i, cur_res)
    psi, s1_nexti = get_s1_neighbor(s1i, cur_res)
    quat_n = hopf_to_quat(
        np.repeat(theta, len(psi)), np.repeat(phi, len(psi)), np.tile(psi, len(theta))
    )
    ind = np.array([np.repeat(s2_nexti, len(psi)), np.tile(s1_nexti, len(theta))])
    ind = ind.T
    # find the 8 nearest neighbors of 16 possible points
    # need to check distance from both +q and -q
    dists = np.minimum(
        np.sum((quat_n - quat) ** 2, axis=1), np.sum((quat_n + quat) ** 2, axis=1)
    )
    ii = np.argsort(dists)[:8]
    return quat_n[ii], ind[ii]


try:
    with open(f"{os.path.dirname(__file__)}/healpy_grid.json") as hf:
        _GRIDS = {int(k): np.array(v).T for k, v in json.load(hf).items()}
except IOError:
    print(
        "WARNING: Couldn't load cached healpy grid; will fall back to importing healpy"
    )
    _GRIDS = None


def pix2ang(Nside, ipix, nest=False, lonlat=False):
    if _GRIDS is not None and Nside in _GRIDS and nest and not lonlat:
        return _GRIDS[Nside][ipix].T
    else:
        try:
            import healpy
        except ImportError:
            raise RuntimeError(
                "You need to `pip install healpy` to run with non-standard grid sizes."
            )
        return healpy.pix2ang(Nside, ipix, nest=nest, lonlat=lonlat)
