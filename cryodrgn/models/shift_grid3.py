import numpy as np


def grid_1d(resol, extent, ngrid):
    Npix = ngrid * 2**resol
    dt = 2 * extent / Npix
    grid = np.arange(Npix, dtype=np.float32) * dt + dt / 2 - extent
    return grid


def grid_3d(resol: int, extent: int, ngrid: int) -> np.ndarray:
    w = grid_1d(resol, extent, ngrid)
    # convention: x is fastest dim, z is slowest dim
    z, y, x = np.meshgrid(w, w, w, indexing="ij")
    grid = np.stack([x, y, z], -1)
    return grid.reshape(-1, 3)


def base_shift_grid(extent: int, ngrid: int) -> np.ndarray:
    return grid_3d(0, extent, ngrid)


# Neighbor Finding


def get_1d_neighbor(mini, curr_res, extent, ngrid):
    Npix = ngrid * 2 ** (curr_res + 1)
    dt = 2 * extent / Npix
    ind = np.array([2 * mini, 2 * mini + 1], dtype=np.float32)
    return dt * ind + dt / 2 - extent, ind


def get_base_id(id_, ngrid):
    xi = id_ % ngrid
    yi = id_ % ngrid**2 // ngrid
    zi = id_ // ngrid**2
    return xi, yi, zi


def get_neighbor(xi, yi, zi, curr_res, extent, ngrid):
    """
    Return the 4 nearest neighbors at the next resolution level
    """
    x_next, xii = get_1d_neighbor(xi, curr_res, extent, ngrid)
    y_next, yii = get_1d_neighbor(yi, curr_res, extent, ngrid)
    z_next, zii = get_1d_neighbor(zi, curr_res, extent, ngrid)
    z, y, x = np.meshgrid(z_next, y_next, x_next, indexing="ij")
    t_next = np.stack((x, y, z), -1).reshape(-1, 3)
    zi, yi, xi = np.meshgrid(zii, yii, xii, indexing="ij")
    id_next = np.stack((xi, yi, zi), -1).reshape(-1, 3)
    return t_next, id_next.astype(int)
