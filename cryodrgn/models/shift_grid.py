import numpy as np


def grid_1d(resol: int, extent: int, ngrid: int, shift: int = 0) -> np.ndarray:
    Npix = ngrid * 2**resol
    dt = 2 * extent / Npix
    grid = np.arange(Npix, dtype=np.float32) * dt + dt / 2 - extent + shift
    return grid


def grid_2d(
    resol: int, extent: int, ngrid: int, xshift: int = 0, yshift: int = 0
) -> np.ndarray:
    x = grid_1d(resol, extent, ngrid, shift=xshift)
    y = grid_1d(resol, extent, ngrid, shift=yshift)
    # convention: x is fast dim, y is slow dim
    grid = np.stack(np.meshgrid(x, y), -1)
    return grid.reshape(-1, 2)


def base_shift_grid(
    resol: int, extent: int, ngrid: int, xshift: int = 0, yshift: int = 0
) -> np.ndarray:
    return grid_2d(resol, extent, ngrid, xshift, yshift)


# Neighbor Finding


def get_1d_neighbor(mini, cur_res, extent, ngrid):
    Npix = ngrid * 2 ** (cur_res + 1)
    dt = 2 * extent / Npix
    ind = np.array([2 * mini, 2 * mini + 1], dtype=np.float32)
    return dt * ind + dt / 2 - extent, ind


def get_base_ind(ind, ngrid):
    xi = ind % ngrid
    yi = ind // ngrid
    return np.stack((xi, yi), axis=1)


def get_neighbor(xi, yi, cur_res, extent, ngrid):
    """
    Return the 4 nearest neighbors at the next resolution level
    """
    x_next, xii = get_1d_neighbor(xi, cur_res, extent, ngrid)
    y_next, yii = get_1d_neighbor(yi, cur_res, extent, ngrid)
    t_next = np.stack(np.meshgrid(x_next, y_next), -1).reshape(-1, 2)
    ind_next = np.stack(np.meshgrid(xii, yii), -1).reshape(-1, 2)
    return t_next, ind_next
