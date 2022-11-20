import numpy as np

try:
    import cupy as cp  # type: ignore
except ImportError:
    cp = np


def fft2_center(img):
    pp = np if isinstance(img, np.ndarray) else cp

    return pp.fft.fftshift(
        pp.fft.fft2(pp.fft.fftshift(img, axes=(-1, -2))), axes=(-1, -2)
    )


def fftn_center(img):
    pp = np if isinstance(img, np.ndarray) else cp

    return pp.fft.fftshift(pp.fft.fftn(pp.fft.fftshift(img)))


def ifftn_center(V):
    pp = np if isinstance(V, np.ndarray) else cp

    V = pp.fft.ifftshift(V)
    V = pp.fft.ifftn(V)
    V = pp.fft.ifftshift(V)
    return V


def ht2_center(img):
    f = fft2_center(img)
    return f.real - f.imag


def htn_center(img):
    pp = np if isinstance(img, np.ndarray) else cp

    f = pp.fft.fftshift(pp.fft.fftn(pp.fft.fftshift(img)))
    return f.real - f.imag


def iht2_center(img):
    img = fft2_center(img)
    img /= img.shape[-1] * img.shape[-2]
    return img.real - img.imag


def ihtn_center(V):
    pp = np if isinstance(V, np.ndarray) else cp

    V = pp.fft.fftshift(V)
    V = pp.fft.fftn(V)
    V = pp.fft.fftshift(V)
    V /= pp.product(V.shape)
    return V.real - V.imag


def symmetrize_ht(ht, preallocated=False):
    pp = np if isinstance(ht, np.ndarray) else cp

    if preallocated:
        D = ht.shape[-1] - 1
        sym_ht = ht
    else:
        if len(ht.shape) == 2:
            ht = ht.reshape(1, *ht.shape)
        assert len(ht.shape) == 3
        D = ht.shape[-1]
        B = ht.shape[0]
        sym_ht = pp.empty((B, D + 1, D + 1), dtype=ht.dtype)
        sym_ht[:, 0:-1, 0:-1] = ht
    assert D % 2 == 0
    sym_ht[:, -1, :] = sym_ht[:, 0]  # last row is the first row
    sym_ht[:, :, -1] = sym_ht[:, :, 0]  # last col is the first col
    sym_ht[:, -1, -1] = sym_ht[:, 0, 0]
    if len(sym_ht) == 1:
        sym_ht = sym_ht[0]
    return sym_ht
