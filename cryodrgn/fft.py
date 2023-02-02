from itertools import repeat
from multiprocessing.pool import ThreadPool as Pool
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


def transform_in_chunks(
    f, data, chunksize=None, n_workers=1, inplace=False, *args, **kwargs
):
    def g(src, dest, indices):
        dest[indices, ...] = f(src[indices, ...], *args, **kwargs)

    ndim = data.ndim
    n = data.shape[0] if ndim == 3 else 1
    indices = np.arange(n)
    chunksize = min(n, chunksize) if chunksize is not None else n
    chunk_indices = np.array_split(indices, n // chunksize)
    n_chunks = len(chunk_indices)

    if data.ndim == 2:
        src = data[np.newaxis, ...]
    elif data.ndim == 3:
        src = data
    else:
        raise NotImplementedError

    dest = src if inplace else src.copy()
    with Pool(n_workers) as pool:
        pool.starmap(
            g, zip(repeat(src, n_chunks), repeat(dest, n_chunks), chunk_indices)
        )
        return dest.reshape(data.shape)


def normalize(
    img, mean=0, std=None, std_n=1000, inplace=False, chunksize=None, n_workers=1
):
    def _normalize(img, mean, std):
        return (img - mean) / std

    if std is None:
        # Since std is a memory consuming process, use the first std_n samples for std determination
        pp = np if isinstance(img, np.ndarray) else cp
        std = pp.std(img[:std_n, ...])

    return transform_in_chunks(
        _normalize,
        img,
        chunksize=chunksize,
        n_workers=n_workers,
        inplace=inplace,
        mean=mean,
        std=std,
    )


def ht2_center(img, inplace=False, chunksize=None, n_workers=1):
    def _ht2_center(img):
        _img = fft2_center(img)
        return (_img.real - _img.imag).astype(img.dtype)

    return transform_in_chunks(
        _ht2_center, img, chunksize=chunksize, n_workers=n_workers, inplace=inplace
    )


def htn_center(img):
    pp = np if isinstance(img, np.ndarray) else cp

    f = pp.fft.fftshift(pp.fft.fftn(pp.fft.fftshift(img)))
    return (f.real - f.imag).astype(img.dtype)


def iht2_center(img):
    img = fft2_center(img)
    img /= img.shape[-1] * img.shape[-2]
    return (img.real - img.imag).astype(img.dtype)


def ihtn_center(V):
    pp = np if isinstance(V, np.ndarray) else cp

    V = pp.fft.fftshift(V)
    V = pp.fft.fftn(V)
    V = pp.fft.fftshift(V)
    V /= pp.product(V.shape)
    return (V.real - V.imag).astype(V.dtype)


def symmetrize_ht(ht, preallocated=False):
    pp = np if isinstance(ht, np.ndarray) else cp

    if preallocated:
        D = ht.shape[-1] - 1
        sym_ht = ht
        if len(sym_ht.shape) == 2:
            sym_ht = sym_ht.reshape(1, *sym_ht.shape)
        assert len(sym_ht.shape) == 3
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
