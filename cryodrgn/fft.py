from itertools import repeat
import logging
from multiprocessing.pool import ThreadPool as Pool
import numpy as np
import torch
from torch.fft import fftshift, ifftshift, fft2, fftn, ifftn


logger = logging.getLogger(__name__)


def fft2_center(img):
    x = fftshift(img, dim=(-1, -2))
    x = fft2(x, dim=(-1, -2))
    x = fftshift(x, dim=(-1, -2))
    return x


def fftn_center(img):
    return fftshift(fftn(fftshift(img)))


def ifftn_center(img):
    return ifftshift(ifftn(ifftshift(img)))


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
        raise RuntimeError

    if inplace:
        dest = src
    else:
        if isinstance(src, np.ndarray):
            dest = src.copy()
        elif isinstance(src, torch.Tensor):
            dest = src.detach()

    with Pool(n_workers) as pool:
        pool.starmap(
            g, zip(repeat(src, n_chunks), repeat(dest, n_chunks), chunk_indices)
        )
        return dest.reshape(data.shape)


def normalize(
    img, mean=0, std=None, std_n=1_000_000, inplace=False, chunksize=None, n_workers=1
):
    def _normalize(img, mean, std):
        return (img - mean) / std

    if std is None:
        # Since std is a memory consuming process, use the first std_n samples for std determination
        std = torch.std(img[:std_n, ...])

    retVal = transform_in_chunks(
        _normalize,
        img,
        chunksize=chunksize,
        n_workers=n_workers,
        inplace=inplace,
        mean=mean,
        std=std,
    )

    logger.info(f"Normalized HT by {mean} +/- {std}")
    return retVal


def ht2_center(img, inplace=False, chunksize=None, n_workers=1):
    def _ht2_center(img):
        _img = fft2_center(img)
        return _img.real - _img.imag

    return transform_in_chunks(
        _ht2_center, img, chunksize=chunksize, n_workers=n_workers, inplace=inplace
    )


def htn_center(img):
    f = fftshift(fftn(fftshift(img)))
    return (f.real - f.imag).astype(img.dtype)


def iht2_center(img):
    img = fft2_center(img)
    img /= img.shape[-1] * img.shape[-2]
    return img.real - img.imag


def ihtn_center(V):
    V = fftshift(V)
    V = fftn(V)
    V = fftshift(V)
    V /= torch.product(V.shape)
    return (V.real - V.imag).astype(V.dtype)


def symmetrize_ht(ht, preallocated=False):
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
        sym_ht = torch.empty((B, D + 1, D + 1), dtype=ht.dtype)
        sym_ht[:, 0:-1, 0:-1] = ht
    assert D % 2 == 0
    sym_ht[:, -1, :] = sym_ht[:, 0]  # last row is the first row
    sym_ht[:, :, -1] = sym_ht[:, :, 0]  # last col is the first col
    sym_ht[:, -1, -1] = sym_ht[:, 0, 0]
    if len(sym_ht) == 1:
        sym_ht = sym_ht[0]
    return sym_ht
