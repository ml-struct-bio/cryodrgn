from itertools import repeat
import logging
from multiprocessing.pool import ThreadPool as Pool
import numpy as np
import torch
from torch.fft import fftshift, ifftshift, fft2, fftn, ifftn


logger = logging.getLogger(__name__)


def transform_in_chunks(
    f, data, chunksize=None, n_workers=1, inplace=False, *args, **kwargs
):
    def g(src, dest, indices):
        logger.debug(f"Transforming chunk of size {len(indices)}")
        dest[indices, ...] = f(src[indices, ...], *args, **kwargs)

    ndim = data.ndim
    if data.ndim == 2:
        n = 1
        src = data[np.newaxis, ...]
    elif data.ndim == 3:
        n = data.shape[0]
        src = data
    else:
        raise RuntimeError(f"Unsupported data ndim {ndim}")

    indices = np.arange(n)
    chunksize = min(n, chunksize) if chunksize is not None else n
    chunk_indices = np.array_split(indices, n // chunksize)
    n_chunks = len(chunk_indices)

    if inplace:
        dest = src
    else:
        if isinstance(src, np.ndarray):
            dest = src.copy()
        elif isinstance(src, torch.Tensor):
            dest = src.detach().clone()

    with Pool(n_workers) as pool:
        pool.starmap(
            g, zip(repeat(src, n_chunks), repeat(dest, n_chunks), chunk_indices)
        )
        return dest.reshape(data.shape)


def normalize(
    img, mean=0, std=None, std_n=None, inplace=False, chunksize=None, n_workers=1
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

    logger.info(f"Normalized by {mean} +/- {std}")
    return retVal


def fft2_center(img):
    return fftshift(fft2(fftshift(img, dim=(-1, -2))), dim=(-1, -2))


def fftn_center(img):
    return fftshift(fftn(fftshift(img)))


def ifftn_center(img):
    if isinstance(img, np.ndarray):
        # Note: We can't just typecast a complex ndarray using torch.Tensor(array) !
        img = torch.complex(torch.Tensor(img.real), torch.Tensor(img.imag))
    x = ifftshift(img)
    y = ifftn(x)
    z = ifftshift(y)
    return z


def ht2_center(img, inplace=False, chunksize=None, n_workers=1):
    def _ht2_center(img):
        _img = fft2_center(img)
        return _img.real - _img.imag

    return transform_in_chunks(
        _ht2_center, img, chunksize=chunksize, n_workers=n_workers, inplace=inplace
    )


def htn_center(img):
    _img = fftshift(fftn(fftshift(img)))
    return _img.real - _img.imag


def iht2_center(img):
    img = fft2_center(img)
    img /= img.shape[-1] * img.shape[-2]
    return img.real - img.imag


def ihtn_center(img):
    img = fftshift(img)
    img = fftn(img)
    img = fftshift(img)
    img /= torch.prod(torch.tensor(img.shape, device=img.device))
    return img.real - img.imag


def symmetrize_ht(ht, preallocated=False):
    if ht.ndim == 2:
        ht = ht[np.newaxis, ...]
    assert ht.ndim == 3
    n = ht.shape[0]

    if preallocated:
        D = ht.shape[-1] - 1
        sym_ht = ht
    else:
        D = ht.shape[-1]
        sym_ht = torch.empty((n, D + 1, D + 1), dtype=ht.dtype)
        sym_ht[:, 0:-1, 0:-1] = ht

    assert D % 2 == 0
    sym_ht[:, -1, :] = sym_ht[:, 0, :]  # last row is the first row
    sym_ht[:, :, -1] = sym_ht[:, :, 0]  # last col is the first col
    sym_ht[:, -1, -1] = sym_ht[:, 0, 0]  # last corner is first corner

    if n == 1:
        sym_ht = sym_ht[0, ...]

    return sym_ht
