import logging
import numpy as np
import torch
from torch.fft import fftshift, ifftshift, fft2, fftn, ifftn


logger = logging.getLogger(__name__)


def normalize(img, mean=0, std=None, std_n=None):
    if std is None:
        # Since std is a memory consuming process, use the first std_n samples for std determination
        std = torch.std(img[:std_n, ...])

    logger.info(f"Normalized by {mean} +/- {std}")
    return (img - mean) / std


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


def ht2_center(img):
    _img = fft2_center(img)
    return _img.real - _img.imag


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


def symmetrize_ht(ht):
    if ht.ndim == 2:
        ht = ht[np.newaxis, ...]
    assert ht.ndim == 3
    n = ht.shape[0]

    D = ht.shape[-1]
    sym_ht = torch.empty((n, D + 1, D + 1), dtype=ht.dtype, device=ht.device)
    sym_ht[:, 0:-1, 0:-1] = ht

    assert D % 2 == 0
    sym_ht[:, -1, :] = sym_ht[:, 0, :]  # last row is the first row
    sym_ht[:, :, -1] = sym_ht[:, :, 0]  # last col is the first col
    sym_ht[:, -1, -1] = sym_ht[:, 0, 0]  # last corner is first corner

    if n == 1:
        sym_ht = sym_ht[0, ...]

    return sym_ht
