import numpy as np
import numpy.fft
import torch
import torch.fft
import pytest


img_np = np.random.random((100, 32, 32))
img_torch = torch.tensor(img_np)


@pytest.mark.fast
@pytest.mark.unit
def test_fft2():
    assert np.allclose(torch.fft.fft2(img_torch).cpu().numpy(), numpy.fft.fft2(img_np))


@pytest.mark.fast
@pytest.mark.unit
def test_fftn():
    assert np.allclose(torch.fft.fftn(img_torch).cpu().numpy(), numpy.fft.fftn(img_np))


@pytest.mark.fast
@pytest.mark.unit
def test_ifftn():
    assert np.allclose(
        torch.fft.ifftn(img_torch).cpu().numpy(), numpy.fft.ifftn(img_np)
    )


@pytest.mark.fast
@pytest.mark.unit
def test_fftshift():
    assert np.allclose(
        torch.fft.fftshift(img_torch, dim=(-1, -2)).cpu().numpy(),
        numpy.fft.fftshift(img_np, axes=(-1, -2)),
    )


@pytest.mark.fast
@pytest.mark.unit
def test_ifftshift():
    assert np.allclose(
        torch.fft.ifftshift(img_torch, dim=(-1, -2)).cpu().numpy(),
        numpy.fft.ifftshift(img_np, axes=(-1, -2)),
    )
