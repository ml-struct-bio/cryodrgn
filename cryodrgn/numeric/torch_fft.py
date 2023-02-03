import torch
from cryodrgn.numeric.base_fft import FFT


class TorchFFT(FFT):
    """
    Define a unified wrapper class for Scipy FFT functions

    To be consistent with Pyfftw, not all arguments are included.
    """

    def fft(self, x, axis=-1, workers=-1):
        return torch.fft.fft(x, dim=axis)

    def ifft(self, x, axis=-1, workers=-1):
        return torch.fft.ifft(x, dim=axis)

    def fft2(self, x, axes=(-2, -1), workers=-1):
        return torch.fft.fft2(x, dim=axes)

    def ifft2(self, x, axes=(-2, -1), workers=-1):
        return torch.fft.ifft2(x, dim=axes)

    def fftn(self, x, axes=None, workers=-1):
        return torch.fft.fftn(x, dim=axes)

    def ifftn(self, x, axes=None, workers=-1):
        return torch.fft.ifftn(x, dim=axes)

    def fftshift(self, x, axes=None):
        return torch.fft.fftshift(x, dim=axes)

    def ifftshift(self, x, axes=None):
        return torch.fft.ifftshift(x, dim=axes)
