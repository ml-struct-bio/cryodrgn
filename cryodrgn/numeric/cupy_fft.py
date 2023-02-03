import cupy as cp

from cryodrgn.numeric.base_fft import FFT


class CupyFFT(FFT):
    """
    Define a unified wrapper class for Cupy FFT functions

    To be consistent with Scipy and Pyfftw, not all arguments are included.
    """

    def fft(self, x, axis=-1, workers=-1):
        return cp.fft.fft(x, axis=axis)

    def ifft(self, x, axis=-1, workers=-1):
        return cp.fft.ifft(x, axis=axis)

    def fft2(self, x, axes=(-2, -1), workers=-1):
        return cp.fft.fft2(x, axes=axes)

    def ifft2(self, x, axes=(-2, -1), workers=-1):
        return cp.fft.ifft2(x, axes=axes)

    def fftn(self, x, axes=None, workers=-1):
        return cp.fft.fftn(x, axes=axes)

    def ifftn(self, x, axes=None, workers=-1):
        return cp.fft.ifftn(x, axes=axes)

    def fftshift(self, x, axes=None):
        return cp.fft.fftshift(x, axes=axes)

    def ifftshift(self, x, axes=None):
        return cp.fft.ifftshift(x, axes=axes)
