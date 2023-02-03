import numpy as np

from cryodrgn.numeric.base_fft import FFT


class NumpyFFT(FFT):
    def fft(self, x, axis=-1, workers=-1):
        return np.fft.fft(x, axis=axis)

    def ifft(self, x, axis=-1, workers=-1):
        return np.fft.ifft(x, axis=axis)

    def fft2(self, x, axes=(-2, -1), workers=-1):
        return np.fft.fft2(x, axes=axes)

    def ifft2(self, x, axes=(-2, -1), workers=-1):
        return np.fft.ifft2(x, axes=axes)

    def fftn(self, x, axes=None, workers=-1):
        return np.fft.fftn(x, axes=axes)

    def ifftn(self, x, axes=None, workers=-1):
        return np.fft.ifftn(x, axes=axes)

    def fftshift(self, x, axes=None):
        return np.fft.fftshift(x, axes=axes)

    def ifftshift(self, x, axes=None):
        return np.fft.ifftshift(x, axes=axes)
