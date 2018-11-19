import numpy as np

def fft2_center(img):
    return np.fft.fftshift(np.fft.fft2(np.fft.fftshift(img)))

def ifftn(V):
    V = np.fft.ifftshift(V)
    V = np.fft.ifftn(V)
    V = np.fft.ifftshift(V)
    return V
