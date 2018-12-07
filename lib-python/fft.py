import numpy as np

def fft2_center(img):
    return np.fft.fftshift(np.fft.fft2(np.fft.fftshift(img)))

def ifftn(V):
    V = np.fft.ifftshift(V)
    V = np.fft.ifftn(V)
    V = np.fft.ifftshift(V)
    return V

def ht2_center(img):
    f = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(img)))
    return f.real-f.imag

def ihtn_center(V):
    V = np.fft.fftshift(V)
    V = np.fft.fftn(V)
    V = np.fft.fftshift(V)
    V /= np.product(V.shape)
    return V.real - V.imag
