import numpy as np

def fft2_center(img):
    return np.fft.fftshift(np.fft.fft2(np.fft.fftshift(img,axes=(-1,-2))),axes=(-1,-2))

def fftn_center(img):
    return np.fft.fftshift(np.fft.fftn(np.fft.fftshift(img)))

def ifftn_center(V):
    V = np.fft.ifftshift(V)
    V = np.fft.ifftn(V)
    V = np.fft.ifftshift(V)
    return V

def ht2_center(img):
    f = fft2_center(img)
    return f.real-f.imag

def htn_center(img):
    f = np.fft.fftshift(np.fft.fftn(np.fft.fftshift(img)))
    return f.real-f.imag

def ihtn_center(V):
    V = np.fft.fftshift(V)
    V = np.fft.fftn(V)
    V = np.fft.fftshift(V)
    V /= np.product(V.shape)
    return V.real - V.imag

def symmetrize_ht(ht):
    if len(ht.shape) == 2:
        D = ht.shape[0]
        ht = ht.reshape(1,*ht.shape)
    else:
        assert len(ht.shape) == 3
        D = ht.shape[1]
    assert D % 2 == 0
    B = ht.shape[0]
    sym_ht = np.empty((B,D+1,D+1),dtype=ht.dtype)
    sym_ht[:,0:-1,0:-1] = ht
    sym_ht[:,-1,:] = sym_ht[:,0] # last row is the first row
    sym_ht[:,:,-1] = sym_ht[:,:,0] # last col is the first col
    sym_ht[:,-1,-1] = sym_ht[:,0,0]
    if len(sym_ht) == 1:
        sym_ht = sym_ht[0]
    return sym_ht
