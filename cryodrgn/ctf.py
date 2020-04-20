import numpy as np
import torch
from . import utils
log = utils.log

def compute_ctf(freqs, dfu, dfv, dfang, volt, cs, w, phase_shift=0, bfactor=None):
    '''
    Compute the 2D CTF
   
    Input: 
        freqs (np.ndarray) Nx2 or BxNx2 tensor of 2D spatial frequencies
        dfu (float or Bx1 tensor): DefocusU (Angstrom)
        dfv (float or Bx1 tensor): DefocusV (Angstrom)
        dfang (float or Bx1 tensor): DefocusAngle (degrees)
        volt (float or Bx1 tensor): accelerating voltage (kV)
        cs (float or Bx1 tensor): spherical aberration (mm)
        w (float or Bx1 tensor): amplitude contrast ratio
        phase_shift (float or Bx1 tensor): degrees 
        bfactor (float or Bx1 tensor): envelope fcn B-factor (Angstrom^2)
    '''
    assert freqs.shape[-1] == 2
    # convert units
    volt = volt * 1000
    cs = cs * 10**7
    dfang = dfang * np.pi / 180
    phase_shift = phase_shift * np.pi / 180
    
    # lam = sqrt(h^2/(2*m*e*Vr)); Vr = V + (e/(2*m*c^2))*V^2
    lam = 12.2639 / (volt + 0.97845e-6 * volt**2)**.5
    x = freqs[...,0]
    y = freqs[...,1]
    ang = torch.atan2(y,x)
    s2 = x**2 + y**2
    df = .5*(dfu + dfv + (dfu-dfv)*torch.cos(2*(ang-dfang)))
    gamma = 2*np.pi*(-.5*df*lam*s2 + .25*cs*lam**3*s2**2) - phase_shift
    ctf = (1-w**2)**.5*torch.sin(gamma) - w*torch.cos(gamma) 
    if bfactor is not None:
        ctf *= torch.exp(-bfactor/4*s2)
    return ctf

def compute_ctf_np(freqs, dfu, dfv, dfang, volt, cs, w, phase_shift=0, bfactor=None):
    '''
    Compute the 2D CTF
   
    Input: 
        freqs (np.ndarray) Nx2 array of 2D spatial frequencies
        dfu (float): DefocusU (Angstrom)
        dfv (float): DefocusV (Angstrom)
        dfang (float): DefocusAngle (degrees)
        volt (float): accelerating voltage (kV)
        cs (float): spherical aberration (mm)
        w (float): amplitude contrast ratio
        phase_shift (float): degrees 
        bfactor (float): envelope fcn B-factor (Angstrom^2)
    '''
    # convert units
    volt = volt * 1000
    cs = cs * 10**7
    dfang = dfang * np.pi / 180
    phase_shift = phase_shift * np.pi / 180
    
    # lam = sqrt(h^2/(2*m*e*Vr)); Vr = V + (e/(2*m*c^2))*V^2
    lam = 12.2639 / np.sqrt(volt + 0.97845e-6 * volt**2)
    x = freqs[:,0]
    y = freqs[:,1]
    ang = np.arctan2(y,x)
    s2 = x**2 + y**2
    df = .5*(dfu + dfv + (dfu-dfv)*np.cos(2*(ang-dfang)))
    gamma = 2*np.pi*(-.5*df*lam*s2 + .25*cs*lam**3*s2**2) - phase_shift
    ctf = np.sqrt(1-w**2)*np.sin(gamma) - w*np.cos(gamma) 
    if bfactor is not None:
        ctf *= np.exp(-bfactor/4*s2)
    return np.require(ctf,dtype=freqs.dtype)

def print_ctf_params(params):
    assert len(params) == 9
    log('Image size (pix)  : {}'.format(int(params[0])))
    log('A/pix             : {}'.format(params[1]))
    log('DefocusU (A)      : {}'.format(params[2]))
    log('DefocusV (A)      : {}'.format(params[3]))
    log('Dfang (deg)       : {}'.format(params[4]))
    log('voltage (kV)      : {}'.format(params[5]))
    log('cs (mm)           : {}'.format(params[6]))
    log('w                 : {}'.format(params[7]))
    log('Phase shift (deg) : {}'.format(params[8]))

def plot_ctf(D,Apix,ctf_params):
    assert len(ctf_params) == 7
    import matplotlib.pyplot as plt
    import seaborn as sns
    freqs = np.stack(np.meshgrid(np.linspace(-.5,.5,D,endpoint=False),np.linspace(-.5,.5,D,endpoint=False)),-1)/Apix
    freqs = freqs.reshape(-1,2)
    c = compute_ctf_np(freqs, *ctf_params)
    sns.heatmap(c.reshape(D, D))

def load_ctf_for_training(D, ctf_params_pkl):
    assert D%2 == 0
    ctf_params = utils.load_pkl(ctf_params_pkl)
    assert ctf_params.shape[1] == 9
    # Replace original image size with current dimensions
    Apix = ctf_params[0,0]*ctf_params[0,1]/D
    ctf_params[:,0] = D
    ctf_params[:,1] = Apix
    print_ctf_params(ctf_params[0])
    # Slice out the first column (D)
    return ctf_params[:,1:]
    


