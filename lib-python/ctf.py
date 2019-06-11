import numpy as np
import utils
log = utils.log

def compute_ctf(freqs, dfu, dfv, dfang, volt, cs, w, phase_shift=0, bfactor=None):
    '''
    Compute the 2d CTF
   
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
    gamma = 2*np.pi*(-.5*df*lam*s2 + .25*cs*lam**3*s2**2) + phase_shift
    ctf = np.sqrt(1-w**2)*np.sin(gamma) - w*np.cos(gamma) 
    if bfactor is not None:
        ctf *= np.exp(-bfactor/4*s2)
    return np.require(ctf,dtype=freqs.dtype)

def print_ctf_params(params):
    assert len(params) == 7
    log('A/pix       : {}'.format(params[0]))
    log('DefocusU (A): {}'.format(params[1]))
    log('DefocusV (A): {}'.format(params[2]))
    log('Dfang (deg) : {}'.format(params[3]))
    log('voltage (kV): {}'.format(params[4]))
    log('cs (mm)     : {}'.format(params[5]))
    log('w           : {}'.format(params[6]))
