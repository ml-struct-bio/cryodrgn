'''
Implementation of Yershova et al. "Generating uniform incremental 
grids on SO(3) using the Hopf fribration"
'''

import numpy as np
import healpy as hp
import lie_tools

#@memoize
def grid_s1(resol):
	Npix = 6*2**resol
	dt = 2*np.pi/Npix
	grid = np.arange(Npix)*dt + dt/2
	return grid

#@memoize
def grid_s2(resol):
	Nside = 2**resol
	Npix = 12*Nside*Nside
	theta, phi = hp.pix2ang(Nside, np.arange(Npix), nest=True)
	return theta, phi

def hopf_to_quat(theta, phi, psi):
	'''
	Hopf coordinates to quaternions
	theta: [0,pi)
	phi: [0, 2pi)
	psi: [0, 2pi)
	'''
	ct = np.cos(theta/2)
	st = np.sin(theta/2)
	quat = np.empty((len(theta),4))
	quat[:,0] = ct*np.cos(psi/2)
	quat[:,1] = ct*np.sin(psi/2)
	quat[:,2] = st*np.cos(phi+psi/2)
	quat[:,3] = st*np.sin(phi+psi/2)
	return quat

def SO3_grid(resol):
	theta, phi = grid_s2(resol)
	psi = grid_s1(resol)
	quat = hopf_to_quat(np.repeat(theta,len(psi)),
						np.repeat(phi,len(psi)),
						np.tile(psi,len(theta)))
	return quat #hmm convert to rot matrix?

def base_SO3_grid():
	return SO3_grid(0)

#### Neighbor finding ####

def get_s1_neighbor(mini, curr_res):
	Npix = 6*2**(curr_res+1)
	dt = 2*np.pi/Npix
	return np.array([2*mini, 2*mini+1])*dt + dt/2

def get_s2_neighbor(mini, curr_res):
	Nside = 2**(curr_res+1)
	return hp.pix2ang(Nside, np.arange(4)+4*mini, nest=True)

def find_ind_base(mini):
	psii = mini%6
	thetai = int(mini/6)
	return thetai, psii

def get_neighbor(s2i, s1i, curr_res):
	theta, phi = get_s2_neighbor(s2i, curr_res)
	s2_nexti = np.arange(4)+4*s2i
	psi = get_s1_neighbor(s1i, curr_res)
	s1_nexti = np.arange(2)+2*s1i
	quat = hopf_to_quat(np.repeat(theta,len(psi)),
						np.repeat(phi,len(psi)),
						np.tile(psi,len(theta)))
	ind = np.array([np.repeat(s2_nexti,len(psi)),
					np.tile(s1_nexti, len(theta))])
	ind = ind.T
	return quat, ind





	
