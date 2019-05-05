'''
Reconstruct 3D density from 2D images with assigned angles

Ellen Zhong
8/8/18
'''

import argparse
import numpy as np
import sys, os
import time
import pickle

sys.path.insert(0,'{}/lib-python'.format(os.path.dirname(os.path.abspath(__file__))))
import utils
import mrc
import fft

log = utils.log

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mrcs', help='Input MRCs stack')
    parser.add_argument('pkl', help='EMAN euler angles')
    parser.add_argument('-o', type=os.path.abspath,help='Output MRC file')
    parser.add_argument('--is-rot',action='store_true',help='Input angles are rotation matrices')
    parser.add_argument('--indices',help='Indices to iterator over (pkl)')
    parser.add_argument('--trans', type=os.path.abspath, help='Optionally provide translations (.pkl)')
    parser.add_argument('--tscale',type=float,help='Scale all translations by this amount')
    return parser

def add_slice(V, counts, ff_coord, ff, D):
    d2 = int(D/2)
    xf, yf, zf = ff_coord.astype(int) # floor
    xc, yc, zc = np.ceil(ff_coord).astype(int) # ceiling
    def add_for_corner(xi,yi,zi):
        #f = (xi < 128) * (xi >= -128) * (yi < 128) * (yi >= -128) * (zi < 128) * (zi >= -128)
        dist = np.array([xi,yi,zi]) - ff_coord
        w = 1 - np.sum(dist**2, axis=0)**.5
        w[w<0]=0
        #w = np.exp(-np.sum(dist**2,axis=0)/2/.05)/(np.pi*2*.05)**1.5 # pdf of 3d gaussian with covariance diag(0.05)
        V[(zi+d2,yi+d2,xi+d2)] += w*ff
        counts[(zi+d2,yi+d2,xi+d2)] += w
    add_for_corner(xf,yf,zf)
    add_for_corner(xc,yf,zf)
    add_for_corner(xf,yc,zf)
    add_for_corner(xf,yf,zc)
    add_for_corner(xc,yc,zf)
    add_for_corner(xf,yc,zc)
    add_for_corner(xc,yf,zc)
    add_for_corner(xc,yc,zc)
    return V, counts

def main(args):
    if not os.path.exists(os.path.dirname(args.o)):
        os.makedirs(os.path.dirname(args.o))

    t1 = time.time()    
    images, _ , _ = mrc.parse_mrc(args.mrcs,lazy=True)
    N = len(images)
    angles = utils.load_angles(args.pkl)
    if len(angles) < N:
        log('Warning: # images != # angles. Backprojecting first {} images'.format(len(angles)))
        N = len(angles)
    if args.trans:
        trans = utils.load_angles(args.trans)
        if args.tscale:
            trans *= args.tscale
    else:
        trans = None

    n, m = images[0].get().shape
    assert n == m, "Image dimensions must be square"
    D = n

    V = np.zeros((D,D,D),dtype=complex)
    counts = np.zeros((D,D,D))

    xx,yy = np.meshgrid(np.arange(-D/2,D/2),np.arange(-D/2,D/2))
    zz = np.zeros(xx.shape)
    COORD = np.array([xx.ravel(), yy.ravel(), zz.ravel()])
    MASK = np.where(np.sum(COORD**2,axis=0)**.5 <=(D/2-1))
    COORD = COORD[:,MASK[0]]

    # we need a 2D lattice (normalized differently) for implementing the fourier shift 
    TCOORD = np.stack([xx, yy],axis=2)/D # DxDx2

    if args.indices:
        iterator = pickle.load(open(args.indices,'rb'))
    else:
        iterator = range(N)
    for ii in iterator:
        if ii%100==0: log('image {}'.format(ii))
        ff = fft.fft2_center(images[ii].get())
        if trans is not None:
            tfilt = np.dot(TCOORD,trans[ii])*-2*np.pi
            tfilt = np.cos(tfilt) + np.sin(tfilt)*1j
            ff *= tfilt
        ff = ff.ravel()[MASK]
        if args.is_rot:
            rot = angles[ii]
        else:
            rot = utils.R_from_eman(angles[ii,0],angles[ii,1],angles[ii,2])
        ff_coord = np.dot(rot.T,COORD)
        add_slice(V,counts,ff_coord,ff,D)
    z = np.where(counts == 0.0)
    td = time.time()-t1
    log('Backprojected {} images in {}s ({}s per image)'.format(N, td, td/N ))
    log('{}% voxels missing data'.format(100*len(z[0])/D**3))
    counts[z] = 1.0
    V /= counts
    V = fft.ifftn_center(V)
    mrc.write(args.o,V.astype('float32'))


if __name__ == '__main__':
    main(parse_args().parse_args())
