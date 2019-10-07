'''
Backproject a stack of images via linear interpolation
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
    parser.add_argument('mrcs', help='Input .mrcs image stack')
    parser.add_argument('pkl', help='Rotation matrices or Euler angles (EMAN convention)')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output .mrc file')
    parser.add_argument('--indices',help='Indices to iterate over (pkl)')
    parser.add_argument('--trans', type=os.path.abspath, help='Optionally provide translations (.pkl)')
    parser.add_argument('--tscale', type=float, help='Scale all translations by this amount')
    parser.add_argument('--first', type=int, help='Backproject the first N images')
    parser.add_argument('--tilt', help='Tilt series .mrcs image stack')
    parser.add_argument('--tilt-deg', type=float, default=45, help='Right-handed x-axis tilt offset in degrees (default: %(default)s)')
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

    if args.tilt is not None:
        images_tilt, _, _ = mrc.parse_mrc(args.tilt,lazy=True)
        assert len(images) == len(images_tilt)
        assert images[0].get().shape == images_tilt[0].get().shape
        tilt = utils.xrot(args.tilt_deg)


    rot = utils.load_pkl(args.pkl)
    if len(rot.shape) == 2:
        log('Converting poses from EMAN euler angles to rotation matrices')
        rot = np.array([utils.R_from_eman(*x) for x in rot])

    if args.trans:
        trans = utils.load_pkl(args.trans)
        if args.tscale is not None:
            log('Scaling translations by {}'.format(args.tscale))
            trans *= args.tscale
    else:
        trans = None

    n, m = images[0].get().shape
    log('Loaded {} {}x{} images'.format(N,n,m))
    assert n == m, "Image dimensions must be square"
    D = n

    V = np.zeros((D,D,D),dtype=complex)
    counts = np.zeros((D,D,D))

    xx,yy = np.meshgrid(np.arange(-D/2,D/2),np.arange(-D/2,D/2))
    zz = np.zeros(xx.shape)
    COORD = np.array([xx.ravel(), yy.ravel(), zz.ravel()])
    MASK = np.where(np.sum(COORD**2,axis=0)**.5 <=(D/2-1))
    COORD = COORD[:,MASK[0]]

    # 2D lattice spanning [-.5,.5) for FT phase shift 
    TCOORD = np.stack([xx, yy],axis=2)/D # DxDx2

    if args.indices:
        iterator = pickle.load(open(args.indices,'rb'))
    elif args.first:
        iterator = range(args.first)
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
        r = rot[ii]
        ff_coord = np.dot(r.T,COORD)
        add_slice(V,counts,ff_coord,ff,D)

        # tilt series
        if args.tilt is not None:
            ff_tilt = fft.fft2_center(images_tilt[ii].get())
            if trans is not None: ff_tilt *= tfilt
            ff_tilt = ff_tilt.ravel()[MASK]
            ff_coord = np.dot(np.dot(r.T, tilt.T), COORD)
            add_slice(V,counts,ff_coord,ff_tilt,D)

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
