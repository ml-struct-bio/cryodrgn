'''
Backproject a stack of images via linear interpolation
'''

import argparse
import numpy as np
import sys, os
import time
import pickle

import torch

sys.path.insert(0,'{}/lib-python'.format(os.path.dirname(os.path.abspath(__file__))))
import utils
import mrc
import fft
from pose import PoseTracker
from lattice import Lattice
import dataset
import ctf

log = utils.log

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mrcs', help='Input .mrcs image stack')
    parser.add_argument('--poses', type=os.path.abspath, nargs='*', required=True, help='Image rotations and optionally translations (.pkl)')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output .mrc file')
    parser.add_argument('--invert-data', action='store_true')
    parser.add_argument('--indices',help='Indices to iterate over (pkl)')
    parser.add_argument('--trans', type=os.path.abspath, help='Optionally provide translations (.pkl)')
    parser.add_argument('--tscale', type=float, help='Scale all translations by this amount')
    parser.add_argument('--first', type=int, help='Backproject the first N images')
    parser.add_argument('--tilt', help='Tilt series .mrcs image stack')
    parser.add_argument('--tilt-deg', type=float, default=45, help='Right-handed x-axis tilt offset in degrees (default: %(default)s)')
    parser.add_argument('--ctf', metavar='pkl', type=os.path.abspath, help='CTF parameters (.pkl) if particle stack is not phase flipped')
    return parser

def add_slice(V, counts, ff_coord, ff, D):
    d2 = int(D/2)
    ff_coord = ff_coord.transpose(0,1)
    xf, yf, zf = ff_coord.floor().long()
    xc, yc, zc = ff_coord.ceil().long()
    def add_for_corner(xi,yi,zi):
        dist = torch.stack([xi,yi,zi]).float() - ff_coord
        w = 1 - dist.pow(2).sum(0).pow(.5)
        w[w<0]=0
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
    t1 = time.time()    
    log(args)
    if not os.path.exists(os.path.dirname(args.o)):
        os.makedirs(os.path.dirname(args.o))

    ## set the device
    use_cuda = torch.cuda.is_available()
    device = torch.device('cuda' if use_cuda else 'cpu')
    log('Use cuda {}'.format(use_cuda))
    if use_cuda:
        torch.set_default_tensor_type(torch.cuda.FloatTensor)

    # load the particles
    if args.tilt is None:
        data = dataset.LazyMRCData(args.mrcs, norm=(0,1), invert_data=args.invert_data)
        tilt = None
    else:
        data = dataset.TiltMRCData(args.mrcs, args.tilt, norm=(0,1), invert_data=args.invert_data)
        tilt = torch.tensor(utils.xrot(args.tilt_deg).astype(np.float32))
    D = data.D
    Nimg = data.N

    lattice = Lattice(D, extent=D//2)

    posetracker = PoseTracker.load(args.poses, Nimg, None, args.tscale, None)

    if args.ctf is not None:
        log('Loading ctf params from {}'.format(args.ctf))
        ctf_params = utils.load_pkl(args.ctf)
        assert ctf_params.shape == (Nimg, 7)
        ctf.print_ctf_params(ctf_params[0])
        ctf_params = torch.tensor(ctf_params)
    else: ctf_params = None

    V = torch.zeros((D,D,D))
    counts = torch.zeros((D,D,D))
    
    mask = lattice.get_circular_mask(D//2)

    if args.indices:
        iterator = pickle.load(open(args.indices,'rb'))
    elif args.first:
        iterator = range(args.first)
    else:
        iterator = range(Nimg)

    for ii in iterator:
        if ii%100==0: log('image {}'.format(ii))
        r, t = posetracker.get_pose(ii)
        ff = data.get(ii)
        if tilt is not None:
            ff, ff_tilt = ff # EW
        ff = torch.tensor(ff)
        ff = ff.view(-1)[mask]
        if ctf_params is not None:
            freqs = lattice.freqs2d/ctf_params[ii,0]
            c = ctf.compute_ctf(freqs, *ctf_params[ii,1:]).view(-1)[mask]
            ff *= c.sign()
        if t is not None:
            ff = lattice.translate_ht(ff.view(1,-1),t.view(1,1,2), mask)
        ff_coord = lattice.coords[mask] @ r
        add_slice(V, counts, ff_coord, ff, D)

        # tilt series
        if args.tilt is not None:
            ff_tilt = torch.tensor(ff_tilt)
            ff_tilt = ff_tilt.view(-1)[mask]
            if ff_tilt is not None:
                ff_tilt *= c.sign()
            if t is not None:
                ff_tilt = lattice.translate_ht(ff_tilt.view(1,-1), t.view(1,1,2), mask)
            ff_coord = lattice.coords[mask] @ tilt @ r
            add_slice(V, counts, ff_coord, ff_tilt, D)

    z = np.where(counts == 0.0)
    td = time.time()-t1
    log('Backprojected {} images in {}s ({}s per image)'.format(Nimg, td, td/Nimg ))
    log('{}% voxels missing data'.format(100*len(z[0])/D**3))
    counts[z] = 1.0
    V /= counts
    V = fft.ihtn_center(V[0:-1,0:-1,0:-1].cpu().numpy())
    mrc.write(args.o,V.astype('float32'))

if __name__ == '__main__':
    main(parse_args().parse_args())
