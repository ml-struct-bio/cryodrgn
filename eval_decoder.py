'''
Evaluate the decoder from bnb_het.py at specified values of z
'''
import numpy as np
import sys, os
import argparse
import pickle
from datetime import datetime as dt
import matplotlib.pyplot as plt 

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable
from torch.distributions import Normal

sys.path.insert(0,os.path.abspath(os.path.dirname(__file__))+'/lib-python')
import mrc
import utils
import fft
import lie_tools
from lattice import Lattice
from models import HetOnlyVAE

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('weights', help='Model weights')
    parser.add_argument('--norm', nargs=2, type=float, required=True)
    parser.add_argument('--dim', nargs=3, default=[64,64,64], type=int)
    parser.add_argument('-z', type=np.float32, nargs='*', help='')
    parser.add_argument('--z-start', type=np.float32, nargs='*', help='')
    parser.add_argument('--z-end', type=np.float32, nargs='*', help='')
    parser.add_argument('--zfile')
    parser.add_argument('-n', type=int, default=10, help='')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output MRC or directory')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')
    parser.add_argument('--l-extent', type=float, default=1.0, help='Coordinate lattice size (default: %(default)s)')

    group = parser.add_argument_group('Architecture parameters')
    group.add_argument('--qlayers', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--qdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--zdim', type=int, default=1, help='Dimension of latent variable')
    group.add_argument('--encode-mode', default='resid', choices=('conv','resid','mlp','tilt'), help='Type of encoder network')
    group.add_argument('--players', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--pdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--enc-mask', type=int, help='Circulask mask of image for encoder')
    group.add_argument('--enc-type', choices=('geom_ft','geom_full','geom_lowf','geom_nohighf','linear_lowf','none'), default='linear_lowf', help='Type of positional encoding')

    return parser

def main(args):
    log(args)
    t1 = dt.now()

    ## set the device
    use_cuda = torch.cuda.is_available()
    log('Use cuda {}'.format(use_cuda))
    if use_cuda:
        torch.set_default_tensor_type(torch.cuda.FloatTensor)

    nz, ny, nx = args.dim
    assert nz == ny == nx
    D = nz+1
    lattice = Lattice(D, extent=args.l_extent)
    if args.enc_mask: args.enc_mask = lattice.get_circular_mask(args.enc_mask)
    model = HetOnlyVAE(lattice, args.qlayers, args.qdim, args.players, args.pdim,
                args.zdim, encode_mode=args.encode_mode, enc_mask=args.enc_mask, enc_type=args.enc_type)

    log('Loading weights from {}'.format(args.weights))
    checkpoint = torch.load(args.weights)
    model.load_state_dict(checkpoint['model_state_dict'])

    model.eval()

    if args.z_start or args.zfile:
        if args.z_start:
            assert args.z_end
            assert not args.z
            assert not args.zfile
            args.z_start = np.array(args.z_start)
            args.z_end = np.array(args.z_end)
            z = np.repeat(np.arange(args.n,dtype=np.float32), args.zdim).reshape((args.n,args.zdim))
            z *= ((args.z_end - args.z_start)/(args.n-1))
            z += args.z_start
        else:
            assert not args.z_start
            z = np.loadtxt(args.zfile).reshape(-1,args.zdim)

        if not os.path.exists(args.o):
            os.makedirs(args.o)

        for i,zz in enumerate(z):
            log(zz)
            vol = model.decoder.eval_volume(lattice.coords, lattice.D, lattice.extent, args.norm, zz) 
            out_mrc = '{}/traj{}.mrc'.format(args.o,i)
            mrc.write(out_mrc, vol.astype(np.float32))

    else:
        z = np.array(args.z)
        log(z)
        vol = model.decoder.eval_volume(lattice.coords, lattice.D, lattice.extent, args.norm, z) 
        mrc.write(args.o, vol.astype(np.float32))

    td = dt.now()-t1
    log('Finsihed in {}'.format(td))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)

