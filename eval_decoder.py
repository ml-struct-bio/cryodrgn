'''
Evaluate VAE decoder on an image stack
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
from models import HetVAE

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('weights', help='Model weights')
    parser.add_argument('--norm', nargs=2, type=float)
    parser.add_argument('--dim', nargs=3, default=[64,64,64], type=int)
    parser.add_argument('-z', type=np.float32, nargs='*', help='')
    parser.add_argument('-z-start', type=np.float32, nargs='*', help='')
    parser.add_argument('-z-end', type=np.float32, nargs='*', help='')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output MRC')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')

    group = parser.add_argument_group('Architecture parameters')
    group.add_argument('-b','--batch-size', type=int, default=100, help='Minibatch size (default: %(default)s)')
    group.add_argument('--qlayers', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--qdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--zdim', type=int, default=1, help='Dimension of latent variable')
    group.add_argument('--encode-mode', default='resid', choices=('conv','resid','mlp'), help='Type of encoder network')
    group.add_argument('--players', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--pdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    return parser

def eval_volume(model, nz, ny, nx, zval, rnorm):
    '''Evaluate the model on a nz x ny x nx lattice'''
    zdim = len(zval)
    z = torch.zeros(nx*ny,zdim,device=model.lattice.device, dtype=model.lattice.dtype)
    zval = torch.Tensor(zval).to(model.lattice.device)
    z = z + zval

    vol_f = np.zeros((nz,ny,nx),dtype=complex)
    assert not model.training
    # evaluate the volume by zslice to avoid memory overflows
    for i, dz in enumerate(np.linspace(-1,1,nz,endpoint=False)):
        x = model.lattice + torch.tensor([0,0,dz], device=model.lattice.device, 
                                                   dtype=model.lattice.dtype)
        x = torch.cat((x,z),dim=-1)
        with torch.no_grad():
            y = model.decoder(x)
            y = y.view(ny, nx).cpu().numpy()
        vol_f[i] = y
    vol = fft.ihtn_center(vol_f*rnorm[1]+rnorm[0])
    return vol, vol_f

def main(args):
    log(args)
    t1 = dt.now()

    ## set the device
    use_cuda = torch.cuda.is_available()
    log('Use cuda {}'.format(use_cuda))
    if use_cuda:
        torch.set_default_tensor_type(torch.cuda.FloatTensor)

    nz, ny, nx = args.dim
    model = HetVAE(nx, ny, nx*ny, args.qlayers, args.qdim, args.players, args.pdim,
                args.zdim, encode_mode=args.encode_mode)

    log('Loading weights from {}'.format(args.weights))
    checkpoint = torch.load(args.weights)
    model.load_state_dict(checkpoint['model_state_dict'])

    model.eval()

    if args.z_start:
        assert args.z_end
        assert not args.z
        args.z_start = np.array(args.z_start)
        args.z_end = np.array(args.z_end)
        z = np.repeat(np.arange(10,dtype=np.float32), args.zdim).reshape((10,args.zdim))
        z *= ((args.z_end - args.z_start)/10.)
        z += args.z_start
        if not os.path.exists(args.o):
            os.makedirs(args.o)

        for i,zz in enumerate(z):
            log(zz)
            vol, _ = eval_volume(model, nz, ny, nx, zz, args.norm) 
            out_mrc = '{}/traj{}.mrc'.format(args.o,i)
            mrc.write(out_mrc, vol.astype(np.float32))

    else:
        z = np.array(args.z)
        log(z)
        vol, _ = eval_volume(model, nz, ny, nx, z, args.norm) 
        mrc.write(args.o, vol.astype(np.float32))

    td = dt.now()-t1
    log('Finsihed in {}'.format(td))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)

