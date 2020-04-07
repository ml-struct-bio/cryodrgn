'''
Evaluate the decoder at specified values of z
'''
import numpy as np
import sys, os
import argparse
import pickle
from datetime import datetime as dt
import matplotlib.pyplot as plt 

import torch

sys.path.insert(0,os.path.abspath(os.path.dirname(__file__))+'/lib-python')
import mrc
import utils
import fft
import lie_tools
import config
from lattice import Lattice
from models import HetOnlyVAE

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('weights', help='Model weights')
    parser.add_argument('-c', '--config', metavar='PKL', help='CryoDRGN configuration')
    parser.add_argument('-z', type=np.float32, nargs='*', help='')
    parser.add_argument('--z-start', type=np.float32, nargs='*', help='')
    parser.add_argument('--z-end', type=np.float32, nargs='*', help='')
    parser.add_argument('-n', type=int, default=10, help='Number of structures between z_start and z_end')
    parser.add_argument('--zfile', help='Text file witth values of z to evaluate')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output .mrc or directory')
    parser.add_argument('--Apix', type=float, default=1, help='Pixel size to add to .mrc header (default: %(default)s A/pix)')
    parser.add_argument('--flip', action='store_true', help='Flip handedness of output volume')
    parser.add_argument('--prefix', default='vol_', help='Prefix when writing out multiple .mrc files (default: %(default)s)')
    parser.add_argument('-d','--downsample', type=int, help='Optionally downsample volumes to this box size')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')

    group = parser.add_argument_group('Overwrite architecture hyperparameters in config.pkl')
    group.add_argument('--norm', nargs=2, type=float)
    group.add_argument('-D', type=int, help='Box size')
    group.add_argument('--qlayers', type=int, help='Number of hidden layers')
    group.add_argument('--qdim', type=int, help='Number of nodes in hidden layers')
    group.add_argument('--zdim', type=int,  help='Dimension of latent variable')
    group.add_argument('--encode-mode', choices=('conv','resid','mlp','tilt'), help='Type of encoder network')
    group.add_argument('--players', type=int, help='Number of hidden layers')
    group.add_argument('--pdim', type=int, help='Number of nodes in hidden layers')
    group.add_argument('--enc-mask', type=int, help='Circular mask radius for image encoder')
    group.add_argument('--pe-type', choices=('geom_ft','geom_full','geom_lowf','geom_nohighf','linear_lowf','none'), help='Type of positional encoding')
    group.add_argument('--domain', choices=('hartley','fourier'))
    parser.add_argument('--l-extent', type=float, help='Coordinate lattice size')
    return parser

def main(args):
    t1 = dt.now()

    ## set the device
    use_cuda = torch.cuda.is_available()
    log('Use cuda {}'.format(use_cuda))
    if use_cuda:
        torch.set_default_tensor_type(torch.cuda.FloatTensor)

    if args.config is not None:
        args = config.load_config(args.config, args)
    log(args)

    if args.downsample:
        assert args.downsample%2==0
        assert args.downsample < args.D, "Must be smaller than original box size"
    D = args.D + 1
    lattice = Lattice(D, extent=args.l_extent)
    if args.enc_mask: 
        args.enc_mask = lattice.get_circular_mask(args.enc_mask)
        in_dim = args.enc_mask.sum()
    else:
        in_dim = lattice.D**2
    model = HetOnlyVAE(lattice, args.qlayers, args.qdim, args.players, args.pdim,
                in_dim, args.zdim, encode_mode=args.encode_mode, enc_mask=args.enc_mask, enc_type=args.pe_type, domain=args.domain)

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
            if args.downsample:
                extent = lattice.extent * (args.downsample/args.D)
                vol = model.decoder.eval_volume(lattice.get_downsample_coords(args.downsample+1), 
                                                args.downsample+1, extent, args.norm, zz) 
            else:
                vol = model.decoder.eval_volume(lattice.coords, lattice.D, lattice.extent, args.norm, zz) 
            out_mrc = '{}/{}{:03d}.mrc'.format(args.o, args.prefix, i)
            if args.flip:
                vol = vol[::-1]
            mrc.write(out_mrc, vol.astype(np.float32), Apix=args.Apix)

    else:
        z = np.array(args.z)
        log(z)
        if args.downsample:
            extent = lattice.extent * (args.downsample/args.D)
            vol = model.decoder.eval_volume(lattice.get_downsample_coords(args.downsample+1), 
                                            args.downsample+1, extent, args.norm, z) 
        else:
            vol = model.decoder.eval_volume(lattice.coords, lattice.D, lattice.extent, args.norm, z) 
        if args.flip:
            vol = vol[::-1]
        mrc.write(args.o, vol.astype(np.float32), Apix=args.Apix)

    td = dt.now()-t1
    log('Finsihed in {}'.format(td))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)

