'''
Evaluate the decoder at specified values of z
'''
import numpy as np
import sys, os
import argparse
import pickle
from datetime import datetime as dt
import matplotlib.pyplot as plt 
import pprint

import torch

from cryodrgn import mrc
from cryodrgn import utils
from cryodrgn import fft
from cryodrgn import lie_tools
from cryodrgn import config
from cryodrgn.lattice import Lattice
from cryodrgn.models import HetOnlyVAE

log = utils.log
vlog = utils.vlog

def add_args(parser):
    parser.add_argument('weights', help='Model weights')
    parser.add_argument('-c', '--config', metavar='PKL', required=True, help='CryoDRGN config.pkl file')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output .mrc or directory')
    parser.add_argument('--prefix', default='vol_', help='Prefix when writing out multiple .mrc files (default: %(default)s)')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')

    group = parser.add_argument_group('Specify z values')
    group.add_argument('-z', type=np.float32, nargs='*', help='Specify one z-value')
    group.add_argument('--z-start', type=np.float32, nargs='*', help='Specify a starting z-value')
    group.add_argument('--z-end', type=np.float32, nargs='*', help='Specify an ending z-value')
    group.add_argument('-n', type=int, default=10, help='Number of structures between [z_start, z_end]')
    group.add_argument('--zfile', help='Text file with z-values to evaluate')

    group = parser.add_argument_group('Volume arguments')
    group.add_argument('--Apix', type=float, default=1, help='Pixel size to add to .mrc header (default: %(default)s A/pix)')
    group.add_argument('--flip', action='store_true', help='Flip handedness of output volume')
    group.add_argument('--invert', action='store_true', help='Invert contrast of output volume')
    group.add_argument('-d','--downsample', type=int, help='Downsample volumes to this box size (pixels)')

    group = parser.add_argument_group('Overwrite architecture hyperparameters in config.pkl')
    group.add_argument('--norm', nargs=2, type=float)
    group.add_argument('-D', type=int, help='Box size')
    group.add_argument('--enc-layers', dest='qlayers', type=int, help='Number of hidden layers')
    group.add_argument('--enc-dim', dest='qdim', type=int, help='Number of nodes in hidden layers')
    group.add_argument('--zdim', type=int,  help='Dimension of latent variable')
    group.add_argument('--encode-mode', choices=('conv','resid','mlp','tilt'), help='Type of encoder network')
    group.add_argument('--dec-layers', dest='players', type=int, help='Number of hidden layers')
    group.add_argument('--dec-dim', dest='pdim', type=int, help='Number of nodes in hidden layers')
    group.add_argument('--enc-mask', type=int, help='Circular mask radius for image encoder')
    group.add_argument('--pe-type', choices=('geom_ft','geom_full','geom_lowf','geom_nohighf','linear_lowf','none'), help='Type of positional encoding')
    group.add_argument('--feat-sigma', type=float, help="Scale for random Gaussian features")
    group.add_argument('--pe-dim', type=int, help='Num sinusoid features in positional encoding (default: D/2)')
    group.add_argument('--domain', choices=('hartley','fourier'))
    group.add_argument('--l-extent', type=float, help='Coordinate lattice size')
    group.add_argument('--activation', choices=('relu','leaky_relu'), default='relu', help='Activation (default: %(default)s)')
    return parser

def check_inputs(args):
    if args.z_start:
        assert args.z_end, "Must provide --z-end with argument --z-start"
    assert sum((bool(args.z), bool(args.z_start), bool(args.zfile))) == 1, "Must specify either -z OR --z-start/--z-end OR --zfile"

def main(args):
    check_inputs(args)
    t1 = dt.now()

    ## set the device
    use_cuda = torch.cuda.is_available()
    device = torch.device('cuda' if use_cuda else 'cpu')
    log('Use cuda {}'.format(use_cuda))
    if not use_cuda:
        log('WARNING: No GPUs detected')

    log(args)
    cfg = config.overwrite_config(args.config, args)
    log('Loaded configuration:')
    pprint.pprint(cfg)

    D = cfg['lattice_args']['D'] # image size + 1
    zdim = cfg['model_args']['zdim']
    norm = cfg['dataset_args']['norm']

    if args.downsample:
        assert args.downsample % 2 == 0, "Boxsize must be even"
        assert args.downsample <= D - 1, "Must be smaller than original box size"
    
    model, lattice = HetOnlyVAE.load(cfg, args.weights, device=device)
    model.eval()

    ### Multiple z ###
    if args.z_start or args.zfile:
        
        ### Get z values
        if args.z_start:
            args.z_start = np.array(args.z_start)
            args.z_end = np.array(args.z_end)
            z = np.repeat(np.arange(args.n,dtype=np.float32), zdim).reshape((args.n, zdim))
            z *= ((args.z_end - args.z_start)/(args.n-1))
            z += args.z_start
        else:
            z = np.loadtxt(args.zfile).reshape(-1, zdim)

        if not os.path.exists(args.o):
            os.makedirs(args.o)

        log(f'Generating {len(z)} volumes')
        for i,zz in enumerate(z):
            log(zz)
            if args.downsample:
                extent = lattice.extent * (args.downsample/(D-1))
                vol = model.decoder.eval_volume(lattice.get_downsample_coords(args.downsample+1),
                                                args.downsample+1, extent, norm, zz)
            else:
                vol = model.decoder.eval_volume(lattice.coords, lattice.D, lattice.extent, norm, zz)
            out_mrc = '{}/{}{:03d}.mrc'.format(args.o, args.prefix, i)
            if args.flip:
                vol = vol[::-1]
            if args.invert:
                vol *= -1
            mrc.write(out_mrc, vol.astype(np.float32), Apix=args.Apix)

    ### Single z ###
    else:
        z = np.array(args.z)
        log(z)
        if args.downsample:
            extent = lattice.extent * (args.downsample/(D-1))
            vol = model.decoder.eval_volume(lattice.get_downsample_coords(args.downsample+1),
                                            args.downsample+1, extent, norm, z)
        else:
            vol = model.decoder.eval_volume(lattice.coords, lattice.D, lattice.extent, norm, z)
        if args.flip:
            vol = vol[::-1]
        if args.invert:
            vol *= -1
        mrc.write(args.o, vol.astype(np.float32), Apix=args.Apix)

    td = dt.now()-t1
    log('Finished in {}'.format(td))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    utils._verbose = args.verbose
    main(args)

