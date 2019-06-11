'''
Compute equivariance loss on an image stack
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
import dataset
from lattice import Lattice
from models import HetOnlyVAE
from losses import EquivarianceLoss

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('particles', help='Particle stack file (.mrc)')
    parser.add_argument('weights', help='Model weights')
    parser.add_argument('--norm', type=float, nargs=2, default=None, help='Data normalization as shift, 1/scale (default: mean, std of dataset)')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output pickle')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')

    group = parser.add_argument_group('Architecture parameters')
    group.add_argument('-b','--batch-size', type=int, default=100, help='Minibatch size (default: %(default)s)')
    group.add_argument('--qlayers', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--qdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--zdim', type=int, default=1, help='Dimension of latent variable')
    group.add_argument('--encode-mode', default='resid', choices=('conv','resid','mlp'), help='Type of encoder network')
    group.add_argument('--players', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--pdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--enc-mask', type=int, help='Circulask mask of image for encoder')
    return parser

def main(args):
    log(args)
    t1 = dt.now()

    ## set the device
    use_cuda = torch.cuda.is_available()
    log('Use cuda {}'.format(use_cuda))
    if use_cuda:
        torch.set_default_tensor_type(torch.cuda.FloatTensor)

    # load the particles
    data = dataset.MRCData(args.particles, norm=args.norm)
    D = data.D 
    Nimg = data.N

    lattice = Lattice(D)
    if args.enc_mask: args.enc_mask = lattice.get_circular_mask(args.enc_mask)
    model = HetOnlyVAE(lattice, args.qlayers, args.qdim, args.players, args.pdim,
                args.zdim, encode_mode=args.encode_mode, enc_mask=args.enc_mask)
    equiv = EquivarianceLoss(model, D)

    log('Loading weights from {}'.format(args.weights))
    checkpoint = torch.load(args.weights)
    model.load_state_dict(checkpoint['model_state_dict'])

    model.eval()
    z_all = []
    z_rot_all = []
    data_iterator = DataLoader(data, batch_size=args.batch_size, shuffle=False)
    for batch in data_iterator:
        y = batch[0].to(device)
        with torch.no_grad():
            mu, logvar = model.encode(y) 
        theta = torch.rand(len(y))*2*np.pi
        y = torch.unsqueeze(y, 1)
        img_rot = equiv.rotate(y, theta)
        img_rot = torch.squeeze(img_rot)
        with torch.no_grad():
            mu_rot, _ = model.encode(img_rot)
        z_all.append(mu.detach().cpu().numpy())
        z_rot_all.append(mu_rot.detach().cpu().numpy())
    z_all = np.vstack(z_all)
    z_rot_all = np.vstack(z_rot_all)
    log(np.mean(z_all))
    log(np.mean(z_rot_all))
    with open(args.o,'wb') as f:
        pickle.dump(z_all, f)
        pickle.dump(z_rot_all, f)

    td = dt.now()-t1
    log('Finsihed in {}'.format(td))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)

