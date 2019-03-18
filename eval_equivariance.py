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
from models import BNBHetOpt, HetVAE
from losses import EquivarianceLoss

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('particles', help='Particle stack file (.mrc)')
    parser.add_argument('weights', help='Model weights')
    parser.add_argument('-N', type=int, help='First N images (default: all images)')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output pickle')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')
    parser.add_argument('--ref',type=int, default=0,help='Reference image to align by')

    group = parser.add_argument_group('Architecture parameters')
    group.add_argument('-b','--batch-size', type=int, default=100, help='Minibatch size (default: %(default)s)')
    group.add_argument('--qlayers', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--qdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--zdim', type=int, default=1, help='Dimension of latent variable')
    group.add_argument('--encode-mode', default='resid', choices=('conv','resid','mlp'), help='Type of encoder network')
    group.add_argument('--players', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--pdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    return parser

def loss_function(recon_y, y, w_eps, z_std, mu, logvar):
    gen_loss = F.mse_loss(recon_y, y)  
    cross_entropy = torch.tensor([np.log(8*np.pi**2)], device=y.device) # cross entropy between gaussian and uniform on SO3
    entropy = lie_tools.so3_entropy(w_eps,z_std)
    kld_rot = cross_entropy - entropy
    kld_conf = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())
    return gen_loss, kld_rot.mean(), kld_conf

def main(args):
    log(args)
    t1 = dt.now()

    ## set the device
    use_cuda = torch.cuda.is_available()
    log('Use cuda {}'.format(use_cuda))
    if use_cuda:
        torch.set_default_tensor_type(torch.cuda.FloatTensor)

    # load the particles
    particles_real, _, _ = mrc.parse_mrc(args.particles,lazy=True)
    particles_real = np.asarray([x.get() for x in particles_real])
    particles_real = particles_real.astype(np.float32)
    Nimg, ny, nx = particles_real.shape
    nz = max(nx,ny)
    log('Loaded {} {}x{} images'.format(Nimg, ny, nx))
    particles_ft = np.asarray([fft.ht2_center(img).astype(np.float32) for img in particles_real])
    assert particles_ft.shape == (Nimg,ny,nx)
    rnorm  = [np.mean(particles_real), np.std(particles_real)]
    rnorm[0] = 0
    rnorm[1] = np.median([np.max(x) for x in particles_real])
    particles_real = (particles_real - rnorm[0])/rnorm[1]
    log('Normalizing particles by mean, std: {} +/- {}'.format(*rnorm))

    rnorm  = [np.mean(particles_ft), np.std(particles_ft)]
    log('Particle stack mean, std: {} +/- {}'.format(*rnorm))
    rnorm[0] = 0
    log('Normalizing FT by mean, std: {} +/- {}'.format(*rnorm))
    particles_ft = (particles_ft - rnorm[0])/rnorm[1]

    if args.N:
        log('Using first {} images'.format(args.N))
        particles_real = particles_real[:args.N]
        particles_ft = particles_ft[:args.N]
        Nimg = args.N

    model = HetVAE(nx, ny, args.qlayers, args.qdim, args.players, args.pdim,
                args.zdim, encode_mode=args.encode_mode)
    bnb = BNBHetOpt(model,ny,nx)
    equiv = EquivarianceLoss(model, ny, nx)

    log('Loading weights from {}'.format(args.weights))
    checkpoint = torch.load(args.weights)
    model.load_state_dict(checkpoint['model_state_dict'])

    model.eval()

    recon_all = []
    z_all = []
    z_rot_all = []

    num_batches = np.ceil(Nimg / args.batch_size).astype(int)
    for minibatch_i in np.array_split(np.arange(Nimg),num_batches):
        # inference with real space image
        y = torch.from_numpy(particles_ft[minibatch_i])
        if use_cuda: 
            y = y.cuda()
        mu, logvar = model.encode(y) 

        n = len(minibatch_i)
        theta = torch.rand(n)*2*np.pi
        y = torch.unsqueeze(y, 1)
        img_rot = equiv.rotate(y, theta)
        img_rot = torch.squeeze(img_rot)
        mu_rot, _ = model.encode(img_rot)

        z_all.append(mu.detach().cpu().numpy())
        z_rot_all.append(mu_rot.detach().cpu().numpy())

    if args.o:
        z_all = np.vstack(z_all)
        z_rot_all = np.vstack(z_rot_all)
        log(np.mean(z_all))
        with open(args.o,'wb') as f:
            pickle.dump(z_all, f)
            pickle.dump(z_rot_all, f)

    td = dt.now()-t1
    log('Finsihed in {}'.format(td))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)

