'''
Evaluate VAE on an image stack
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
    parser.add_argument('particles', help='Particle stack file (.mrc)')
    parser.add_argument('particles_tilt', help='Particle stack file (.mrc)')
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
    group.add_argument('--encode-mode', default='resid', choices=('conv','resid','mlp','tilt'), help='Type of encoder network')
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
    particles, _, _ = mrc.parse_mrc(args.particles)
    particles = particles.astype(np.float32)
    Nimg, ny, nx = particles.shape
    nz = max(nx,ny)
    log('Loaded {} {}x{} images'.format(Nimg, ny, nx))
    particles = np.asarray([fft.ht2_center(img).astype(np.float32) for img in particles])
    rnorm  = [np.mean(particles), np.std(particles)]
    log('Particle FT stack mean, std: {} +/- {}'.format(*rnorm))
    rnorm[0] = 0
    log('Normalizing FT by mean, std: {} +/- {}'.format(*rnorm))
    particles = (particles - rnorm[0])/rnorm[1]

    # load particles from tilt series
    particles_tilt, _, _ = mrc.parse_mrc(args.particles_tilt)
    assert particles_tilt.shape == (Nimg, ny, nx), 'Tilt series pair must have same dimensions as untilted particles'
    particles_tilt = np.asarray([fft.ht2_center(img).astype(np.float32) for img in particles_tilt])
    particles_tilt = (particles_tilt - rnorm[0])/rnorm[1]

    if args.N:
        log('Using first {} images'.format(args.N))
        particles = particles[:args.N]
        particles_tilt = particles_tilt[:args.N]
        Nimg = args.N

    lattice = Lattice(nx)
    model = HetOnlyVAE(lattice, ny*nx, args.qlayers, args.qdim, args.players, args.pdim,
                args.zdim, encode_mode=args.encode_mode)

    log('Loading weights from {}'.format(args.weights))
    checkpoint = torch.load(args.weights)
    model.load_state_dict(checkpoint['model_state_dict'])

    model.eval()

    recon_all = []
    z_all = []

    num_batches = np.ceil(Nimg / args.batch_size).astype(int)
    for minibatch_i in np.array_split(np.arange(Nimg),num_batches):
        # inference with real space image
        y = torch.from_numpy(particles[minibatch_i])
        yt = torch.from_numpy(particles_tilt[minibatch_i])
        if use_cuda: 
            y = y.cuda()
            yt = yt.cuda()
        mu, logvar = model.encode((y,yt))

        z_all.append(mu.detach().cpu().numpy())

    if args.o:
        z_all = np.vstack(z_all)
        log(np.mean(z_all))
        with open(args.o,'wb') as f:
            pickle.dump(z_all, f)

    td = dt.now()-t1
    log('Finsihed in {}'.format(td))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)

