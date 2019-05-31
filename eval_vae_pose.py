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
from torch.utils.data import DataLoader

sys.path.insert(0,os.path.abspath(os.path.dirname(__file__))+'/lib-python')
import mrc
import utils
import fft
import lie_tools
import dataset

from models import VAE
from lattice import Lattice

from scipy.linalg import logm
def geodesic_so3(A,B):
    return np.sum(logm(np.dot(A.T,B))**2)**.5

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('particles', help='Particle stack file (.mrc)')
    parser.add_argument('weights', help='Model weights')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output pickle for rotations')
    parser.add_argument('--out-trans', type=os.path.abspath, help='Output pickle for translations (optional)')
    parser.add_argument('--tilt', help='Particle stack file (.mrcs)')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')

    group = parser.add_argument_group('Training parameters')
    group.add_argument('-b','--batch-size', type=int, default=100, help='Minibatch size (default: %(default)s)')
    group.add_argument('--qlayers', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--qdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--encode-mode', default='resid', choices=('conv','resid','mlp','tilt'), help='Type of encoder network (default: %(default)s)')
    group.add_argument('--enc-mask', type=int, help='Circulask mask of image for encoder')
    group.add_argument('--players', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--pdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    return parser

def loss_function(recon_y, y, w_eps, z_std):
    gen_loss = F.mse_loss(recon_y, y)  
    cross_entropy = torch.tensor([np.log(8*np.pi**2)], device=y.device) # cross entropy between gaussian and uniform on SO3
    entropy = lie_tools.so3_entropy(w_eps,z_std)
    kld = cross_entropy - entropy
    #assert kld > 0
    return gen_loss, kld.mean()

def fast_dist(a,b):
    return np.sum((a-b)**2)

def align_rot(rotA,rotB,i,flip=False):
    x = -1*np.ones((3,3))
    x[2] = 1
    if flip:
        rotB = [x*r for r in rotB]
                            
    ref = rotA[i]
    rot = np.array([np.dot(x, ref.T) for x in rotA])
    ref = rotB[i]
    rot_hat = np.array([np.dot(x, ref.T) for x in rotB])
    return rot, rot_hat


def main(args):
    log(args)
    t1 = dt.now()

    ## set the device
    use_cuda = torch.cuda.is_available()
    device = torch.device('cuda' if use_cuda else 'cpu')
    log('Use cuda {}'.format(use_cuda))
    if use_cuda:
        torch.set_default_tensor_type(torch.cuda.FloatTensor)

    # load the particles
    if args.tilt is None:
        data = dataset.MRCData(args.particles)
    else:
        data = dataset.TiltMRCData(args.particles, args.tilt)
    Nimg = data.N
    D = data.D

    lattice = Lattice(D)
    if args.enc_mask: args.enc_mask = lattice.get_circular_mask(args.enc_mask)
    model = VAE(lattice, args.qlayers, args.qdim, args.players, args.pdim,
                encode_mode=args.encode_mode, no_trans=not bool(args.out_trans), enc_mask=args.enc_mask)
    log('Loading weights from {}'.format(args.weights))
    checkpoint = torch.load(args.weights)
    model.load_state_dict(checkpoint['model_state_dict'])

    model.eval()

    rot_all = []
    trans_all = []
    recon_all = []
    data_generator = DataLoader(data, batch_size=args.batch_size, shuffle=True)

    for minibatch in data_generator:
        if args.tilt is None:
            minibatch = minibatch[0].to(device)
            z_mu, z_std, t_mu, t_logvar = model.encode(minibatch)
        else:
            y = minibatch[0].to(device)
            yt = minibatch[1].to(device)
            z_mu, z_std, t_mu, t_logvar = model.encode(y,yt)
        rot_all.append(z_mu.detach().cpu().numpy())
        if args.out_trans:
            trans_all.append(t_mu.detach().cpu().numpy())

    rot_all = np.vstack(rot_all)
    with open(args.o,'wb') as f:
        pickle.dump(rot_all, f)
    if args.out_trans:
        trans_all = np.vstack(trans_all)
        with open(args.out_trans,'wb') as f:
            pickle.dump(trans_all, f)

    td = dt.now()-t1
    log('Finsihed in {}'.format(td))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)

