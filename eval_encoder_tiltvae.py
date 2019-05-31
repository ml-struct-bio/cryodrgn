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
from models import TiltVAE
from lattice import Lattice
import dataset

from scipy.linalg import logm
def geodesic_so3(A,B):
    return np.sum(logm(np.dot(A.T,B))**2)**.5

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('particles', help='Particle stack file (.mrc)')
    parser.add_argument('particles_tilt', help='Particle stack file (.mrc)')
    parser.add_argument('weights', help='Particle stack file (.mrc)')
    parser.add_argument('--tilt', type=float, default=45, help='Right-handed x-axis tilt offset in degrees (default: %(default)s)')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output pickle')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')
    parser.add_argument('--ref',type=int, default=0,help='Reference image to align by')
    parser.add_argument('--angles',help='euler angles (.pkl)')
    parser.add_argument('--flip-hand', action='store_true', help='Flip hand of reference euler angles')
    parser.add_argument('--save-recon')

    group = parser.add_argument_group('Training parameters')
    group.add_argument('-b','--batch-size', type=int, default=100, help='Minibatch size (default: %(default)s)')
    group.add_argument('--qlayers', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--qdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--players', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--pdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--enc-mask', type=int, help='Circulask mask of image for encoder')
    group.add_argument('--no-trans', action='store_true', help="Don't translate images")
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
    assert args.no_trans, "Not implemented"
    t1 = dt.now()

    ## set the device
    use_cuda = torch.cuda.is_available()
    device = torch.device('cuda' if use_cuda else 'cpu')
    log('Use cuda {}'.format(use_cuda))
    if use_cuda:
        torch.set_default_tensor_type(torch.cuda.FloatTensor)

    # load the particles
    data = dataset.TiltMRCData(args.particles, args.particles_tilt)
    Nimg = data.N
    D = data.D

    theta = args.tilt*np.pi/180
    tilt = np.array([[1.,0.,0.],
                    [0, np.cos(theta), -np.sin(theta)],
                    [0, np.sin(theta), np.cos(theta)]]).astype(np.float32)

    lattice = Lattice(D)
    if args.enc_mask: args.enc_mask = lattice.get_circular_mask(args.enc_mask)
    model = TiltVAE(lattice, tilt, args.qlayers, args.qdim, args.players, args.pdim, no_trans=args.no_trans, enc_mask=args.enc_mask)

    log('Loading weights from {}'.format(args.weights))
    checkpoint = torch.load(args.weights)
    model.load_state_dict(checkpoint['model_state_dict'])

    model.eval()

    rot_all = []
    data_iterator = DataLoader(data, batch_size=args.batch_size, shuffle=False)
    for batch in data_iterator:
        y = batch[0].to(device)
        yt = batch[1].to(device)
        rot, _, _, _, _, _, _ = model.encode(y,yt)
        rot_all.append(rot.detach().cpu().numpy())

    rot_all = np.vstack(rot_all)
    if args.o:
        with open(args.o,'wb') as f:
            pickle.dump(rot_all, f)
    td = dt.now()-t1
    log('Finsihed in {}'.format(td))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)

