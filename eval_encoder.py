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
from models import VAE

from scipy.linalg import logm
def geodesic_so3(A,B):
    return np.sum(logm(np.dot(A.T,B))**2)**.5

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('particles', help='Particle stack file (.mrc)')
    parser.add_argument('weights', help='Particle stack file (.mrc)')
    parser.add_argument('-N', type=int, default=1000, help='First N images (default: %(default)s)')
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
    group.add_argument('--encode-mode', default='resid', choices=('conv','resid','mlp'), help='Type of encoder network')
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
    log('Use cuda {}'.format(use_cuda))

    # load the particles
    particles_real, _, _ = mrc.parse_mrc(args.particles,lazy=True)
    particles_real = np.asarray([x.get() for x in particles_real[:args.N]])
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

    model = VAE(nx, ny, args.qlayers, args.qdim, args.players, args.pdim,
                encode_mode=args.encode_mode)
    if use_cuda:
        model.cuda()
        model.lattice = model.lattice.cuda()

    if args.weights:
        log('Loading weights from {}'.format(args.weights))
        checkpoint = torch.load(args.weights)
        model.load_state_dict(checkpoint['model_state_dict'])

    model.eval()
    torch.no_grad().__enter__()

    rot_all = []
    recon_all = []
    for epoch in range(1):
        ii = 0
        num_batches = np.ceil(Nimg / args.batch_size).astype(int)
        for minibatch_i in np.array_split(np.arange(Nimg),num_batches):
            ii += 1

            # inference with real space image
            y = Variable(torch.from_numpy(np.asarray([particles_real[i] for i in minibatch_i])))
            if use_cuda: y = y.cuda()

            rot, z_std = model.latent_encoder(model.encoder(y))

            flip = torch.tensor([[-1,-1,-1],[-1,-1,-1],[1,1,1]],dtype=torch.float32, device=y.device)
            rot2 = rot*flip
            align = lie_tools.s2s2_to_SO3(model.align[:3], model.align[3:])
            rot2 = rot2 @ align

            x = model.lattice @ rot # R.T*x
            y_hat = model.decoder(x)
            y_hat = y_hat.view(-1, ny, nx)
            
            x = model.lattice @ rot2 # R.T*x
            y_hat2 = model.decoder(x)
            y_hat2 = y_hat2.view(-1, ny, nx)
    
            y = Variable(torch.from_numpy(np.asarray([particles_ft[i] for i in minibatch_i])))
            if use_cuda: y = y.cuda()

            B = y.size(0)
            gen_loss = (y_hat - y).pow(2).view(B, -1).mean(-1)
            gen_loss_mirror = (y_hat2 - y).pow(2).view(B, -1).mean(-1)
            total = torch.min(gen_loss, gen_loss_mirror)
            log('# [Batch {}/{}] gen loss={:4f}'.format(ii, num_batches, total.mean().item()))

            rot=rot.detach().cpu().numpy()
            rot2=rot2.detach().cpu().numpy()
            gen_loss=gen_loss.detach().cpu().numpy()
            gen_loss_mirror=gen_loss_mirror.detach().cpu().numpy()
            minn=np.argmin(np.vstack((gen_loss,gen_loss_mirror)),axis=0)
            print(np.sum(minn))
            tmp=np.stack((rot,rot2))
            rot3 = np.array([tmp[x,i] for i,x in enumerate(minn)])
            rot_all.append(rot3)
            if args.save_recon:
                y_hat=y_hat.detach().cpu().numpy() 
                y_hat2=y_hat2.detach().cpu().numpy() 
                tmp = np.stack((y_hat,y_hat2))
                y_hat3 = np.array([tmp[x,i] for i,x in enumerate(minn)])
                recon_all.append(y_hat3)

    rot_all = np.vstack(rot_all)
    if args.o:
        with open(args.o,'wb') as f:
            pickle.dump(rot_all, f)
    if args.save_recon:
        recon_all = np.vstack(recon_all)
        with open(args.save_recon,'wb') as f:
            pickle.dump(recon_all, f)

    td = dt.now()-t1
    log('Finsihed in {}'.format(td))

    if args.angles:

        # compare angles with ground truth  
        angles = pickle.load(open(args.angles,'rb'))
        angles = angles[:Nimg]
    
        rot = [utils.R_from_eman(*x) for x in angles] # ground truth
        rot_hat = rot_all # predicted
        rot, rot_hat = align_rot(rot, rot_hat, args.ref, args.flip_hand)

        #dists = np.asarray([geodesic_so3(a,b) for a,b in zip(rot, rot_hat)])
        dists = np.asarray([fast_dist(a,b) for a,b in zip(rot, rot_hat)])

        print('Median: ', np.median(dists))
        w = np.where(np.asarray(dists)<1)
        print('Dist < 1: ', w[:10])

        with open(args.o+'.ind.pkl','wb') as f:
            pickle.dump(w[0], f)
        w = np.where(np.asarray(dists)>4)
        print('Dist > 4:', w[:10])

        i = w[0][-1]
        print('Image ',i)
        print('Ground truth rot: ',rot[i])
        print('Predicted rot: ',rot_hat[i])

        i = w[0][-2]
        print('Image ',i)
        print('Ground truth rot: ',rot[i])
        print('Predicted rot: ',rot_hat[i])
        plt.hist(dists,10)
        plt.show()

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)

