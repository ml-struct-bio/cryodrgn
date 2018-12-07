'''
Train a NN to model a 3D EM density map
    given 2D projections with angular assignments
'''
import numpy as np
import sys, os
import argparse
import pickle
from datetime import datetime as dt

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable

sys.path.insert(0,os.path.abspath(os.path.dirname(__file__))+'/lib-python')
import mrc
import utils
import fft

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('particles', help='Particle stack file (.mrc)')
    parser.add_argument('metadata', help='Eman euler angles (.pkl)')
    parser.add_argument('-o', '--outdir', type=os.path.abspath, required=True, help='Output directory to save model')
    parser.add_argument('-d', '--device', type=int, default=-2, help='Compute device to use')
    parser.add_argument('--load-weights', type=os.path.abspath, help='Initialize network with existing weights')
    parser.add_argument('--save-intermediates', action='store_true', help='Save out structure each epoch')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')

    group = parser.add_argument_group('NN parameters')
    group.add_argument('-n', '--num-epochs', type=int, default=10, help='Number of training epochs (default: %(default)s)')
    group.add_argument('--layers', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--dim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--wd', type=float, default=0, help='Weight decay in Adam optimizer (default: %(default)s)')
    group.add_argument('--lr', type=float, default=1e-3, help='Learning rate in Adam optimizer (default: %(default)s)')
    return parser

class ResidLinear(nn.Module):
    def __init__(self, nin, nout):
        super(ResidLinear, self).__init__()
        self.linear = nn.Linear(nin, nout)

    def forward(self, x):
        z = self.linear(x) + x
        return z

def define_model(nlayers, hidden_dim, activation):
    layers = [nn.Linear(3, hidden_dim), activation()]
    for n in range(nlayers):
        layers.append(ResidLinear(hidden_dim, hidden_dim))
        layers.append(activation())
    layers.append(nn.Linear(hidden_dim,1))
    return nn.Sequential(*layers)

def eval_volume(model, nz, ny, nx, x_coord, use_cuda, rnorm):
    vol_f = np.zeros((nz,ny,nx),dtype=complex)
    assert not model.training
    # evaluate the volume by zslice to avoid memory overflows
    for i, z in enumerate(np.linspace(-1,1,nz,endpoint=False)):
        x = Variable(torch.from_numpy(x_coord + np.array([0,0,z],dtype=np.float32)))
        if use_cuda: x = x.cuda()
        with torch.no_grad():
            y = model(x)
            y = y.view(ny, nx).cpu().numpy()
        vol_f[i] = y*rnorm[1]+rnorm[0]
    vol = fft.ihtn_center(vol_f)
    return vol, vol_f

def main(args):
    t1 = dt.now()
    if args.outdir is not None and not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # load the particles
    particles, _, _ = mrc.parse_mrc(args.particles)
    Nimg, ny, nx = particles.shape
    log('Loaded {} {}x{} images'.format(Nimg, ny, nx))
    particles = np.stack([fft.ht2_center(img).astype(np.float32) for img in particles],0)
    assert particles.shape == (Nimg,ny,nx)
    rnorm  = np.mean(particles), np.std(particles)
    log('Normalizing by mean, std: {} +/- {}'.format(*rnorm))
    particles -= rnorm[0]
    particles /= rnorm[1]

    # load the metadata
    eman_eulers = utils.load_angles(args.metadata)
    assert len(eman_eulers) == len(particles)

    # centered and scaled xy plane, values between -1 and 1
    nz = max(nx,ny)
    x0, x1 = np.meshgrid(np.linspace(-1, 1, nx, endpoint=False), np.linspace(-1, 1, ny, endpoint=False))
    x_coord = np.stack([x0.ravel(),x1.ravel(),np.zeros(ny*nx)],1).astype(np.float32)

    model = define_model(args.layers, args.dim, nn.ReLU)

    if args.load_weights:
        log('Initializing weights from {}'.format(args.load_weights))
        with open(args.load_weights,'rb') as f:
            model.load_state_dict(pickle.load(f))

    ## set the device
    d = args.device
    use_cuda = (d != -1) and torch.cuda.is_available()
    if d >= 0:
        torch.cuda.set_device(d)
    log('Use cuda {}'.format(use_cuda))
    if use_cuda:
        model.cuda()

    optim = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.wd)

    num_epochs = args.num_epochs
    for epoch in range(num_epochs):
        loss_accum = 0
        ii = 0
        for i in np.random.permutation(Nimg):
            rot = utils.R_from_eman(*eman_eulers[i]).astype(np.float32)
            x = Variable(torch.from_numpy(np.dot(x_coord, rot))) 
            y = Variable(torch.from_numpy(particles[i]))
            if use_cuda:
                x = x.cuda()
                y = y.cuda()
            y_hat = model(x)
            y_hat = y_hat.view(ny, nx)
            loss = F.mse_loss(y_hat, y)

            loss.backward()
            optim.step()
            optim.zero_grad()

            loss_accum += loss.item()
            if ii > 0 and ii % 100 == 0:
                log('# [{}/{} epochs] [{}/{} images] training mse={:.8f}'.format(epoch+1, num_epochs, ii, Nimg, loss_accum/100))
                loss_accum = 0
            ii += 1

        if args.save_intermediates:
            model.eval()
            vol, vol_f = eval_volume(model, nz, ny, nx, x_coord, use_cuda, rnorm)
            with open('{}/reconstruct.{}.mrc'.format(args.outdir,epoch), 'wb') as f:
                mrc.write(f, vol.astype(np.float32))
            with open('{}/reconstruct.{}.pkl'.format(args.outdir,epoch), 'wb') as f:
                pickle.dump(vol_f, f)
            model.train()

    ## save model weights and evaluate the image model 
    model.eval()
    with open('{}/weights.pkl'.format(args.outdir), 'wb') as f:
        pickle.dump(model.state_dict(), f)

    vol, vol_f = eval_volume(model, nz, ny, nx, x_coord, use_cuda, rnorm)
    mrc.write('{}/reconstruct.mrc'.format(args.outdir), vol.astype(np.float32))
    with open('{}/reconstruct.pkl'.format(args.outdir), 'wb') as f:
        pickle.dump(vol_f, f)
    
    td = dt.now()-t1
    log('Finsihed in {} ({} per epoch)'.format(td, td/num_epochs))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)
