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

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('particles', help='Particle stack file (.mrc)')
    parser.add_argument('metadata', help='Eman euler angles (.pkl)')
    parser.add_argument('-o', '--outdir', type=os.path.abspath, required=True, help='Output directory to save model')
    parser.add_argument('-d', '--device', type=int, default=-2, help='Compute device to use')
    parser.add_argument('--load-weights', type=os.path.abspath, help='Initialize network with existing weights')

    group = parser.add_argument_group('NN parameters')
    group.add_argument('-n', '--num-epochs', type=int, default=10, help='Number of training epochs (default: %(default)s)')
    group.add_argument('--dim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    group.add_argument('--wd', type=float, default=0, help='Weight decay in Adam optimizer (default: %(default)s)')
    group.add_argument('--lr', type=float, default=1e-5, help='Learning rate in Adam optimizer (default: %(default)s)')
    return parser

class ResidLinear(nn.Module):
    def __init__(self, nin, nout):
        super(ResidLinear, self).__init__()
        self.linear = nn.Linear(nin, nout)

    def forward(self, x):
        z = self.linear(x) + x
        return z

def R_from_eman(a,b,y):
    a *= np.pi/180.
    b *= np.pi/180.
    y *= np.pi/180.
    ca, sa = np.cos(a), np.sin(a)
    cb, sb = np.cos(b), np.sin(b)
    cy, sy = np.cos(y), np.sin(y)
    Ra = np.array([[ca,-sa,0],[sa,ca,0],[0,0,1]])
    Rb = np.array([[1,0,0],[0,cb,-sb],[0,sb,cb]])
    Ry = np.array(([cy,-sy,0],[sy,cy,0],[0,0,1]))
    return np.dot(np.dot(Ry,Rb),Ra)

def load_angles(pkl):
    with open(pkl,'rb') as f:
        ang = pickle.load(f)
    return ang

def fft2_center(img):
    '''input: centered realspace image (M x N)
    output: (M x N x 2) real and imag components of image ft'''
    ff = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(img)))
    return np.stack([ff.real, ff.imag], 2).astype(np.float32)

def eval_volume(model, nz, ny, nx, x_coord, use_cuda, rnorm, inorm):
    vol_f = np.zeros((nz,ny,nx),dtype=complex)
    assert not model.training
    # evaluate the volume by zslice to avoid memory overflows
    for i, z in enumerate(np.linspace(-1,1,nz,endpoint=False)):
        x = Variable(torch.from_numpy(x_coord + np.array([0,0,z],dtype=np.float32)))
        if use_cuda: x = x.cuda()
        with torch.no_grad():
            y = model(x)
            y = y.view(ny, nx, 2).cpu().numpy()
        vol_f[i] = (y[:,:,0]*rnorm[1]+rnorm[0]) + 1j*(y[:,:,1]*inorm[1]+inorm[0])
    vol = np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(vol_f)))
    vol = np.asarray([x[::-1] for x in vol])
    return vol, vol_f

def weighted_mse_loss(y_hat, y, w):
    '''
    y_hat: FloatTensor of dim ny x nx x 2
    y: FloatTensor of dim ny x nx x 2
    w: FloatTensor of dim ny x nx
    '''
    return (((y_hat-y) ** 2).sum(dim=2) * w).sum() / w.sum()

def main(args):
    t1 = dt.now()
    if args.outdir is not None and not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # load the particles
    with open(args.particles, 'rb') as f:
        content = f.read()
    particles, _, _ = mrc.parse(content)
    Nimg, ny, nx = particles.shape
    log('Loaded {} {}x{} images'.format(Nimg, ny, nx))
    particles = np.stack([fft2_center(img[::-1]) for img in particles],0)
    assert particles.shape == (Nimg,ny,nx,2)
    rnorm  = np.mean(particles[:,:,:,0]), np.std(particles[:,:,:,0])
    inorm = np.mean(particles[:,:,:,1]), np.std(particles[:,:,:,1])
    log('Normalizing by mean, std of real component: {} +/- {}'.format(*rnorm))
    log('Normalizing by mean, std of imag component: {} +/- {}'.format(*inorm))
    particles[:,:,:,0] = (particles[:,:,:,0]-rnorm[0])/rnorm[1]
    particles[:,:,:,1] = (particles[:,:,:,1]-inorm[0])/inorm[1]

    # load the metadata
    eman_eulers = load_angles(args.metadata)
    assert len(eman_eulers) == len(particles)

    # centered and scaled xy plane, values between -1 and 1
    # use a left-handed coordinate system -- this is necessary for compatibility
    # with eman2 projection data.
    nz = max(nx,ny)
    x1, x0 = np.meshgrid(np.linspace(-1, 1, ny, endpoint=False), np.linspace(-1, 1, nx, endpoint=False), indexing='ij')
    x_coord = np.stack([x0.ravel(),x1.ravel(),np.zeros(ny*nx)],1).astype(np.float32)

    # weights for weighted mse loss
    x_w = x0**2 + x1**2
    x_w[int(ny/2), int(nx/2)] = x_w[int(ny/2), int(nx/2)+1] # replace the 0 weight at the origin
    x_w = Variable(torch.from_numpy(x_w.astype(np.float32)))

    Rs = [R_from_eman(*ang).astype(np.float32) for ang in eman_eulers]

    hidden_dim = args.dim # 512
    activation = nn.ReLU
    model = nn.Sequential( nn.Linear(3, hidden_dim)
                         #, nn.BatchNorm1d(hidden_dim)
                         , activation()
                         , ResidLinear(hidden_dim, hidden_dim)
                         , activation()
                         , ResidLinear(hidden_dim, hidden_dim)
                         , activation()
                         , ResidLinear(hidden_dim, hidden_dim)
                         , activation()
                         , ResidLinear(hidden_dim, hidden_dim)
                         , activation()
                         , ResidLinear(hidden_dim, hidden_dim)
                         , activation()
                         , ResidLinear(hidden_dim, hidden_dim)
                         , activation()
                         , ResidLinear(hidden_dim, hidden_dim)
                         , activation()
                         , ResidLinear(hidden_dim, hidden_dim)
                         , activation()
                         , ResidLinear(hidden_dim, hidden_dim)
                         , activation()
                         , ResidLinear(hidden_dim, hidden_dim)
                         , activation()
                         , nn.Linear(hidden_dim, 2)
                         )

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
        x_w = x_w.cuda()

    optim = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.wd)

    num_epochs = args.num_epochs
    for epoch in range(num_epochs):
        loss_accum = 0
        ii = 0
        for i in np.random.permutation(Nimg):
            rot = Rs[i]
            x = Variable(torch.from_numpy(np.dot(x_coord, rot))) # this is R.T * x, if x is a column vec
            y = Variable(torch.from_numpy(particles[i]))
            if use_cuda:
                x = x.cuda()
                y = y.cuda()
            y_hat = model(x)
            y_hat = y_hat.view(ny, nx, 2)
            #loss = weighted_mse_loss(y_hat, y, x_w)
            loss = F.mse_loss(y_hat, y)

            loss.backward()
            optim.step()
            optim.zero_grad()

            loss_accum += loss.item()
            if ii > 0 and ii % 100 == 0:
                log('# [{}/{} epochs] [{}/{} images] training mse={:.8f}'.format(epoch+1, num_epochs, ii, Nimg, loss_accum/100))
                loss_accum = 0
            ii += 1

        model.eval()
        vol, vol_f = eval_volume(model, nz, ny, nx, x_coord, use_cuda, rnorm, inorm)
        with open('{}/reconstruct.{}.mrc'.format(args.outdir,epoch), 'wb') as f:
            mrc.write(f, vol.astype(np.float32))
        with open('{}/reconstruct.{}.pkl'.format(args.outdir,epoch), 'wb') as f:
            pickle.dump(vol_f, f)
        model.train()

    ## save model weights and evaluate the image model 
    model.eval()
    with open('{}/weights.pkl'.format(args.outdir), 'wb') as f:
        pickle.dump(model.state_dict(), f)

    vol, vol_f = eval_volume(model, nz, ny, nx, x_coord, use_cuda, rnorm, inorm)
    with open('{}/reconstruct.mrc'.format(args.outdir), 'wb') as f:
        mrc.write(f, vol.astype(np.float32))
    with open('{}/reconstruct.pkl'.format(args.outdir,epoch), 'wb') as f:
        pickle.dump(vol_f, f)
    
    log('Finsihed in {}'.format(dt.now()-t1))

if __name__ == '__main__':
    main(parse_args().parse_args())
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)
