'''
VAE-based neural reconstruction with orientational inference

Ellen Zhong
12/7/2018
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
from torch.distributions import Normal

sys.path.insert(0,os.path.abspath(os.path.dirname(__file__))+'/lib-python')
import mrc
import utils
import fft
import lie_tools

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('particles', help='Particle stack file (.mrc)')
    parser.add_argument('-o', '--outdir', type=os.path.abspath, required=True, help='Output directory to save model')
    parser.add_argument('-d', '--device', type=int, default=-2, help='Compute device to use')
    parser.add_argument('--load-weights', type=os.path.abspath, help='Initialize network with existing weights')
    parser.add_argument('--save-intermediates', action='store_true', help='Save out structure each epoch')
    parser.add_argument('--log-interval', type=int, default=2, help='Logging interval in batches')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')

    group = parser.add_argument_group('Training parameters')
    group.add_argument('-n', '--num-epochs', type=int, default=10, help='Number of training epochs (default: %(default)s)')
    group.add_argument('-b','--batch-size', type=int, default=100, help='Minibatch size (default: %(default)s)')
    group.add_argument('--wd', type=float, default=0, help='Weight decay in Adam optimizer (default: %(default)s)')
    group.add_argument('--lr', type=float, default=1e-3, help='Learning rate in Adam optimizer (default: %(default)s)')
    group.add_argument('--beta', type=float, default=1.0, help='Weight of latent loss (default: %(default)s)')

    group = parser.add_argument_group('Encoder Network')
    group.add_argument('--qlayers', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--qdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')

    group = parser.add_argument_group('Decoder Network')
    group.add_argument('--players', type=int, default=10, help='Number of hidden layers (default: %(default)s)')
    group.add_argument('--pdim', type=int, default=128, help='Number of nodes in hidden layers (default: %(default)s)')
    return parser

class ResidLinear(nn.Module):
    def __init__(self, nin, nout):
        super(ResidLinear, self).__init__()
        self.linear = nn.Linear(nin, nout)

    def forward(self, x):
        z = self.linear(x) + x
        return z

class VAE(nn.Module):
    def __init__(self, nx, ny, encode_layers, encode_dim, decode_layers, decode_dim):
        super(VAE, self).__init__()
        self.nx = nx
        self.ny = ny
        self.in_dim = nx*ny
        self.encoder = self.get_encoder(nx*ny, encode_layers, encode_dim, nn.ReLU) #in_dim -> R9
        self.decoder = self.get_decoder(decode_layers, decode_dim, nn.ReLU) #R3 -> R1
        
        # centered and scaled xy plane, values between -1 and 1
        x0, x1 = np.meshgrid(np.linspace(-1, 1, nx, endpoint=False), np.linspace(-1, 1, ny, endpoint=False))
        lattice = np.stack([x0.ravel(),x1.ravel(),np.zeros(ny*nx)],1).astype(np.float32)
        self.lattice = torch.from_numpy(lattice)
    
    def get_encoder(self, in_dim, nlayers, hidden_dim, activation):
        '''
        Return a NN model mapping an image to mean and covariance of an element in so3
        '''
        layers = [nn.Linear(in_dim, hidden_dim), activation()]
        for n in range(nlayers):
            layers.append(ResidLinear(hidden_dim, hidden_dim))
            layers.append(activation())
        layers.append(nn.Linear(hidden_dim,9)) 
        return nn.Sequential(*layers)
    
    def get_decoder(self, nlayers, hidden_dim, activation):
        '''
        Return a NN model mapping cartesian coordinates to electron density
        (represented in Hartley reciprocal space)
        '''
        layers = [nn.Linear(3, hidden_dim), activation()]
        for n in range(nlayers):
            layers.append(ResidLinear(hidden_dim, hidden_dim))
            layers.append(activation())
        layers.append(nn.Linear(hidden_dim,1))
        return nn.Sequential(*layers)

    def reparameterize(self, z):
        '''
        Reparameterize SO(3) latent variable
        # z represents mean on S2xS2 and variance on so3, which enocdes a Gaussian distribution on SO3
        # See section 2.5 of http://ethaneade.com/lie.pdf
        # convert so3 to SO3 with exponential map
        '''
        z_mu = lie_tools.s2s2_to_SO3(z[:, :3], z[:, 3:6])
        z_log_var = z[:, 6:]
        z_std = torch.exp(.5*z_log_var) # or could do softplus z[6:]
        
        # resampling trick
        eps = torch.randn_like(z_std)
        w_eps = eps*z_std
        rot_eps = lie_tools.expmap(w_eps)
        rot_sampled = z_mu @ rot_eps
        return rot_sampled, w_eps, z_std

    def forward(self, img):
        z = self.encoder(img.view(-1,self.in_dim))
        rot, w_eps, z_std = self.reparameterize(z)

        # transform lattice by rot
        x = self.lattice @ rot # R.T*x
        y_hat = self.decoder(x)
        y_hat = y_hat.view(-1, self.ny, self.nx)
        return y_hat, w_eps, z_std

def loss_function(recon_y, y, w_eps, z_std):
    gen_loss = F.mse_loss(recon_y, y)  
    cross_entropy = torch.tensor([np.log(8*np.pi**2)], device=y.device) # cross entropy between gaussian and uniform on SO3
    entropy = lie_tools.so3_entropy(w_eps,z_std)
    kld = cross_entropy - entropy
    #assert kld > 0
    return gen_loss, kld.mean()

def eval_volume(model, nz, ny, nx, rnorm):
    '''Evaluate the model on a nz x ny x nx lattice'''
    vol_f = np.zeros((nz,ny,nx),dtype=complex)
    assert not model.training
    # evaluate the volume by zslice to avoid memory overflows
    for i, z in enumerate(np.linspace(-1,1,nz,endpoint=False)):
        x = model.lattice + torch.tensor([0,0,z], device=model.lattice.device, dtype=model.lattice.dtype)
        with torch.no_grad():
            y = model.decoder(x)
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
    nz = max(nx,ny)
    log('Loaded {} {}x{} images'.format(Nimg, ny, nx))
    particles = np.asarray([fft.ht2_center(img).astype(np.float32) for img in particles])
    assert particles.shape == (Nimg,ny,nx)
    rnorm  = [np.mean(particles), np.std(particles)]
    rnorm[0] = 0
    log('Normalizing by mean, std: {} +/- {}'.format(*rnorm))
    particles = (particles - rnorm[0])/rnorm[1]

    model = VAE(nx, ny, args.qlayers, args.qdim, args.players, args.pdim)

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
        model.lattice = model.lattice.cuda()

    optim = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.wd)
    
    # training loop
    num_epochs = args.num_epochs
    for epoch in range(num_epochs):
        loss_accum = 0
        ii = 0
        num_batches = np.ceil(Nimg / args.batch_size).astype(int)
        for minibatch_i in np.array_split(np.random.permutation(Nimg),num_batches):
            y = Variable(torch.from_numpy(np.asarray([particles[i] for i in minibatch_i])))
            if use_cuda:
                y = y.cuda()

            y_recon, w_eps, z_std = model(y) 
            gen_loss, kld = loss_function(y_recon, y, w_eps, z_std)
            loss = gen_loss + args.beta*kld/(nx*ny)
            loss.backward()
            optim.step()
            optim.zero_grad()

            loss_accum += loss.item()*len(minibatch_i)
            if ii > 0 and ii % args.log_interval == 0:
                log('# [Train Epoch: {}/{}] [{}/{} images] gen loss={:.4f}, kld={:.4f}, loss={:.4f}'.format(epoch+1, num_epochs, ii*args.batch_size, Nimg, gen_loss.item(), kld.item(), loss.item()))
            ii += 1
        log('# =====> Epoch: {} Average loss: {:.4f}'.format(epoch+1, loss_accum/Nimg))

        if args.save_intermediates:
            model.eval()
            vol, vol_f = eval_volume(model, nz, ny, nx, rnorm)
            mrc.write('{}/reconstruct.{}.mrc'.format(args.outdir,epoch), vol.astype(np.float32))
            with open('{}/reconstruct.{}.pkl'.format(args.outdir,epoch), 'wb') as f:
                pickle.dump(vol_f, f)
            with open('{}/weights.{}.pkl'.format(args.outdir,epoch), 'wb') as f:
                pickle.dump(model.state_dict(), f)
            model.train()

    ## save model weights and evaluate the model on 3D lattice
    model.eval()
    with open('{}/weights.pkl'.format(args.outdir), 'wb') as f:
        pickle.dump(model.state_dict(), f)

    vol, vol_f = eval_volume(model, nz, ny, nx, rnorm)
    mrc.write('{}/reconstruct.mrc'.format(args.outdir), vol.astype(np.float32))
    with open('{}/reconstruct.pkl'.format(args.outdir), 'wb') as f:
        pickle.dump(vol_f, f)
    
    td = dt.now()-t1
    log('Finsihed in {} ({} per epoch)'.format(td, td/num_epochs))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)

