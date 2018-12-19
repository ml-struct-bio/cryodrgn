'''Pytorch models'''

import numpy as np
import torch
import torch.nn as nn

import lie_tools

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


