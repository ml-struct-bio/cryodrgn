'''Pytorch models'''

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

import lie_tools

class ResidLinear(nn.Module):
    def __init__(self, nin, nout):
        super(ResidLinear, self).__init__()
        self.linear = nn.Linear(nin, nout)

    def forward(self, x):
        z = self.linear(x) + x
        return z

class VAE(nn.Module):
    def __init__(self, 
            nx, ny, 
            encode_layers, encode_dim, 
            decode_layers, decode_dim,
            group_reparam_in_dims=128
            ):
        super(VAE, self).__init__()
        self.nx = nx
        self.ny = ny
        self.in_dim = nx*ny
        self.encoder = self.get_encoder(nx*ny, 
                            encode_layers, 
                            encode_dim, 
                            nn.ReLU) #in_dim -> hidden_dim
        self.latent_encoder = SO3reparameterize(encode_dim, 
                                group_reparam_in_dims, 
                                nn.ReLU) # hidden_dim -> SO(3) latent variable
        self.decoder = self.get_decoder(decode_layers, 
                            decode_dim, 
                            nn.ReLU) #R3 -> R1
        
        # centered and scaled xy plane, values between -1 and 1
        x0, x1 = np.meshgrid(np.linspace(-1, 1, nx, endpoint=False), # FT is not symmetric around origin
                             np.linspace(-1, 1, ny, endpoint=False))
        lattice = np.stack([x0.ravel(),x1.ravel(),np.zeros(ny*nx)],1).astype(np.float32)
        self.lattice = torch.from_numpy(lattice)
    
    def get_encoder(self, in_dim, nlayers, hidden_dim, activation):
        '''Return an image encoder NN'''
        layers = [nn.Linear(in_dim, hidden_dim), activation()]
        for n in range(nlayers):
            layers.append(ResidLinear(hidden_dim, hidden_dim))
            layers.append(activation())
        return nn.Sequential(*layers)
    
    def get_decoder(self, nlayers, hidden_dim, activation):
        '''
        Return a NN mapping R3 cartesian coordinates to R1 electron density
        (represented in Hartley reciprocal space)
        '''
        layers = [nn.Linear(3, hidden_dim), activation()]
        for n in range(nlayers):
            layers.append(ResidLinear(hidden_dim, hidden_dim))
            layers.append(activation())
        layers.append(nn.Linear(hidden_dim,1))
        return nn.Sequential(*layers)

    def forward(self, img):
        z = self.encoder(img.view(-1,self.in_dim))
        rot, w_eps, z_std = self.latent_encoder(z)

        # transform lattice by rot
        x = self.lattice @ rot # R.T*x
        y_hat = self.decoder(x)
        y_hat = y_hat.view(-1, self.ny, self.nx)
        return y_hat, w_eps, z_std


class SO3reparameterize(nn.Module):
    '''Reparameterize R^N encoder output to SO(3) latent variable'''
    def __init__(self, input_dims, group_reparam_dims, activation):
        super().__init__()
        self.linear = nn.Linear(input_dims, group_reparam_dims)
        self.activation = activation()
        self.s2s2map = nn.Linear(group_reparam_dims, 6)
        self.so3var = nn.Linear(group_reparam_dims, 3)

        # start with big outputs
        self.s2s2map.weight.data.uniform_(-5,5)
        self.s2s2map.bias.data.uniform_(-5,5)

    def reparameterize(self, z, logvar):
        '''
        Reparameterize SO(3) latent variable
        # z represents mean on S2xS2 and variance on so3, which enocdes a Gaussian distribution on SO3
        # See section 2.5 of http://ethaneade.com/lie.pdf
        '''
        z_mu = lie_tools.s2s2_to_SO3(z[:, :3], z[:, 3:]).float()
        z_std = torch.exp(.5*logvar) # or could do softplus
        
        # resampling trick
        eps = torch.randn_like(z_std)
        w_eps = eps*z_std
        rot_eps = lie_tools.expmap(w_eps)
        rot_sampled = z_mu @ rot_eps
        return rot_sampled, w_eps, z_std

    def forward(self, x):
        x = self.activation(self.linear(x))
        z = self.s2s2map(x).double()
        logvar = self.so3var(x)
        return self.reparameterize(z, logvar)

        

