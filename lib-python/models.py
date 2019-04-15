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
            encode_mode = 'mlp'
            ):
        super(VAE, self).__init__()
        self.nx = nx
        self.ny = ny
        self.in_dim = nx*ny
        if encode_mode == 'conv':
            self.encoder = ConvEncoder(encode_dim, encode_dim)
        elif encode_mode == 'resid':
            self.encoder = ResidLinearEncoder(nx*ny, 
                            encode_layers, 
                            encode_dim,  # hidden_dim
                            encode_dim, # out_dim
                            nn.ReLU) #in_dim -> hidden_dim
        elif encode_mode == 'mlp':
            self.encoder = MLPEncoder(nx*ny, 
                            encode_layers, 
                            encode_dim, # hidden_dim
                            encode_dim, # out_dim
                            nn.ReLU) #in_dim -> hidden_dim
        else:
            raise RuntimeError('Encoder mode {} not recognized'.format(encode_mode))
        self.latent_encoder = SO3reparameterize(encode_dim) # hidden_dim -> SO(3) latent variable
        self.decoder = self.get_decoder(decode_layers, 
                            decode_dim, 
                            nn.ReLU) #R3 -> R1
        
        # centered and scaled xy plane, values between -1 and 1
        x0, x1 = np.meshgrid(np.linspace(-1, 1, nx, endpoint=False), # FT is not symmetric around origin
                             np.linspace(-1, 1, ny, endpoint=False))
        lattice = np.stack([x0.ravel(),x1.ravel(),np.zeros(ny*nx)],1).astype(np.float32)
        self.lattice = torch.from_numpy(lattice)
    
   
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
        z_mu, z_std = self.latent_encoder(self.encoder(img))
        rot, w_eps = self.latent_encoder.sampleSO3(z_mu, z_std)

        # transform lattice by rot
        x = self.lattice @ rot # R.T*x
        y_hat = self.decoder(x)
        y_hat = y_hat.view(-1, self.ny, self.nx)
        return y_hat, z_mu, z_std, w_eps

class TiltVAE(nn.Module):
    def __init__(self, 
            nx, ny, tilt,
            encode_layers, encode_dim, 
            decode_layers, decode_dim
            ):
        super(TiltVAE, self).__init__()
        self.nx = nx
        self.ny = ny
        self.encoder = TiltEncoder(nx*ny,
                            encode_layers,
                            encode_dim, # hidden dim
                            encode_dim, # output dim
                            nn.ReLU)
        self.latent_encoder = SO3reparameterize(encode_dim) # hidden_dim -> SO(3) latent variable
        self.decoder = self.get_decoder(decode_layers, 
                            decode_dim, 
                            nn.ReLU) #R3 -> R1
        
        # centered and scaled xy plane, values between -1 and 1
        x0, x1 = np.meshgrid(np.linspace(-1, 1, nx, endpoint=False), # FT is not symmetric around origin
                             np.linspace(-1, 1, ny, endpoint=False))
        lattice = np.stack([x0.ravel(),x1.ravel(),np.zeros(ny*nx)],1).astype(np.float32)
        self.lattice = torch.tensor(lattice)
        assert tilt.shape == (3,3), 'Rotation matrix input required'
        self.tilt = torch.tensor(tilt)
    
   
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

    def forward(self, img, img_tilt):
        z_mu, z_std = self.latent_encoder(self.encoder((img, img_tilt)))
        rot, w_eps = self.latent_encoder.sampleSO3(z_mu, z_std)

        # transform lattice by rot
        x = self.lattice @ rot # R.T*x
        y_hat = self.decoder(x)
        y_hat = y_hat.view(-1, self.ny, self.nx)

        # tilt series pair
        x = self.lattice @ self.tilt @ rot
        y_hat2 = self.decoder(x)
        y_hat2 = y_hat2.view(-1, self.ny, self.nx)
        return y_hat, y_hat2, z_mu, z_std, w_eps

class TiltEncoder(nn.Module):
    def __init__(self, in_dim, nlayers, hidden_dim, out_dim, activation):
        super(TiltEncoder, self).__init__()
        assert nlayers > 2
        self.encoder1 = ResidLinearEncoder(in_dim, nlayers-2, hidden_dim, hidden_dim, activation)
        self.encoder2 = ResidLinearEncoder(hidden_dim*2, 2, hidden_dim, out_dim, activation)

    def forward(self, img):
        x, x_tilt = img
        x_enc = self.encoder1(x)
        x_tilt_enc = self.encoder1(x_tilt)
        z = self.encoder2(torch.cat((x_enc,x_tilt_enc),-1))
        return z

class ResidLinearEncoder(nn.Module):
    def __init__(self, in_dim, nlayers, hidden_dim, out_dim, activation):
        super(ResidLinearEncoder, self).__init__()
        self.in_dim = in_dim
        # define network
        layers = [nn.Linear(in_dim, hidden_dim), activation()]
        for n in range(nlayers-1):
            layers.append(ResidLinear(hidden_dim, hidden_dim))
            layers.append(activation())
        layers.append(ResidLinear(hidden_dim, out_dim))
        layers.append(activation()) # include this?
        self.main = nn.Sequential(*layers)

    def forward(self, img):
        return self.main(img.view(-1,self.in_dim))

class MLPEncoder(nn.Module):
    def __init__(self, in_dim, nlayers, hidden_dim, out_dim, activation):
        super(MLPEncoder, self).__init__()
        self.in_dim = in_dim
        # define network
        layers = [nn.Linear(in_dim, hidden_dim), activation()]
        for n in range(nlayers-1):
            layers.append(nn.Linear(hidden_dim, hidden_dim))
            layers.append(activation())
        layers.append(nn.Linear(hidden_dim, out_dim))
        layers.append(activation()) # include this?
        self.main = nn.Sequential(*layers)

    def forward(self, img):
        return self.main(img.view(-1,self.in_dim))
      
# Adapted from soumith DCGAN
class ConvEncoder(nn.Module):
    def __init__(self, hidden_dim, out_dim):
        super(ConvEncoder, self).__init__()
        ndf = hidden_dim
        self.main = nn.Sequential(
            # input is 1 x 64 x 64
            nn.Conv2d(1, ndf, 4, 2, 1, bias=False),
            nn.LeakyReLU(0.2, inplace=True),
            # state size. (ndf) x 32 x 32
            nn.Conv2d(ndf, ndf * 2, 4, 2, 1, bias=False),
            nn.BatchNorm2d(ndf * 2),
            nn.LeakyReLU(0.2, inplace=True),
            # state size. (ndf*2) x 16 x 16
            nn.Conv2d(ndf * 2, ndf * 4, 4, 2, 1, bias=False),
            nn.BatchNorm2d(ndf * 4),
            nn.LeakyReLU(0.2, inplace=True),
            # state size. (ndf*4) x 8 x 8
            nn.Conv2d(ndf * 4, ndf * 8, 4, 2, 1, bias=False),
            nn.BatchNorm2d(ndf * 8),
            nn.LeakyReLU(0.2, inplace=True),
            # state size. (ndf*8) x 4 x 4
            nn.Conv2d(ndf * 8, out_dim, 4, 1, 0, bias=False),
            nn.LeakyReLU(0.2, inplace=True), #include this?
            # state size. out_dims x 1 x 1
        )
    def forward(self, x):
        x = torch.unsqueeze(x,1)
        x = self.main(x)
        return x.view(x.size(0), -1) # flatten

class SO3reparameterize(nn.Module):
    '''Reparameterize R^N encoder output to SO(3) latent variable'''
    def __init__(self, input_dims):
        super().__init__()
        self.s2s2map = nn.Linear(input_dims, 6)
        self.so3var = nn.Linear(input_dims, 3)

        # start with big outputs
        #self.s2s2map.weight.data.uniform_(-5,5)
        #self.s2s2map.bias.data.uniform_(-5,5)

    def sampleSO3(self, z_mu, z_std):
        '''
        Reparameterize SO(3) latent variable
        # z represents mean on S2xS2 and variance on so3, which enocdes a Gaussian distribution on SO3
        # See section 2.5 of http://ethaneade.com/lie.pdf
        '''
        # resampling trick
        if not self.training:
            return z_mu, z_std
        eps = torch.randn_like(z_std)
        w_eps = eps*z_std
        rot_eps = lie_tools.expmap(w_eps)
        rot_sampled = z_mu @ rot_eps
        return rot_sampled, w_eps

    def forward(self, x):
        z = self.s2s2map(x).double()
        logvar = self.so3var(x)
        z_mu = lie_tools.s2s2_to_SO3(z[:, :3], z[:, 3:]).float()
        z_std = torch.exp(.5*logvar) # or could do softplus
        return z_mu, z_std 

        

