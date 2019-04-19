'''Pytorch models'''

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

import lie_tools
import so3_grid
import utils

log = utils.log

class HetOnlyVAE(nn.Module):
    def __init__(self, lattice, # Lattice object
            in_dim, # nx*ny for single image
            encode_layers, encode_dim, 
            decode_layers, decode_dim,
            z_dim = 1, 
            encode_mode = 'resid'):
        super(HetOnlyVAE, self).__init__()
        self.lattice = lattice
        self.in_dim = in_dim 
        self.z_dim = z_dim
        if encode_mode == 'conv':
            assert lattice.D == 64, "CNN encoder is hard-coded for 64x64 images"
            self.encoder = ConvEncoder(encode_dim, z_dim*2)
        elif encode_mode == 'resid':
            self.encoder = ResidLinearMLP(in_dim, 
                            encode_layers, # nlayers
                            encode_dim,  # hidden_dim
                            z_dim*2, # out_dim
                            nn.ReLU) 
        elif encode_mode == 'mlp':
            self.encoder = MLP(in_dim, 
                            encode_layers, 
                            encode_dim, # hidden_dim
                            z_dim*2, # out_dim
                            nn.ReLU) #in_dim -> hidden_dim
        elif encode_mode == 'tilt':
            self.encoder = TiltEncoder(in_dim,
                            encode_layers,
                            encode_dim,
                            z_dim*2,
                            nn.ReLU)
        else:
            raise RuntimeError('Encoder mode {} not recognized'.format(encode_mode))
        self.decoder = FTSliceDecoder(3+z_dim, # input dim
                            lattice.D, # lattice size
                            decode_layers, # nlayers
                            decode_dim, # hidden dim
                            nn.ReLU) #R3 -> R1
   
    def reparameterize(self, mu, logvar):
        if not self.training:
            return mu
        std = torch.exp(.5*logvar)
        eps = torch.randn_like(std)
        return eps*std + mu

    def encode(self, img):
        z = self.encoder(img)
        return z[:,:self.z_dim], z[:,self.z_dim:]

    def cat_z(self, coords, z):
        assert coords.shape[-1] == 3
        z = z.view(z.size(0), *([1]*(coords.ndimension()-1)))
        z = torch.cat((coords,z.expand(*coords.shape[:-1],1)),dim=-1)
        return z

    def decode(self, rot, z):
        '''
        rot: Bx3x3 rotation matrices
        z: Bxzdim latent coordinate
        '''
        x = self.lattice.coords @ rot # R.T*x
        y_hat = self.decoder(self.cat_z(x,z))
        return y_hat

class FTSliceDecoder(nn.Module):
    '''
    Evaluate a central slice out of a 3D FT of a model, returns representation in
    Hartley reciprocal space

    Exploits the symmetry of the FT where F*(x,y) = F(-x,-y) and only
    evaluates half of the lattice. The decoder is f(x,y,z) => real, imag
    '''
    def __init__(self, in_dim, D, nlayers, hidden_dim, activation):
        '''D: image width or height'''
        super(FTSliceDecoder, self).__init__()
        self.decoder = ResidLinearMLP(in_dim, nlayers, hidden_dim, 2, activation)
        D2 = int(D/2)

        ### various pixel indices to keep track of 
        self.center = D2*D + D2 
        self.extra = np.arange((D2+1)*D, D**2, D) # bottom-left column without conjugate pair
        # evalute the top half of the image up through the center pixel 
        # and extra bottom-left column (todo: just evaluate a D-1 x D-1 image so 
        # we don't have to worry about this)
        self.all_eval = np.concatenate((np.arange(self.center+1), self.extra))
        
        # pixel indices for the top half of the image up to (but not incl) 
        # the center pixel and excluding the top row and left-most column
        i, j = np.meshgrid(np.arange(1,D),np.arange(1,D2+1))
        self.top = (j*D+i).ravel()[:-D2]

        # pixel indices for bottom half of the image after the center pixel
        # excluding left-most column and given in reverse order
        i, j =np.meshgrid(np.arange(1,D),np.arange(D2,D))
        self.bottom_rev = (j*D+i).ravel()[D2:][::-1].copy()

        self.D = D
        self.D2 = D2

    def forward_symmetric(self, lattice, c):
        '''
        central slices with a symmetrizing mask

        lattice: -1 x (2*c+1) x 3+zdim
        c: index of center pixel
        '''
        image = torch.empty((*lattice.shape[:-1],2)) 
        top_half = self.decode(lattice[...,0:c+1,:])
        image[..., 0:c+1,:] = top_half 
        # the bottom half of the image is the complex conjugate of the top half
        image[...,c+1:,:] = top_half[...,np.arange(c-1,-1,-1)]*torch.tensor([1.,-1.])
        return image

    def forward(self, lattice):
        '''Call forward on DxD central slices only'''
        image = torch.empty((*lattice.shape[:-1],2))
        top_half = self.decode(lattice[...,self.all_eval,:])
        image[..., self.all_eval, :] = top_half
        # the bottom half of the image is the complex conjugate of the top half
        image[...,self.bottom_rev, :] = top_half[...,self.top,:]*torch.tensor([1.,-1.])
        return image

    def decode(self, lattice):
        '''Return FT transform'''
        # convention: only evalute the -z points
        w = lattice[...,2] > 0.0
        lattice[...,0:3][w] = -lattice[...,0:3][w] # negate lattice coordinates where z > 0
        result = self.decoder(lattice)
        result[...,1][w] *= -1 # replace with complex conjugate to get correct values for original lattice positions
        return result

    def translate(self, coords, img, t):
        '''
        Translate an image by phase shifting its Fourier transform
        
        coords: wavevectors between [-.5,.5] (img_dims x 2)
        img: FT of image (B x img_dims x 2)
        t: shift in pixels (B x 2)

        img_dims can either be 2D or 1D (unraveled image) 
        '''
        t = t.unsqueeze(2) # Bx2x1 to be able to do bmm
        tfilt = coords @ -t * -2 * np.pi
        tfilt = tfilt.squeeze()
        c = torch.cos(tfilt) # BxN
        s = torch.sin(tfilt) # BxN
        return torch.stack([img[...,0]*c-img[...,1]*s,img[...,0]*s+img[...,1]*c],-1)

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
        assert encode_layers > 2
        if encode_mode == 'conv':
            self.encoder = ConvEncoder(encode_dim, encode_dim)
        elif encode_mode == 'resid':
            self.encoder = ResidLinearMLP(nx*ny, 
                            encode_layers-2, 
                            encode_dim,  # hidden_dim
                            encode_dim, # out_dim
                            nn.ReLU) #in_dim -> hidden_dim
        elif encode_mode == 'mlp':
            self.encoder = MLP(nx*ny, 
                            encode_layers-2, 
                            encode_dim, # hidden_dim
                            encode_dim, # out_dim
                            nn.ReLU) #in_dim -> hidden_dim
        else:
            raise RuntimeError('Encoder mode {} not recognized'.format(encode_mode))
        self.so3_encoder = SO3reparameterize(encode_dim, 1) # hidden_dim -> SO(3) latent variable
        self.trans_encoder = ResidLinearMLP(encode_dim, 1, encode_dim, 4, nn.ReLU)
        self.decoder = FTSliceDecoder(3, nx, decode_layers, decode_dim, nn.ReLU)
        
        # centered and scaled xy plane, values between -1 and 1
        x0, x1 = np.meshgrid(np.linspace(-1, 1, nx, endpoint=False), # FT is not symmetric around origin
                             np.linspace(-1, 1, ny, endpoint=False))
        lattice = np.stack([x0.ravel(),x1.ravel(),np.zeros(ny*nx)],1).astype(np.float32)
        self.lattice = torch.tensor(lattice)

    def reparameterize(self, mu, logvar):
        if not self.training:
            return mu
        std = torch.exp(.5*logvar)
        eps = torch.randn_like(std)
        return eps*std + mu

    def forward(self, img):
        img = img[...,0] - img[...,1]
        enc = nn.ReLU()(self.encoder(img.view(-1,self.in_dim)))
        z_mu, z_std = self.so3_encoder(enc)
        rot, w_eps = self.so3_encoder.sampleSO3(z_mu, z_std)
        z = self.trans_encoder(enc)
        tmu, tlogvar = z[:,:2], z[:,2:]
        t = self.reparameterize(tmu, tlogvar)
        # transform lattice by rot
        x = self.lattice @ rot # R.T*x
        y_hat = self.decoder(x)
        # translate image
        y_hat = self.decoder.translate(self.lattice[:,0:2]/2, y_hat, t)
        y_hat = y_hat.view(-1, self.ny, self.nx, 2)
        return y_hat, z_mu, z_std, w_eps, tmu, tlogvar

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
        #self.decoder = ResidLinearMLP(3, decode_layers, decode_dim, 2, nn.ReLU)
        self.decoder = FTSliceDecoder(3, nx, decode_layers, decode_dim, nn.ReLU)
        
        # centered and scaled xy plane, values between -1 and 1
        x0, x1 = np.meshgrid(np.linspace(-1, 1, nx, endpoint=False), # FT is not symmetric around origin
                             np.linspace(-1, 1, ny, endpoint=False))
        lattice = np.stack([x0.ravel(),x1.ravel(),np.zeros(ny*nx)],1).astype(np.float32)
        self.lattice = torch.tensor(lattice)
        assert tilt.shape == (3,3), 'Rotation matrix input required'
        self.tilt = torch.tensor(tilt)

    def forward(self, img, img_tilt):
        img = img[...,0] - img[...,1]
        img_tilt = img_tilt[...,0] - img_tilt[...,1]
        z_mu, z_std = self.latent_encoder(nn.ReLU()(self.encoder((img, img_tilt))))
        rot, w_eps = self.latent_encoder.sampleSO3(z_mu, z_std)

        # transform lattice by rot
        x = self.lattice @ rot # R.T*x
        y_hat = self.decoder(x)
        y_hat = y_hat.view(-1, self.ny, self.nx, 2)

        # tilt series pair
        x = self.lattice @ self.tilt @ rot
        y_hat2 = self.decoder(x)
        y_hat2 = y_hat2.view(-1, self.ny, self.nx, 2)
        return y_hat, y_hat2, z_mu, z_std, w_eps

class TiltEncoder(nn.Module):
    def __init__(self, in_dim, nlayers, hidden_dim, out_dim, activation):
        super(TiltEncoder, self).__init__()
        assert nlayers > 2
        self.encoder1 = ResidLinearMLP(in_dim, nlayers-2, hidden_dim, hidden_dim, activation)
        self.encoder2 = ResidLinearMLP(hidden_dim*2, 2, hidden_dim, out_dim, activation)
        self.in_dim = in_dim

    def forward(self, img):
        x, x_tilt = img
        x_enc = self.encoder1(x.view(-1,self.in_dim))
        x_tilt_enc = self.encoder1(x_tilt.view(-1,self.in_dim))
        z = self.encoder2(torch.cat((x_enc,x_tilt_enc),-1))
        return z

class ResidLinearMLP(nn.Module):
    def __init__(self, in_dim, nlayers, hidden_dim, out_dim, activation):
        super(ResidLinearMLP, self).__init__()
        layers = [ResidLinear(in_dim, hidden_dim) if in_dim == hidden_dim else nn.Linear(in_dim, hidden_dim), activation()]
        for n in range(nlayers):
            layers.append(ResidLinear(hidden_dim, hidden_dim))
            layers.append(activation())
        layers.append(ResidLinear(hidden_dim, out_dim) if out_dim == hidden_dim else nn.Linear(hidden_dim, out_dim))
        self.main = nn.Sequential(*layers)

    def forward(self, x):
        return self.main(x)

class ResidLinear(nn.Module):
    def __init__(self, nin, nout):
        super(ResidLinear, self).__init__()
        self.linear = nn.Linear(nin, nout)

    def forward(self, x):
        z = self.linear(x) + x
        return z

class MLP(nn.Module):
    def __init__(self, in_dim, nlayers, hidden_dim, out_dim, activation):
        super(MLP, self).__init__()
        layers = [nn.Linear(in_dim, hidden_dim), activation()]
        for n in range(nlayers):
            layers.append(nn.Linear(hidden_dim, hidden_dim))
            layers.append(activation())
        layers.append(nn.Linear(hidden_dim, out_dim))
        self.main = nn.Sequential(*layers)

    def forward(self, x):
        return self.main(x)
      
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
            # state size. out_dims x 1 x 1
        )
    def forward(self, x):
        x = x.view(-1,1,64,64)
        x = self.main(x)
        return x.view(x.size(0), -1) # flatten

class SO3reparameterize(nn.Module):
    '''Reparameterize R^N encoder output to SO(3) latent variable'''
    def __init__(self, input_dims, nlayers=None):
        super().__init__()
        if nlayers is not None:
            self.main = ResidLinearMLP(input_dims, nlayers, input_dims, 9, nn.ReLU)
        else:
            self.main = nn.Linear(input_dims, 9)

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
        z = self.main(x)
        z1 = z[:,:3].double()
        z2 = z[:,3:6].double()
        logvar = z[:,6:]
        z_mu = lie_tools.s2s2_to_SO3(z1,z2).float()
        z_std = torch.exp(.5*logvar) # or could do softplus
        return z_mu, z_std 

        

