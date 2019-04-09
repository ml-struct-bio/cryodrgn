'''Pytorch models'''

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

import lie_tools
import so3_grid
import utils

log = utils.log

class ResidLinear(nn.Module):
    def __init__(self, nin, nout):
        super(ResidLinear, self).__init__()
        self.linear = nn.Linear(nin, nout)

    def forward(self, x):
        z = self.linear(x) + x
        return z

class Lattice:
    def __init__(self, D):
        # centered and scaled xy plane, values between -1 and 1
        # endpoint=False since FT is not symmetric around origin
        x0, x1 = np.meshgrid(np.linspace(-1, 1, D, endpoint=False), 
                             np.linspace(-1, 1, D, endpoint=False))
        coords = np.stack([x0.ravel(),x1.ravel(),np.zeros(D**2)],1).astype(np.float32)
        self.coords = torch.tensor(coords)
        self.D = D
        self.D2 = int(D/2)
        
        self.square_mask = {}
        self.circle_mask = {}

    def get_square_lattice(self, L):
        b,e = self.D2-L, self.D2+L+1
        center_lattice = self.coords.view(self.D,self.D,3)[b:e,b:e,:].contiguous().view(-1,3)
        return center_lattice

    def get_square_mask(self, L):
        '''Return a binary mask for self.coords which restricts coordinates to a centered square lattice'''
        if L in self.square_mask:
            return self.square_mask[L]
        assert 2*L+1 < self.D, 'Mask with size {} too large for lattice with size {}'.format(L,D)
        log('Using square lattice of size {}x{}'.format(2*L+1,2*L+1))
        b,e = self.D2-L, self.D2+L
        c1 = self.coords.view(self.D,self.D,3)[b,b]
        c2 = self.coords.view(self.D,self.D,3)[e,e]
        m1 = self.coords[:,0] >= c1[0]
        m2 = self.coords[:,0] <= c2[0]
        m3 = self.coords[:,1] >= c1[1]
        m4 = self.coords[:,1] <= c2[1]
        mask = m1*m2*m3*m4
        self.square_mask[L] = mask
        return mask

    def get_circular_mask(self, R):
        '''Return a binary mask for self.coords which restricts coordinates to a centered circular lattice'''
        if R in self.circle_mask:
            return self.circle_mask[R]
        assert 2*R+1 < self.D, 'Mask with radius {} too large for lattice with size {}'.format(R,D)
        log('Using circular lattice with radius {}'.format(R))
        r = 2*R/self.D
        mask = self.coords.pow(2).sum(-1) < r**2
        self.circle_mask[R] = mask
        return mask
        
class HetVAE(nn.Module):
    def __init__(self, lattice, # Lattice object
            in_dim, # nx*ny for single image or 2*nx*ny for tilt series
            encode_layers, encode_dim, 
            decode_layers, decode_dim,
            z_dim = 1,
            encode_mode = 'mlp',
            ):
        super(HetVAE, self).__init__()
        self.lattice = lattice
        self.in_dim = in_dim 
        self.z_dim = z_dim
        if encode_mode == 'conv':
            self.encoder = ConvEncoder(encode_dim, z_dim*2)
        elif encode_mode == 'resid':
            self.encoder = ResidLinearEncoder(in_dim, 
                            encode_layers, 
                            encode_dim,  # hidden_dim
                            z_dim*2, # out_dim
                            nn.ReLU) #in_dim -> hidden_dim
        elif encode_mode == 'mlp':
            self.encoder = MLPEncoder(in_dim, 
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
        #self.decoder = ResidLinearDecoder(3+z_dim, 1, decode_layers, 
        self.decoder = FTSliceDecoder(3+z_dim, lattice.D, decode_layers, 
                            decode_dim, 
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

    def forward(self, rot, z):
        '''
        rot: Bx3x3 rotation matrices
        z: Bxzdim latent coordinate
        '''
        x = self.lattice.coords @ rot # R.T*x
        y_hat = self.decoder(self.cat_z(x,z))
        return y_hat

class BNNBHomo:
    '''Branch and no bound for homogeneous reconstruction'''
    def __init__(self, decoder, lattice, tilt=None):
        self.decoder = decoder
        self.lattice = lattice
        self.base_quat = so3_grid.base_SO3_grid()
        self.base_rot = lie_tools.quaternions_to_SO3(torch.tensor(self.base_quat))
        self.nbase = len(self.base_quat)
        assert self.nbase == 576, "Base resolution changed?"
        self.tilt = tilt

    def mask_image(self, images, mask):
        return images[:,mask]

    def eval_grid(self, images, rot, NQ, L, images_tilt=None):
        '''
        images: B x NY x NX 
        rot: (NxQ) x 3 x 3 rotation matrics (N=1 for base grid, N=B for incremental grid)
        NQ: number of slices evaluated for each image
        L: radius of fourier components to evaluate
        '''
        B = images.size(0)
        mask = self.lattice.get_square_mask(L)
        coords = self.lattice.coords[mask]
        c = int(coords.size(-2)/2)
        YX = coords.size(-2)
        assert YX == (2*L+1)**2
        def compute_err(images, rot):
            images = images.view(B,-1)[:,mask]
            y_hat = self.decoder.forward_symmetric(coords @ rot, c)
            y_hat = y_hat.view(-1,NQ,YX) #1xQxYX for base grid, Bx8xYX for incremental grid
            images = images.unsqueeze(1) # Bx1xYX
            err = torch.sum((images-y_hat).pow(2),-1) # BxQ
            return err
        err = compute_err(images, rot)
        if images_tilt is not None:
            err_tilt = compute_err(images_tilt, self.tilt @ rot)
            err += err_tilt
        mini = torch.argmin(err,1) # B
        return mini.cpu().numpy()

    def opt_theta(self, images, L, images_tilt=None, niter=5):
        B = images.size(0)
        assert not self.decoder.training
        with torch.no_grad():
            min_i = self.eval_grid(images, self.base_rot, self.nbase, L, images_tilt=images_tilt)
            min_quat = self.base_quat[min_i]
            s2i, s1i = so3_grid.get_base_indr(min_i)
            for iter_ in range(1,niter+1):
                neighbors = [so3_grid.get_neighbor(min_quat[i], s2i[i], s1i[i], iter_) for i in range(B)]
                quat = np.array([x[0] for x in neighbors])
                ind = np.array([x[1] for x in neighbors])
                rot = lie_tools.quaternions_to_SO3(torch.tensor(quat))
                min_i = self.eval_grid(images, rot, 8, L, images_tilt=images_tilt)
                min_ind = ind[np.arange(B), min_i]
                s2i, s1i = min_ind.T
                min_quat = quat[np.arange(B),min_i]
        return lie_tools.quaternions_to_SO3(torch.tensor(min_quat))

class BNNBHet:
    def __init__(self, model, lattice, tilt=None):
        self.model = model
        self.lattice = lattice
        self.base_quat = so3_grid.base_SO3_grid()
        self.nbase = len(self.base_quat)
        self.base_rot = lie_tools.quaternions_to_SO3(torch.tensor(self.base_quat))
        self.tilt = tilt

    def eval_grid(self, images, rot, z, NQ, L, images_tilt=None):
        B = z.size(0)
        mask = self.lattice.get_square_mask(L)
        coords = self.lattice.coords[mask]
        c = int(coords.size(-2)/2)
        YX = coords.size(-2)
        def compute_err(images, rot):
            images = images.view(B,-1)[:,mask]
            x = self.model.cat_z(coords @ rot, z)
            y_hat = self.model.decoder.forward_symmetric(x,c)
            y_hat = y_hat.view(B, NQ, YX) # BxQxYX
            images = images.unsqueeze(1) # Bx1xYX
            err = torch.sum((images-y_hat).pow(2),-1) # BxQ
            return err
        err = compute_err(images,rot)
        if images_tilt is not None:
            err_tilt = compute_err(images_tilt, self.tilt @ rot)
            err += err_tilt
        mini = torch.argmin(err,1)
        return mini.cpu().numpy()
        
    def opt_theta(self, images, z, images_tilt=None, niter=5, L=None):
        B = images.size(0)
        assert not self.model.training
        with torch.no_grad():
            # expand the base grid B times since each image has a different z
            base_rot = self.base_rot.expand(B,*self.base_rot.shape) # B x 576 x 3 x 3
            min_i = self.eval_grid(images, base_rot, z, self.nbase, L, images_tilt=images_tilt) # B x 576 slices
            min_quat = self.base_quat[min_i]
            s2i, s1i = so3_grid.get_base_indr(min_i)
            for iter_ in range(1,niter+1):
                neighbors = [so3_grid.get_neighbor(min_quat[i], s2i[i], s1i[i], iter_) for i in range(B)]
                quat = np.array([x[0] for x in neighbors])
                ind = np.array([x[1] for x in neighbors])
                rot = lie_tools.quaternions_to_SO3(torch.tensor(quat))
                min_i = self.eval_grid(images, rot, z, 8, L, images_tilt=images_tilt)
                min_ind = ind[np.arange(B), min_i]
                s2i, s1i = min_ind.T
                min_quat = quat[np.arange(B),min_i]
        return lie_tools.quaternions_to_SO3(torch.tensor(min_quat))

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
        self.decoder = ResidLinearDecoder(in_dim, 2, nlayers, hidden_dim, activation)
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
        image = torch.empty(lattice.shape[:-1]) 
        top_half = self.decode(lattice[...,0:c+1,:])
        image[..., 0:c+1] = top_half[...,0] - top_half[...,1] # hartley transform
        # the bottom half of the image is the complex conjugate of the top half
        image[...,c+1:] = (top_half[...,0] + top_half[...,1])[...,np.arange(c-1,-1,-1)]
        return image

    def forward(self, lattice):
        '''Call forward on DxD central slices only'''
        image = torch.empty(lattice.shape[:-1])
        top_half = self.decode(lattice[...,self.all_eval,:])
        image[..., self.all_eval] = top_half[...,0] - top_half[...,1] # hartley transform
        # the bottom half of the image is the complex conjugate of the top half
        image[...,self.bottom_rev] = top_half[...,self.top,0] + top_half[...,self.top,1]
        return image

    def decode(self, lattice):
        '''Return FT transform'''
        # convention: only evalute the -z points
        w = lattice[...,2] > 0.0
        lattice[...,0:3][w] = -lattice[...,0:3][w] # negate lattice coordinates where z > 0
        result = self.decoder(lattice)
        result[...,1][w] *= -1 # replace with complex conjugate to get correct values for original lattice positions
        return result

class ResidLinearDecoder(nn.Module):
    '''
    A NN mapping R3 cartesian coordinates to R1 electron density
    (represented in Hartley reciprocal space)
    '''
    def __init__(self, in_dim, out_dim, nlayers, hidden_dim, activation):
        super(ResidLinearDecoder, self).__init__()
        layers = [nn.Linear(in_dim, hidden_dim), activation()]
        for n in range(nlayers):
            layers.append(ResidLinear(hidden_dim, hidden_dim))
            layers.append(activation())
        layers.append(nn.Linear(hidden_dim,out_dim))
        self.main = nn.Sequential(*layers)

    def forward(self, x):
        return self.main(x)

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
        layers.append(nn.Linear(hidden_dim, out_dim))
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
            # state size. out_dims x 1 x 1
        )
    def forward(self, x):
        x = torch.unsqueeze(x,1)
        x = self.main(x)
        return x.view(x.size(0), -1) # flatten

### not used in this branch ###

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

        

