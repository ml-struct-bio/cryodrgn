'''Pytorch models'''

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

import lie_tools
import so3_grid

class ResidLinear(nn.Module):
    def __init__(self, nin, nout):
        super(ResidLinear, self).__init__()
        self.linear = nn.Linear(nin, nout)

    def forward(self, x):
        z = self.linear(x) + x
        return z

class HetVAE(nn.Module):
    def __init__(self, 
            nx, ny, in_dim,
            encode_layers, encode_dim, 
            decode_layers, decode_dim,
            z_dim = 1,
            encode_mode = 'mlp',
            ):
        super(HetVAE, self).__init__()
        self.nx = nx
        self.ny = ny
        self.in_dim = in_dim # nx*ny for single image or 2*nx*ny for tilt series
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
        else:
            raise RuntimeError('Encoder mode {} not recognized'.format(encode_mode))
        self.decoder = FTSliceDecoder(3+z_dim, nx, decode_layers, 
        #self.decoder = ResidLinearDecoder(3+z_dim, 1, decode_layers, 
                            decode_dim, 
                            nn.ReLU) #R3 -> R1
        
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

    def encode(self, img):
        z = self.encoder(img)
        return z[:,:self.z_dim], z[:,self.z_dim:]

    def decode(self, coords, z):
        '''coords much be lattice(s) of central slices'''
        z = z.view(z.size(0), *([1]*(coords.ndimension()-1)))
        z = torch.cat((coords,z.expand(*coords.shape[:-1],1)),dim=-1)
        y_hat = self.decoder(z)
        y_hat = y_hat.view(-1, self.ny, self.nx)
        return y_hat

    def forward(self, rot, z):
        x = self.lattice @ rot # R.T*x
        y_hat = self.decode(x,z)
        return y_hat

class BNBOpt():
    def __init__(self, model, ny, nx, tilt=None):
        super(BNBOpt, self).__init__()
        self.ny = ny
        self.nx = nx
        self.model = model
        x0, x1 = np.meshgrid(np.linspace(-1, 1, nx, endpoint=False), # FT is not symmetric around origin
                             np.linspace(-1, 1, ny, endpoint=False))
        lattice = np.stack([x0.ravel(),x1.ravel(),np.zeros(ny*nx)],1).astype(np.float32)
        self.lattice = torch.tensor(lattice)
        self.base_quat = so3_grid.base_SO3_grid()
        self.base_rot = lie_tools.quaternions_to_SO3(torch.tensor(self.base_quat))
        self.nbase = len(self.base_quat)
        assert self.nbase == 576, "Base resolution changed?"
        
        if tilt is not None:
            assert tilt.shape == (3,3)
            self.tilt = torch.tensor(tilt)

        # FT is not symmetric around origin
        assert ny == nx
        D = ny
        c = 2/(D-1)*(D/2) -1
        self.center = torch.tensor([c,c]) # pixel coordinate for img[D/2, D/2]

        L = min(int(nx/4),12)
        b, e = int(nx/2-L), int(nx/2+L)
        center_lattice = self.lattice.view(ny,nx,3)[b:e,b:e,:].contiguous().view(-1,3)
        self.L = L
        self.center_slice = (b,e)
        self.center_lattice = center_lattice

    def rotate(self, images, theta):
        '''
        images: BxYxX
        theta: Q, in radians
        '''
        images = images.expand(len(theta), *images.shape) # QxBxYxX
        cos = torch.cos(theta)
        sin = torch.sin(theta)
        rot = torch.stack([cos, sin, -sin, cos], 1).view(-1, 2, 2)
        grid = self.lattice[:,0:2] @ rot
        grid = grid.view(len(rot), self.ny, self.nx, 2) # QxYxXx2
        offset = self.center - grid[:,int(self.ny/2),int(self.nx/2)] # Qx2
        grid += offset[:,None,None,:]
        rotated = F.grid_sample(images, grid) # QxBxYxX
        return rotated.transpose(0,1) # BxQxYxX

    def eval_base_grid(self, images, images_tilt=None):
        B = images.size(0)
        b,e = self.center_slice
        def compute_err(rot, rot_center, images):
            lattice_eval = self.lattice @ rot # 48 x (Y*X) x 3
            center_lattice_eval = self.center_lattice @ rot_center # 48x11x(L*L)x3
            center_lattice_eval = center_lattice_eval.view(-1,3)
            y_hat = self.model(torch.cat((lattice_eval.view(-1,3),center_lattice_eval),0))
            y_hat_center = y_hat[-self.L**2*4*528:].view(1, 48, 11, 2*self.L, 2*self.L) # 1x 48 x 11 x L x L
            y_hat = y_hat[:-self.L**2*4*528].view(1,48,1,self.ny,self.nx) # 1x48x1xYxX
    
            theta = torch.arange(1,12,dtype=torch.float32)*2*np.pi/12 # 11 angles
            img_rot = self.rotate(images, theta) # Bx11xYxX
            images = images.unsqueeze(1) # Bx1xYxX
            images = torch.cat((images, img_rot), 1) # Bx12xYxX
            images = images.view(B,1,12,self.ny,self.nx) # Bx1x12xYxX
            err = torch.sum((images-y_hat).pow(2),(-1,-2)) # Bx48x12
            center_err_rot = torch.sum((images[:,:,1:,b:e,b:e]-y_hat[...,b:e,b:e]).pow(2),(-1,-2)) # Bx48x11
            center_err = torch.sum((images[:,:,0:1,b:e,b:e]-y_hat_center).pow(2),(-1,-2)) # Bx48x11
            err[:,:,1:] += (center_err-center_err_rot)
            err = err.view(B,self.nbase) # Bx576
            return err
        base_rot = self.base_rot.view(48,12,3,3)
        rot = base_rot[:,0,...] 
        rot_center = base_rot[:,1:,...]
        err = compute_err(rot, rot_center, images)
        if images_tilt is not None:
            err_tilt = compute_err(self.tilt @ rot, self.tilt @ rot_center, images_tilt)
            err += err_tilt
        mini = torch.argmin(err,1)
        return mini.cpu().numpy()

    def eval_grid(self, images, rot, NQ, images_tilt=None):
        '''
        images: B x NY x NX 
        rot: (NxQ) x 3 x 3 rotation matrics (N=1 for base grid, N=B for incremental grid)
        NQ: number of slices evaluated for each image
        '''
        B = images.size(0)
        y_hat = self.model(self.lattice @ rot)
        y_hat = y_hat.view(-1,NQ,self.ny,self.nx) #1xQxYxX for base grid, Bx8xYxX for incremental grid
        images = images.view(B,1,self.ny,self.nx) # Bx1xYxX
        err = torch.sum((images-y_hat).pow(2),(-1,-2)) # BxQ

        if images_tilt is not None:
            y_hat = self.model(self.lattice @ self.tilt @ rot)
            y_hat = y_hat.view(-1,NQ,self.ny,self.nx) #1xQxYxX for base grid, Bx8xYxX for incremental grid
            images_tilt = images_tilt.view(B,1,self.ny,self.nx) # Bx1xYxX
            err_tilt = torch.sum((images_tilt-y_hat).pow(2),(-1,-2)) # BxQ
            mini = torch.argmin(err+err_tilt,1) # B
        else:
            mini = torch.argmin(err,1) # B
        return mini.cpu().numpy()

    def opt_theta(self, images, images_tilt=None, niter=5):
        B = images.size(0)
        assert not self.model.training
        with torch.no_grad():
            #min_i2 = self.eval_grid(images, self.base_rot, self.nbase,
            #                          images_tilt=images_tilt)
        
            min_i = self.eval_base_grid(images, images_tilt=images_tilt) # 576 slices
            min_quat = self.base_quat[min_i]
            s2i, s1i = so3_grid.get_base_indr(min_i)
            for iter_ in range(1,niter+1):
                neighbors = [so3_grid.get_neighbor(min_quat[i], s2i[i], s1i[i], iter_) for i in range(B)]
                quat = np.array([x[0] for x in neighbors])
                ind = np.array([x[1] for x in neighbors])
                rot = lie_tools.quaternions_to_SO3(torch.tensor(quat))
                min_i = self.eval_grid(images, rot, 8, images_tilt=images_tilt)
                min_ind = ind[np.arange(B), min_i]
                s2i, s1i = min_ind.T
                min_quat = quat[np.arange(B),min_i]
        return lie_tools.quaternions_to_SO3(torch.tensor(min_quat))

class BNBHetOpt():
    def __init__(self, model, ny, nx, tilt=None):
        super(BNBHetOpt, self).__init__()
        self.ny = ny
        self.nx = nx
        self.model = model # this is the VAE module
        self.base_quat = so3_grid.base_SO3_grid()
        self.nbase = len(self.base_quat)
        self.base_rot = lie_tools.quaternions_to_SO3(torch.tensor(self.base_quat))

        if tilt is not None:
            assert tilt.shape == (3,3)
            self.tilt = torch.tensor(tilt)

        # FT is not symmetric around origin
        assert ny == nx
        D = ny
        c = 2/(D-1)*(D/2) -1
        self.center = torch.tensor([c,c]) # pixel coordinate for img[D/2, D/2]

    def rotate(self, images, theta):
        '''
        images: BxYxX
        theta: Q, in radians
        '''
        images = images.expand(len(theta), *images.shape) # QxBxYxX
        cos = torch.cos(theta)
        sin = torch.sin(theta)
        rot = torch.stack([cos, sin, -sin, cos], 1).view(-1, 2, 2)
        grid = self.model.lattice[:,0:2] @ rot
        grid = grid.view(len(rot), self.ny, self.nx, 2) # QxYxXx2
        offset = self.center - grid[:,int(self.ny/2),int(self.nx/2)] # Qx2
        grid += offset[:,None,None,:]
        rotated = F.grid_sample(images, grid) # QxBxYxX
        return rotated.transpose(0,1) # BxQxYxX

    def eval_base_grid(self, images, z, images_tilt=None):
        B = images.size(0)
        def compute_err(rot, images):
            lattice_eval = (self.model.lattice @ rot) # 48 x (NxY) x 3
            lattice_eval = lattice_eval.expand(B,*lattice_eval.shape) # B x 48 x (XxY) x 3
            y_hat = self.model.decode(lattice_eval, z) # (Bx48)xYxX
            y_hat = y_hat.view(B,1,len(rot),self.ny,self.nx) # Bx1x48xYxX

            theta = torch.arange(1,12,dtype=torch.float32)*2*np.pi/12 # 11 angles
            img_rot = self.rotate(images, theta) # Bx11xYxX
            images = images.unsqueeze(1) # Bx1xYxX
            images = torch.cat((images, img_rot), 1) # Bx12xYxX
            images = images.view(B,12,1,self.ny,self.nx) # Bx12x1xYxX

            err = torch.sum((images-y_hat).pow(2),(-1,-2)) # Bx12x48
            err = err.transpose(1,2).contiguous() # Bx48x12, since psi's are the fast dimension
            err = err.view(B,self.nbase) # Bx576
            return err
        rot = self.base_rot[::12] # 576/12 = 48 images
        err = compute_err(rot, images)
        if images_tilt is not None:
            err_tilt = compute_err(self.tilt @ rot, images_tilt)
            err += err_tilt
        mini = torch.argmin(err,1)
        return mini.cpu().numpy()

    def eval_grid(self, images, rot, z, NQ, images_tilt=None):
        B = z.size(0)
        images = images.view(B,1,self.ny,self.nx) # Bx1xYxX
        y_hat = self.model.decode(self.model.lattice @ rot, z)
        y_hat = y_hat.view(B, NQ, self.ny, self.nx) # don't think we need to do this
        err = torch.sum((images-y_hat).pow(2),(-1,-2)) # B x Q

        if images_tilt is not None:
            images_tilt = images_tilt.view(B,1,self.ny,self.nx) # Bx1xYxX
            y_hat = self.model.decode(self.model.lattice @ self.tilt @ rot, z)
            y_hat = y_hat.view(B, NQ, self.ny, self.nx) # don't think we need to do this
            err_tilt = torch.sum((images_tilt-y_hat).pow(2),(-1,-2)) # B x Q
            mini = torch.argmin(err+err_tilt, 1)
        else:
            mini = torch.argmin(err,1)
        return mini.cpu().numpy()
        
    def opt_theta(self, images, z, images_tilt=None, niter=5):
        B = images.size(0)
        assert not self.model.training
        with torch.no_grad():
            # expand the base grid B times since each image has a different z
            #base_rot = self.base_rot.expand(B,*self.base_rot.shape) # B x 576 x 3 x 3
            #min_i = self.eval_grid(images, base_rot, z, self.nbase,
            #                       images_tilt=images_tilt) # B x 576 slices
            min_i = self.eval_base_grid(images, z, images_tilt)
            min_quat = self.base_quat[min_i]
            s2i, s1i = so3_grid.get_base_indr(min_i)
            for iter_ in range(1,niter+1):
                neighbors = [so3_grid.get_neighbor(min_quat[i], s2i[i], s1i[i], iter_) for i in range(B)]
                quat = np.array([x[0] for x in neighbors])
                ind = np.array([x[1] for x in neighbors])
                rot = lie_tools.quaternions_to_SO3(torch.tensor(quat))
                min_i = self.eval_grid(images, rot, z, 8, images_tilt=images_tilt)
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

    def forward(self, lattice):
        '''Call forward on central slices only'''
        #assert lattice.shape[-2:] == (self.D**2,3)
        #assert torch.nonzero(lattice[...,self.center,:]).size(0) == 0
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

        

