'''Branch and bound classes'''

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

import lie_tools
import so3_grid
import utils

log = utils.log

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
        mask = self.lattice.get_circular_mask(L)
        coords = self.lattice.coords[mask]
        c = int(coords.size(-2)/2)
        YX = coords.size(-2)
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
        mask = self.lattice.get_circular_mask(L)
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
        
    def eval_base_grid(self, images, rotated_images, z, L):
        '''
        images: BxYxX
        rotated_images:Bx11xYxX
        '''
        B = z.size(0)
        mask = self.lattice.get_circular_mask(L)
        coords = self.lattice.coords[mask]
        c = int(coords.size(-2)/2)
        YX = coords.size(-2)

        # only evaluate every 12 points (different points on the sphere)
        rot = self.base_rot[::12] # 48x3x3
        rot = rot.expand(B,*rot.shape) # Bx48x3x3
        x = self.model.cat_z(coords @ rot, z)
        y_hat = self.model.decoder.forward_symmetric(x,c) # Bx48xYX
        y_hat = y_hat.unsqueeze(2) # Bx48x1xYX

        images = images.unsqueeze(1) # Bx1xYxX
        images = torch.cat((images,rotated_images),1) # Bx12xYxX
        images = images.view(B,1,12,-1)[...,mask] # Bx1x12xYX

        err = torch.sum((images-y_hat).pow(2),-1).view(B,self.nbase)
        mini = torch.argmin(err,1)
        return mini.cpu().numpy()

    def opt_theta_rot(self, images, rotated_images, z, niter=5, L=None):
        B = images.size(0)
        assert not self.model.training
        with torch.no_grad():
            min_i = self.eval_base_grid(images, rotated_images, z, L)
            min_quat = self.base_quat[min_i]
            s2i, s1i = so3_grid.get_base_indr(min_i)
            for iter_ in range(1,niter+1):
                neighbors = [so3_grid.get_neighbor(min_quat[i], s2i[i], s1i[i], iter_) for i in range(B)]
                quat = np.array([x[0] for x in neighbors])
                ind = np.array([x[1] for x in neighbors])
                rot = lie_tools.quaternions_to_SO3(torch.tensor(quat))
                min_i = self.eval_grid(images, rot, z, 8, L)
                min_ind = ind[np.arange(B), min_i]
                s2i, s1i = min_ind.T
                min_quat = quat[np.arange(B),min_i]
        return lie_tools.quaternions_to_SO3(torch.tensor(min_quat))


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


