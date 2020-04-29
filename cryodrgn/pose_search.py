
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from cryodrgn import lie_tools
from cryodrgn import so3_grid
from cryodrgn import shift_grid
from cryodrgn import utils

log = utils.log
vlog = utils.vlog

#def unpack(self, packed, NQ, pad_val=float('inf')):
#    '''
#    Unpack a variable batch size data tensor into a padded B x max(NQ) tensor
#    
#    Inputs:
#        packed: (sum(NQ) x ...) data tensor
#        NQ: Bx1 # TODO : check Bx1 or B,
#
#    Returns: B x max(NQ)
#    '''
#    B = len(NQ)
#    new = torch.empty(B,max(NQ),*squashed.shape[1:],dtype=squashed.dtype)
#    prev = 0
#    for i in range(B):
#        new[i,0:NQ[i],...] = squashed[prev:(prev+NQ[i]),...]
#        new[i,NQ[i]:,...] = pad_val 
#        prev += NQ[i]
#    return new

class PoseSearch:
    '''Pose search'''
    def __init__(self, model, lattice, Lmin, Lmax, tilt=None, base_healpy=1, t_extent=5, t_ngrid=7, nkeptposes=24):
        if base_healpy != 1:
            raise NotImplementedError() # TODO

        self.model = model
        self.lattice = lattice
        self.base_quat = so3_grid.base_SO3_grid()
        self.base_rot = lie_tools.quaternions_to_SO3(torch.tensor(self.base_quat))
        self.nbase = len(self.base_quat)
        assert self.nbase == 576, "Base resolution changed?"
        self.base_shifts = torch.tensor(shift_grid.base_shift_grid(t_extent, t_ngrid)).float()
        self.t_extent = t_extent
        self.t_ngrid = t_ngrid
        self.Lmin = Lmin
        self.Lmax = Lmax
        self.tilt = tilt
        self.nkeptposes = nkeptposes

    def eval_grid(self, images, rot, z, NQ, L, images_tilt=None):
        '''
        images: B x T x Npix
        rot: (NxQ) x 3 x 3 rotation matrics (N=1 for base grid, N=B for incremental grid)
        NQ: number of slices evaluated for each image
        L: radius of fourier components to evaluate
        '''
        B = images.size(0)
        mask = self.lattice.get_circular_mask(L)
        coords = self.lattice.coords[mask]
        YX = coords.size(-2)
        def compute_err(images, rot):
            x = coords @ rot
            if z is not None:
                x = self.model.cat_z(x, z)
            y_hat = self.model(x)
            y_hat = y_hat.view(-1, 1, NQ, YX)  #1x1xNQxYX for base grid, Bx1x8xYX for incremental grid
            images = images.unsqueeze(2)  # BxTx1xYX
            err = torch.sum((images - y_hat).pow(2), -1)  # BxTxQ
            return err
        err = compute_err(images, rot)
        if images_tilt is not None:
            err_tilt = compute_err(images_tilt, self.tilt @ rot)
            err += err_tilt
        return err # BxTxQ

    def mask_images(self, images, L):
        '''
        images: B x NY x NX x 2
        Returns: B x Npix at resolution L
        '''
        B = images.size(0)
        mask = self.lattice.get_circular_mask(L)
        return images.view(B,-1)[:,mask]
    
    def translate_images(self, images, shifts, L):
        '''
        images: B x NY x NX
        shifts: B x T x 2 or B
        Returns: B x T x Npix at resolution L
        '''
        B = images.size(0)
        mask = self.lattice.get_circular_mask(L)
        return self.lattice.translate_ht(images.view(B,-1)[:,mask], shifts, mask)


    def subdivide(self, quat, q_ind, t_ind, cur_res):
        '''
        Subdivides poses for next resolution level

        Inputs:
            quat (N x 4 tensor): quaternions
            q_ind (N x 2 np.array): index of current S2xS1 grid
            t_ind (N x 2 np.array): index of current trans grid
            cur_res (int): Current resolution level
        
        Returns:
            quat  (N x 8 x 4) np.array
            q_ind (N x 8 x 2) np.array
            t_ind (N x 4 x 2) np.array
            rot   (N*8 x 3 x 3) tensor
            trans (N*4 x 2) tensor
        '''
        # get neighboring SO3 elements at next resolution level -- todo: make this an array operation
        neighbors = [so3_grid.get_neighbor(quat[i], q_ind[i][0], q_ind[i][1], cur_res) for i in range(len(quat))]
        quat = np.array([x[0] for x in neighbors]) # Bx8x4
        q_ind = np.array([x[1] for x in neighbors]) # Bx8x2
        rot = lie_tools.quaternions_to_SO3(torch.tensor(quat).view(-1,4))
        # get neighboring shifts at next resolution level -- todo: make this an array operation
        neighbors = [shift_grid.get_neighbor(xx,yy, cur_res-1, self.t_extent, self.t_ngrid) for xx,yy in t_ind]
        trans = torch.tensor(np.array([x[0] for x in neighbors]).reshape(-1,2))
        t_ind = np.array([x[1] for x in neighbors]) # Bx4x2
        #quat = np.tile(quat,(1,4,1)) # Bx8x4 -> Bx32x4
        #q_ind = np.tile(q_ind,(1,4,1)) # Bx8x2 -> Bx32x2
        #t_ind = np.repeat(t_ind,8,1) # Bx4x2 -> Bx32x2
        return quat, q_ind, t_ind, rot, trans
   
    # TODO: refactor -- no more L_star, just max_poses
    def keep_matrix(self, bound, max_poses):
        '''
        Returns:
            keep (bound.shape): bool tensor of poses to keep
        '''
        assert bound.ndimension() == 2
        B = bound.shape[0]
        keep_idx = bound.topk(max_poses, dim=-1)[1]
        keep = torch.zeros_like(bound).scatter_(dim=-1, index=keep_idx, value=1)
        return keep

    def opt_theta_trans(self, images, z=None, images_tilt=None, niter=5):
        B = images.size(0)
        assert not self.model.training

        # Expand the base grid B times if each image has a different z
        if z is not None:
            base_rot = self.base_rot.expand(B,*self.base_rot.shape) # B x 576 x 3 x 3
        else:
            base_rot = self.base_rot # 576 x 3 x 3

        # Compute the loss for all poses
        loss = self.eval_grid(self.translate_images(images, self.base_shifts, self.Lmin), 
                                        base_rot, z, self.nbase, self.Lmin, 
                                        images_tilt=self.translate_images(images_tilt, self.base_shifts, self.Lmin) if images_tilt is not None else None) # BxTxQ

        keep = self.keep_matrix(loss.view(B, -1), self.nkeptposes) # B x -1
        keep = keep.view(B, len(self.base_shifts), self.nbase) # BxTxQ
        Np = keep.sum((-1,-2)) # per image # this now constant per image
        w = keep.nonzero().cpu() # sum(Np) x 3; (B,T,Q)
        quat = self.base_quat[w[:,2]]
        s2i, s1i = so3_grid.get_base_ind(w[:,2])
        q_ind = np.stack((s2i, s1i), 1) # Np x 2
        trans = self.base_shifts[w[:,1]]
        xi, yi = shift_grid.get_base_ind(w[:,1], self.t_ngrid)
        t_ind = np.stack((xi, yi), 1) #Np x 2
        batch_ind = w[:,0]

        k = int((self.Lmax-self.Lmin)/(niter-1))
        for iter_ in range(1,niter+1):
            vlog(iter_); vlog(Np)
            L = min(self.Lmin +k*iter_, self.Lmax)
            quat, q_ind, t_ind, rot, trans = self.subdivide(quat, q_ind, t_ind, iter_)
            batch_ind4 = batch_ind.unsqueeze(1).repeat(1,4).view(-1) # repeat each element 4 times
            batch_ind8 = batch_ind.unsqueeze(1).repeat(1,8).view(-1) # repeat each element 8 times
            zb = z[batch_ind8] if z is not None else None
            loss = self.eval_grid(self.translate_images(images[batch_ind4], trans.unsqueeze(1), L).view(len(batch_ind),4,-1), rot, zb, 8, L,
                                        images_tilt=self.translate_images(images_tilt[batch_ind4],trans.unsqueeze(1), L).view(len(batch_ind),4,-1) if images_tilt is not None else None) # sum(NP),4x8
            keep = self.keep_matrix(loss.view(B, -1), self.nkeptposes)  # B x (self.Nkeptposes*32)
            keep = keep.view(B, self.nkeptposes, 4, 8)  # todo -- nkeptposes could be set to something else
            w = keep.nonzero() # NP x 4; (0-B , 0-24, 0-4, 0-8)
            assert len(w) == B*self.nkeptposes, f"{len(w)} != {B} x {self.nkeptposes} at iter {iter_}"
            ww = w[:,1] + w[:,0] * self.nkeptposes
            quat = quat[ww, w[:,3]]
            q_ind = q_ind[ww, w[:,3]]
            t_ind = t_ind[ww, w[:,2]]
        best = self.keep_matrix(loss.view(B,-1), 1)
        best_i = best.view(B, self.nkeptposes, 4, 8).nonzero()
        assert len(best_i) == B
        best_rot = rot[best_i[:,0]*self.nkeptposes + best_i[:,2]]
        best_trans = trans[best_i[:,0]*self.nkeptposes + best_i[:,1]]
        return best_rot, best_trans


