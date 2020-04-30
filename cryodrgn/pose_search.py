
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

    def eval_grid(self, *, images, rot, z, NQ, L, images_tilt=None):
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
            err = (images - y_hat).pow(2).sum(-1)  # BxTxQ  # FIXME switch to correlation
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
        assert len(quat.shape) == 2 and quat.shape[1] == 4
        assert len(q_ind.shape) == 2 and q_ind.shape[1] == 2
        assert len(t_ind.shape) == 2 and t_ind.shape[1] == 2

        # get neighboring SO3 elements at next resolution level -- todo: make this an array operation
        neighbors = [so3_grid.get_neighbor(quat[i], q_ind[i][0], q_ind[i][1], cur_res) for i in range(len(quat))]
        quat = np.array([x[0] for x in neighbors]) # Bx8x4
        q_ind = np.array([x[1] for x in neighbors]) # Bx8x2
        rot = lie_tools.quaternions_to_SO3(torch.from_numpy(quat).view(-1,4))
       
        # get neighboring translations at next resolution level -- todo: make this an array operation
        neighbors = [shift_grid.get_neighbor(xx, yy, cur_res-1, self.t_extent, self.t_ngrid) for xx, yy in t_ind]
        trans = torch.tensor(np.array([x[0] for x in neighbors]).reshape(-1,2))
        t_ind = np.array([x[1] for x in neighbors]) # Bx4x2

        assert len(quat.shape) == 3 and quat.shape[1:] == (8, 4), quat.shape
        assert len(q_ind.shape) == 3 and q_ind.shape[1:] == (8, 2), q_ind.shape
        assert len(t_ind.shape) == 3 and t_ind.shape[1:] == (4, 2), t_ind.shape

        return quat, q_ind, t_ind, rot, trans
   
    def keep_matrix(self, loss, B, max_poses):
        '''
        Returns:
            keep (loss.shape): bool tensor of poses to keep
        '''
        flat_loss = loss.view(B, -1)
        keep_idx = flat_loss.topk(max_poses, dim=-1)[1]
        keep = torch.zeros_like(flat_loss).scatter_(dim=-1, index=keep_idx, value=1)
        return keep.view(loss.shape)

    def opt_theta_trans(self, images, z=None, images_tilt=None, niter=5):
        do_tilt = images_tilt is not None

        B = images.size(0)
        assert not self.model.training

        # Expand the base grid B times if each image has a different z
        if z is not None:
            base_rot = self.base_rot.expand(B,*self.base_rot.shape) # B x 576 x 3 x 3
        else:
            base_rot = self.base_rot # 576 x 3 x 3

        # Compute the loss for all poses 
        loss = self.eval_grid(
            images=self.translate_images(images, self.base_shifts, self.Lmin),
            rot=base_rot,
            z=z,
            NQ=self.nbase,
            L=self.Lmin, 
            images_tilt=self.translate_images(images_tilt, self.base_shifts, self.Lmin) if do_tilt else None
        )
        keep = self.keep_matrix(loss, B, self.nkeptposes) # B x -1
        keepB, keepT, keepQ = keep.nonzero().cpu().t()
        quat = self.base_quat[keepQ]
        q_ind = so3_grid.get_base_ind(keepQ)  # Np x 2
        trans = self.base_shifts[keepT]
        t_ind = shift_grid.get_base_ind(keepT, self.t_ngrid) #  Np x 2
        
        keepB4 = torch.arange(B).repeat_interleave(4 * self.nkeptposes)
        keepB8 = torch.arange(B).repeat_interleave(8 * self.nkeptposes)
        zb = z[keepB8] if z is not None else None
        k = int((self.Lmax - self.Lmin) / (niter - 1))
        for iter_ in range(1,niter+1):
            vlog(iter_); # vlog(Np)
            L = min(self.Lmin +k*iter_, self.Lmax)
            quat, q_ind, t_ind, rot, trans = self.subdivide(quat, q_ind, t_ind, iter_)

            loss = self.eval_grid(
                images=self.translate_images(images[keepB4], trans.unsqueeze(1), L).view(len(keepB),4,-1),  # (B*24, 4, Npoints)
                rot=rot,
                z=zb,
                NQ=8,
                L=L,
                images_tilt=self.translate_images(images_tilt[keepB4],trans.unsqueeze(1), L).view(len(keepB),4,-1) if do_tilt else None # (B*24, 4, Npoints)
            ) # sum(NP),4x8
            keep = self.keep_matrix(loss, B, self.nkeptposes)  # B x (self.Nkeptposes*32)
            keepB, keepT, keepQ = keep.nonzero().t() # NP x 4; (0-B * 24, 0-4, 0-8)
            assert len(keepB) == B * self.nkeptposes, f"{len(w)} != {B} x {self.nkeptposes} at iter {iter_}"
            quat = quat[keepB, keepQ]
            q_ind = q_ind[keepB, keepQ]
            t_ind = t_ind[keepB, keepT]
        
        best = self.keep_matrix(loss, B, 1)
        bestBN, bestT, bestQ = best.nonzero().t()
        assert len(bestBN) == B
        best_rot = rot.view(B * self.nkeptposes, 8, 3, 3)[bestBN, bestQ]
        best_trans = trans.view(B * self.nkeptposes, 4, 2)[bestBN, bestT]
        return best_rot, best_trans


