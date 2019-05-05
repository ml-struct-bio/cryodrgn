'''Branch and bound classes'''

# TODO: Update docstrings in this module!

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

import lie_tools
import so3_grid
import shift_grid
import utils

log = utils.log
vlog = utils.vlog

class BNBHomoRot:
    '''Branch and bound of rotations for homogeneous reconstruction'''
    def __init__(self, decoder, lattice, Lmin, Lmax, tilt=None):
        self.decoder = decoder
        self.lattice = lattice
        self.base_quat = so3_grid.base_SO3_grid()
        self.base_rot = lie_tools.quaternions_to_SO3(torch.tensor(self.base_quat))
        self.nbase = len(self.base_quat)
        assert self.nbase == 576, "Base resolution changed?"
        self.Lmin = Lmin
        self.Lmax = Lmax
        self.tilt = tilt

    def compute_B1(self, images, L, images_tilt=None):
        mask = self.lattice.get_circular_mask(self.Lmax) - self.lattice.get_circular_mask(L)
        power = images.view(images.size(0),-1)[:,mask].pow(2).sum(-1)
        if images_tilt:
            power += images_tilt.view(images_tilt.size(0),-1)[:,mask].pow(2).sum(-1)
        return power

    def compute_B2(self, max_slice, L):
        mask = self.lattice.get_circular_mask(self.Lmax) - self.lattice.get_circular_mask(L)
        power = max_slice[mask].pow(2).sum()
        B2 = -power - 4*(power)**.5
        return B2

    def estimate_max_slice(self):
        rot = self.base_rot[::24] # 48 slices
        y_hat = self.decoder(self.lattice.coords @ rot)
        return y_hat[y_hat.pow(2).sum(-1).argmax()] # YX

    def eval_grid(self, images, rot, NQ, L, images_tilt=None):
        '''
        images: B x NY x NX 
        rot: N x 3 x 3 rotation matrics (N=576 for base grid, N=B*8 for incremental grid)
        NQ: number of slices evaluated for each image, int or array of len B
        L: Max radius of fourier components to evaluate
        '''
        B = images.size(0)
        mask = self.lattice.get_circular_mask(L)
        coords = self.lattice.coords[mask]
        YX = coords.size(-2)
        def compute_err(images, rot):
            images = images.view(B,-1)[:,mask]
            y_hat = self.decoder(coords @ rot) # len(rot) x YX
            if type(NQ) == int: # constant NQ per image
                y_hat = y_hat.view(-1,NQ,YX) #1xQxYX for base grid, Bx8xYX for incremental grid
            else:
                # variable number of poses for each image - create a tensor of B x max(NQ) x YX
                # and tile y_hat into this tensor
                y_hat = self.tile(y_hat, NQ)
            images = images.unsqueeze(1) # Bx1xYX
            err = torch.sum((images-y_hat).pow(2),-1) # BxQ
            return err
        err = compute_err(images, rot)
        if images_tilt is not None:
            err_tilt = compute_err(images_tilt, self.tilt @ rot)
            err += err_tilt
        return err # B x Q

    def tile(self, squashed, NQ, nan_val=float('inf')):
        '''Tile a squashed, variable batch size data tensor into Bxmax(NQ) tensor'''
        B = len(NQ)
        new = torch.empty((B,max(NQ),squashed.shape[-1]),dtype=squashed.dtype)
        prev = 0
        for i in range(B):
            new[i,0:NQ[i],:] = squashed[prev:(prev+NQ[i]),:]
            new[i,NQ[i]:,:] = nan_val 
            prev += NQ[i]
        return new
    
    def tile_np(self, squashed, NQ, nan_val=float('inf')):
        '''Tile a squashed, variable batch size data tensor into Bxmax(NQ) tensor'''
        B = len(NQ)
        new = np.empty((B,max(NQ),squashed.shape[-1]),dtype=squashed.dtype)
        prev = 0
        for i in range(B):
            new[i,0:NQ[i],:] = squashed[prev:(prev+NQ[i]),:]
            new[i,NQ[i]:,:] = nan_val 
            prev += NQ[i]
        return new

    def compute_bound(self, images, rot, NQ, L, images_tilt=None, probabilistic=False):
        '''Compute the lower bound'''
        A = self.eval_grid(images, rot, NQ, L, images_tilt) # B x Q
        if L == self.Lmax: return A
        if probabilistic:
            B1 = self.compute_B1(images, L, images_tilt).unsqueeze(1) # Bx1
            B2 = self.compute_B2(self.max_slice, L)  # 1
        else:
            B1 = B2 = 0
        if images_tilt is not None: B2 *= 2
        return A+B1+B2

    # todo: clean up and rename this func
    def keep_matrix(self, bound, Lstar, max_poses, probabilistic=False):
        keep = bound <= Lstar # B x Q array of 0/1s
        NQ = keep.sum(1) # B, array of Nposes per batch element
        if (NQ > max_poses).any():
            # filter keep with percentile heuristic
            vlog('Warning: More than {} poses below upper bound')
            w = bound.argsort()[:,max_poses:]
            B,Q = bound.shape
            keep[np.arange(B).repeat(Q-max_poses), w.contiguous().view(-1)] *= 0
            NQ = keep.sum(1) # B, array of Nposes per batch element
        if (NQ == 0).any():
            assert probabilistic
            w = bound.argmin(-1)
            keep[np.arange(B), w] = 1     
            NQ = keep.sum(1) # B, array of Nposes per batch element
        return keep, NQ
       
    def bound_base(self, bound, Lstar, max_poses=24, probabilistic=False):
        '''Helper function to filter next poses to try'''
        keep, NQ = self.keep_matrix(bound, Lstar, max_poses, probabilistic)
        ### get squashed list of all poses to keep ###
        # keep.nonzero() returns Nx2 with N 2D indices of nonzero elements
        # we want the second index since that tells us the ind of the quat to keep
        # there may be duplicates in this list...
        min_i = keep.nonzero()[:,1].cpu().numpy() 
        quat = self.base_quat[min_i] # Qx4, where Q=sum(NQ)
        s2i, s1i = so3_grid.get_base_indr(min_i)
        grid_ind = np.stack((s2i, s1i),1) # Qx2, where Q=sum(NQ)
        assert len(quat) == NQ.sum()
        return NQ, quat, grid_ind

    def bound(self, bound, Lstar, quat_block, grid_ind_block, probabilistic=False):
        '''Helper function to filter next poses to try'''
        max_poses = int(bound.shape[1]*.125)
        keep, NQ = self.keep_matrix(bound, Lstar, max_poses, probabilistic)
        ### get squashed list of all poses to keep ###
        w = np.where(keep.cpu())
        quat = quat_block[w]
        grid_ind = grid_ind_block[w]
        assert len(quat) == NQ.sum()
        return NQ, quat, grid_ind

    def subdivide(self, NQ, quat, grid_ind, curr_res):
        '''Get quaternions and grid indices for next resolution level'''
        assert NQ.sum() == len(quat)
        assert len(quat) == len(grid_ind)
        NQ *= 8
        neighbors = [so3_grid.get_neighbor(quat[i], grid_ind[i][0], grid_ind[i][1], curr_res) for i in range(len(quat))]
        quat = np.array([x[0] for x in neighbors]).reshape(-1,4) # Qx8x4->8Qx4
        grid_ind = np.array([x[1] for x in neighbors]).reshape(-1,2) # Qx8x2->8Qx2
        assert NQ.sum() == len(quat)
        quat_block = self.tile_np(quat, NQ)
        grid_ind_block = self.tile_np(grid_ind, NQ)
        return NQ, quat, grid_ind, quat_block, grid_ind_block

    def opt_theta(self, images, L, images_tilt=None, niter=5, probabilistic=False):
        assert L is None # ew: to keep API consistent with BNNB
        B = images.size(0)
        assert not self.decoder.training
        with torch.no_grad():
            if probabilistic:
                self.max_slice = self.estimate_max_slice() # where to put this?
            bound = self.compute_bound(images, self.base_rot, self.nbase, self.Lmin, images_tilt=images_tilt, probabilistic=probabilistic)
            Lstar = self.eval_grid(images, self.base_rot[torch.argmin(bound,1)], 1, self.Lmax, images_tilt=images_tilt) # expensive objective function
            NQ, quat, grid_ind = self.bound_base(bound, Lstar, probabilistic=probabilistic)
            k = int((self.Lmax-self.Lmin)/(niter-1))
            for iter_ in range(1,niter+1): # resolution level
                L = min(self.Lmin + k*iter_, self.Lmax)
                NQ, quat, grid_ind, quat_block, grid_ind_block = self.subdivide(NQ, quat, grid_ind, iter_)
                rot = lie_tools.quaternions_to_SO3(torch.tensor(quat))
                bound = self.compute_bound(images, rot, NQ, L, images_tilt=images_tilt, probabilistic=probabilistic) # Bxmax(NQ)
                min_i = torch.argmin(bound, 1)
                min_quat = quat_block[np.arange(B),min_i.cpu()]
                Lstar = self.eval_grid(images, lie_tools.quaternions_to_SO3(torch.tensor(min_quat)), 1, self.Lmax, images_tilt=images_tilt)
                NQ, quat, grid_ind = self.bound(bound, Lstar, quat_block, grid_ind_block, probabilistic=probabilistic)
        return lie_tools.quaternions_to_SO3(torch.tensor(min_quat))

class BNBHomo:
    '''Branch and bound of rotation and translation for homogeneous reconstruction'''
    def __init__(self, decoder, lattice, Lmin, Lmax, tilt=None, t_extent=5):
        self.decoder = decoder
        self.lattice = lattice
        self.base_quat = so3_grid.base_SO3_grid()
        self.base_rot = lie_tools.quaternions_to_SO3(torch.tensor(self.base_quat))
        self.nbase = len(self.base_quat)
        assert self.nbase == 576, "Base resolution changed?"
        self.base_shifts = torch.tensor(shift_grid.base_shift_grid(t_extent)).float()
        self.t_extent = t_extent
        self.Lmin = Lmin
        self.Lmax = Lmax
        self.tilt = tilt

    def compute_B1(self, images, L, images_tilt=None):
        mask = self.lattice.get_circular_mask(self.Lmax) - self.lattice.get_circular_mask(L)
        power = images.view(images.size(0),-1)[:,mask].pow(2).sum(-1)
        if images_tilt:
            power += images_tilt.view(images_tilt.size(0),-1)[:,mask].pow(2).sum(-1)
        return power

    def compute_B2(self, max_slice, L):
        mask = self.lattice.get_circular_mask(self.Lmax) - self.lattice.get_circular_mask(L)
        power = max_slice[mask].pow(2).sum()
        B2 = -power - 4*(power)**.5
        return B2

    def estimate_max_slice(self):
        rot = self.base_rot[::24] # 48 slices
        y_hat = self.decoder(self.lattice.coords @ rot)
        return y_hat[y_hat.pow(2).sum(-1).argmax()] # YX

    def eval_grid(self, images, rot, NQ, L, images_tilt=None):
        '''
        images: B x T x Npix x 2
        rot: (NxQ) x 3 x 3 rotation matrics (N=1 for base grid, N=B for incremental grid)
        NQ: number of slices evaluated for each image
        L: radius of fourier components to evaluate
        '''
        B = images.size(0)
        mask = self.lattice.get_circular_mask(L)
        coords = self.lattice.coords[mask]
        YX = coords.size(-2)
        def compute_err(images, rot):
            y_hat = self.decoder(coords @ rot)
            y_hat = y_hat.view(-1,1,NQ,YX) #1x1xQxYX for base grid, Bx1x8xYXx2 for incremental grid
            images = images.unsqueeze(2) # BxTx1xYX
            err = torch.sum((images-y_hat).pow(2),-1) # BxTxQ
            return err
        err = compute_err(images, rot)
        if images_tilt is not None:
            err_tilt = compute_err(images_tilt, self.tilt @ rot)
            err += err_tilt
        return err # BxTxQ

    def shift_images(self, images, shifts, L):
        '''
        images: B x NY x NX x 2
        shifts: B x T x 2 or B x 2
        Returns: B x T x Npix x 2 at resolution L
        '''
        B = images.size(0)
        mask = self.lattice.get_circular_mask(L)
        coords = self.lattice.coords[mask,0:2]/2# 2D wavevector between -.5 and .5
        return self.decoder.translate_ht(coords, images.view(B,-1)[:,mask], shifts)

    def tile(self, squashed, NQ, nan_val=float('inf')):
        '''Tile a squashed, variable batch size data tensor into Bxmax(NQ) tensor'''
        B = len(NQ)
        new = torch.empty(B,max(NQ),*squashed.shape[1:],dtype=squashed.dtype)
        prev = 0
        for i in range(B):
            new[i,0:NQ[i],:] = squashed[prev:(prev+NQ[i]),:]
            new[i,NQ[i]:,:] = nan_val 
            prev += NQ[i]
        return new
 
    def subdivide(self, quat, q_ind, t_ind, curr_res):
        # get neighboring SO3 elements at next resolution level -- todo: make this an array operation
        neighbors = [so3_grid.get_neighbor(quat[i], q_ind[i][0], q_ind[i][1], curr_res) for i in range(len(quat))]
        quat = np.array([x[0] for x in neighbors]) # Bx8x4
        q_ind = np.array([x[1] for x in neighbors]) # Bx8x2
        rot = lie_tools.quaternions_to_SO3(torch.tensor(quat).view(-1,4))
        # get neighboring shifts at next resolution level -- todo: make this an array operation
        neighbors = [shift_grid.get_neighbor(xx,yy, curr_res-1, self.t_extent) for xx,yy in t_ind]
        trans = torch.tensor(np.array([x[0] for x in neighbors]).reshape(-1,2))
        t_ind = np.array([x[1] for x in neighbors]) # Bx4x2

        quat = np.tile(quat,(1,4,1)) # Bx8x4 -> Bx32x4
        q_ind = np.tile(q_ind,(1,4,1)) # Bx8x2 -> Bx32x2
        t_ind = np.repeat(t_ind,8,1) # Bx4x2 -> Bx32x2
        return quat, q_ind, rot, trans, t_ind

    def keep_matrix(self, bound, Lstar, max_poses, probabilistic):
        if probabilistic: raise NotImplementedError
        B = bound.shape[0]
        bound = bound.view(B,-1)
        keep = bound <= Lstar
        Np = keep.sum(1)
        if (Np > max_poses).any():
            w = bound.argsort(1)[:,max_poses:].cpu().contiguous().view(-1)
            maxQ = keep.shape[1]
            keep[np.arange(B).repeat(maxQ-max_poses), w] *= 0
        if (Np == 0).any():
            assert probabilistic
            w = bound.argmin(1)
            keep[np.arange(B),w] = 1
        return keep
            
    def opt_theta_trans(self, images, L, images_tilt=None, niter=5, probabilistic=False):
        B = images.size(0)
        assert not self.decoder.training
        with torch.no_grad():
            if probabilistic:
                self.max_slice = self.estimate_max_slice()
            # todo: check shift images
            bound = self.eval_grid(self.shift_images(images, self.base_shifts, self.Lmin), 
                                            self.base_rot, self.nbase, self.Lmin, 
                                            images_tilt=self.shift_images(images_tilt, self.base_shifts, self.Lmin) if images_tilt is not None else None) # BxTxQ
            if probabilistic: # bound += B1 + B2
                raise NotImplementedError
            mini = torch.argmin(bound.view(B,-1),1)
            qi = mini % self.nbase
            ti = mini // self.nbase
            min_trans = self.base_shifts[ti] # unsqueeze by 1 to get right dim for translate_ht func
            min_rot = self.base_rot[qi]
            Lstar = self.eval_grid(self.shift_images(images, min_trans.unsqueeze(1), self.Lmax),
                                        min_rot, 1, self.Lmax,
                                        images_tilt = self.shift_images(images_tilt, min_trans.unsqueeze(1), self.Lmax) if images_tilt is not None else None)
            keep = self.keep_matrix(bound, Lstar.view(B,1), 24, probabilistic) # B
            keep = keep.view(B,len(self.base_shifts), self.nbase) # BxTxQ
            Np = keep.sum((-1,-2)) # per image # todo: filter keep by Np to max_poses
            w = keep.nonzero().cpu() # sum(Np) x 3
            # todo: compare len(w) vs. len(set(w[:,2])) ... if much smaller then we should think about potential speedups
            quat = self.base_quat[w[:,2]]
            s2i, s1i = so3_grid.get_base_ind(w[:,2])
            q_ind = np.stack((s2i,s1i),1) # Np x 2
            trans = self.base_shifts[w[:,1]]
            xi, yi = shift_grid.get_base_ind(w[:,1])
            t_ind = np.stack((xi,yi),1) #Np x 2
            batch_ind = w[:,0]

            k = int((self.Lmax-self.Lmin)/(niter-1))
            for iter_ in range(1,niter+1):
                vlog(iter_); vlog(Np)
                L = min(self.Lmin +k*iter_, self.Lmax)
                quat, q_ind, rot, trans, t_ind = self.subdivide(quat, q_ind, t_ind, iter_)
                batch_ind4 = batch_ind.unsqueeze(1).repeat(1,4).view(-1) # repeat each element 4 times
                bound = self.eval_grid(self.shift_images(images[batch_ind4], trans.unsqueeze(1), L).view(len(batch_ind),4,-1), rot, 8, L,
                                            images_tilt=self.shift_images(images_tilt[batch_ind4],trans.unsqueeze(1),L).view(len(batch_ind),4,-1) if images_tilt is not None else None) # sum(NP),4x8
                bound2 = self.tile(bound, Np) # Bxmax(Np)x4x8
                min_i = bound2.view(B,-1).argmin(1) 
                min_i[1:] += 32*Np.cumsum(0)[0:-1]
                min_rot = rot[min_i//32*8+min_i%8]
                min_trans = trans[min_i/8]
                Lstar = self.eval_grid(self.shift_images(images, min_trans.unsqueeze(1), self.Lmax), min_rot, 1, self.Lmax,
                                        images_tilt=self.shift_images(images_tilt, min_trans.unsqueeze(1), self.Lmax) if images_tilt is not None else None) # Bx1x1
                keep = self.keep_matrix(bound2, Lstar.view(B,1), bound2.shape[1], probabilistic) # Bx(max(Np)*32)
                w = keep.nonzero() # sum(Np) x 2
                batch_ind = w[:,0]
                tmp = 32*Np.cumsum(0)
                tmp[-1] = 0
                w[:,1] += tmp[(w[:,0]-1)]
                w=w.cpu()
                quat = quat.reshape(-1,4)[w[:,1]].reshape(-1,4) # final reshape needed to keep 2D array if len(w) == 1
                q_ind = q_ind.reshape(-1,2)[w[:,1]].reshape(-1,2)
                t_ind = t_ind.reshape(-1,2)[w[:,1]].reshape(-1,2)
                Np = keep.sum(1)
                assert Np.sum() == len(t_ind)
        return min_rot, min_trans

class BNNBHomo:
    '''Branch and no bound for homogeneous reconstruction'''
    def __init__(self, decoder, lattice, tilt=None, t_extent=5):
        self.decoder = decoder
        self.lattice = lattice
        self.base_quat = so3_grid.base_SO3_grid()
        self.base_rot = lie_tools.quaternions_to_SO3(torch.tensor(self.base_quat))
        self.nbase = len(self.base_quat)
        self.base_shifts = torch.tensor(shift_grid.base_shift_grid(t_extent)).float()
        self.t_extent = t_extent
        assert self.nbase == 576, "Base resolution changed?"
        self.tilt = tilt

    def eval_grid(self, images, rot, NQ, L, images_tilt=None):
        '''
        images: B x T x Npix x 2
        rot: (NxQ) x 3 x 3 rotation matrics (N=1 for base grid, N=B for incremental grid)
        NQ: number of slices evaluated for each image
        L: radius of fourier components to evaluate
        '''
        B = images.size(0)
        mask = self.lattice.get_circular_mask(L)
        coords = self.lattice.coords[mask]
        YX = coords.size(-2)
        def compute_err(images, rot):
            y_hat = self.decoder(coords @ rot)
            y_hat = y_hat.view(-1,1,NQ,YX) #1x1xQxYXx2 for base grid, Bx1x8xYXx2 for incremental grid
            images = images.unsqueeze(2) # BxTx1xYXx2
            err = torch.sum((images-y_hat).pow(2),-1) # BxTxQ
            return err
        err = compute_err(images, rot)
        if images_tilt is not None:
            err_tilt = compute_err(images_tilt, self.tilt @ rot)
            err += err_tilt
        mini = torch.argmin(err.view(B,-1),1) # B
        qi = mini % NQ
        ti = mini / NQ # integer division
        return qi.cpu().numpy(), ti.cpu().numpy()

    def opt_theta(self, images, L, images_tilt=None, niter=5):
        B = images.size(0)
        # mask images
        images = self.mask_images(images, L)
        images = images.unsqueeze(1) # not testing any translations
        if images_tilt is not None:
            images_tilt = self.mask_images(images_tilt,L)
            images_tilt = images_tilt.unsqueeze(1) # not testing any translations
        assert not self.decoder.training
        with torch.no_grad():
            min_i, _ = self.eval_grid(images, self.base_rot, self.nbase, L, images_tilt=images_tilt)
            min_quat = self.base_quat[min_i]
            s2i, s1i = so3_grid.get_base_indr(min_i)
            for iter_ in range(1,niter+1):
                neighbors = [so3_grid.get_neighbor(min_quat[i], s2i[i], s1i[i], iter_) for i in range(B)]
                quat = np.array([x[0] for x in neighbors])
                ind = np.array([x[1] for x in neighbors])
                rot = lie_tools.quaternions_to_SO3(torch.tensor(quat))
                min_i, _ = self.eval_grid(images, rot, 8, L, images_tilt=images_tilt)
                min_ind = ind[np.arange(B), min_i]
                s2i, s1i = min_ind.T
                min_quat = quat[np.arange(B),min_i]
        return lie_tools.quaternions_to_SO3(torch.tensor(min_quat))

    def mask_images(self, images, L):
        mask = self.lattice.get_circular_mask(L)
        images = images.view(images.size(0),-1)[:,mask]
        return images

    def shift_images(self, images, shifts, L):
        '''
        images: B x N x 2
        shifts: B x T x 2 or Tx2
        Returns: B x T x N x 2
        '''
        mask = self.lattice.get_circular_mask(L)
        coords = self.lattice.coords[mask,0:2]/2# 2D wavevector between -.5 and .5
        return self.decoder.translate_ht(coords, images, shifts)

    def opt_theta_trans(self, images, L, images_tilt=None, niter=5):
        B = images.size(0)
        # mask images
        images = self.mask_images(images, L)
        if images_tilt is not None:
            images_tilt = self.mask_images(images_tilt,L)
        assert not self.decoder.training
        with torch.no_grad():
            min_qi, min_ti = self.eval_grid(self.shift_images(images, self.base_shifts, L), 
                                            self.base_rot, self.nbase, L, 
                                            images_tilt=self.shift_images(images_tilt, self.base_shifts, L) if images_tilt is not None else None)
            min_quat = self.base_quat[min_qi]
            min_t = self.base_shifts[min_ti]
            s2i, s1i = so3_grid.get_base_ind(min_qi)
            xi, yi = shift_grid.get_base_ind(min_ti) # check if we need to cast as int
            for iter_ in range(1,niter+1):
                # get neighboring SO3 elements at next resolution level -- todo: make this an array operation
                neighbors = [so3_grid.get_neighbor(min_quat[i], s2i[i], s1i[i], iter_) for i in range(B)]
                quat = np.array([x[0] for x in neighbors])
                ind = np.array([x[1] for x in neighbors])
                rot = lie_tools.quaternions_to_SO3(torch.tensor(quat))
                # get neighboring shifts at next resolution level -- todo: make this an array operation
                neighbors = [shift_grid.get_neighbor(xx,yy, iter_, self.t_extent) for xx,yy in zip(xi,yi)]
                trans = torch.tensor(np.array([x[0] for x in neighbors]))
                t_ind = np.array([x[1] for x in neighbors])
                
                min_qi, min_ti = self.eval_grid(self.shift_images(images, trans, L), rot, 8, L, 
                                                images_tilt=self.shift_images(images_tilt, trans, L) if images_tilt is not None else None)

                # set new s2i, s1i, min_quat, xi, yi, min_trans
                min_ind = ind[np.arange(B), min_qi]
                s2i, s1i = min_ind.T
                min_quat = quat[np.arange(B),min_qi]
                xi, yi = t_ind[np.arange(B),min_ti].T # check this
                min_trans = trans[np.arange(B), min_ti] # check this, Bx2?
        return lie_tools.quaternions_to_SO3(torch.tensor(min_quat)), min_trans

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
        YX = coords.size(-2)
        def compute_err(images, rot):
            images = images.view(B,-1)[:,mask]
            x = self.model.cat_z(coords @ rot, z)
            y_hat = self.model.decoder(x)
            y_hat = y_hat.view(B, NQ, YX) # BxQxYXx2
            images = images.unsqueeze(1) # Bx1xYXx2
            err = torch.sum((images-y_hat).pow(2),-1) # BxQ
            return err
        err = compute_err(images,rot)
        if images_tilt is not None:
            err_tilt = compute_err(images_tilt, self.tilt @ rot)
            err += err_tilt
        mini = torch.argmin(err,1)
        return mini.cpu().numpy()
        
    # for rotation speed-up trick
    def eval_base_grid(self, images, rotated_images, z, L):
        '''
        images: BxYxX
        rotated_images:Bx11xYxX
        '''
        B = z.size(0)
        mask = self.lattice.get_circular_mask(L)
        coords = self.lattice.coords[mask]
        YX = coords.size(-2)

        # only evaluate every 12 points (different points on the sphere)
        rot = self.base_rot[::12] # 48x3x3
        rot = rot.expand(B,*rot.shape) # Bx48x3x3
        x = self.model.cat_z(coords @ rot, z)
        y_hat = self.model.decoder(x) # Bx48xYX
        y_hat = y_hat.unsqueeze(2) # Bx48x1xYX

        images = images.unsqueeze(1) # Bx1xYxX
        images = torch.cat((images,rotated_images),1) # Bx12xYxX
        images = images.view(B,1,12,-1)[...,mask] # Bx1x12xYX

        err = torch.sum((images-y_hat).pow(2),-1).view(B,self.nbase)
        mini = torch.argmin(err,1)
        return mini.cpu().numpy()

    # for rotation speed-up trick
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

    def opt_theta_trans(self, *args, **kwargs):
        raise NotImplementedError

class BNBHet:
    '''Branch and bound over rotation and translation for heterogeneous reconstruction'''
    def __init__(self, model, lattice, Lmin, Lmax, tilt=None, t_extent=5):
        self.model = model
        self.lattice = lattice
        self.base_quat = so3_grid.base_SO3_grid()
        self.base_rot = lie_tools.quaternions_to_SO3(torch.tensor(self.base_quat))
        self.nbase = len(self.base_quat)
        assert self.nbase == 576, "Base resolution changed?"
        self.base_shifts = torch.tensor(shift_grid.base_shift_grid(t_extent)).float()
        self.t_extent = t_extent
        self.Lmin = Lmin
        self.Lmax = Lmax
        self.tilt = tilt

    def compute_B1(self, images, L, images_tilt=None):
        mask = self.lattice.get_circular_mask(self.Lmax) - self.lattice.get_circular_mask(L)
        power = images.view(images.size(0),-1)[:,mask].pow(2).sum(-1)
        if images_tilt:
            power += images_tilt.view(images_tilt.size(0),-1)[:,mask].pow(2).sum(-1)
        return power

    def compute_B2(self, max_slice, L):
        mask = self.lattice.get_circular_mask(self.Lmax) - self.lattice.get_circular_mask(L)
        power = max_slice[mask].pow(2).sum()
        B2 = -power - 4*(power)**.5
        return B2

    def estimate_max_slice(self):
        rot = self.base_rot[::24] # 48 slices
        y_hat = self.model.decoder(self.lattice.coords @ rot)
        return y_hat[y_hat.pow(2).sum(-1).argmax()] # YX

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
            x = self.model.cat_z(coords @ rot, z)
            y_hat = self.model.decoder(x)
            y_hat = y_hat.view(-1,1,NQ,YX) #1x1xQxYX for base grid, Bx1x8xYXx2 for incremental grid
            images = images.unsqueeze(2) # BxTx1xYX
            err = torch.sum((images-y_hat).pow(2),-1) # BxTxQ
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
    
    def shift_images(self, images, shifts, L):
        '''
        images: B x NY x NX
        shifts: B x T x 2 or B
        Returns: B x T x Npix at resolution L
        '''
        B = images.size(0)
        mask = self.lattice.get_circular_mask(L)
        coords = self.lattice.coords[mask,0:2]/2# 2D wavevector between -.5 and .5
        return self.model.decoder.translate_ht(coords, images.view(B,-1)[:,mask], shifts)

    def tile(self, squashed, NQ, nan_val=float('inf')):
        '''Tile a squashed, variable batch size data tensor into Bxmax(NQ) tensor'''
        B = len(NQ)
        new = torch.empty(B,max(NQ),*squashed.shape[1:],dtype=squashed.dtype)
        prev = 0
        for i in range(B):
            new[i,0:NQ[i],:] = squashed[prev:(prev+NQ[i]),:]
            new[i,NQ[i]:,:] = nan_val 
            prev += NQ[i]
        return new
 
    def subdivide(self, quat, q_ind, t_ind, curr_res):
        # get neighboring SO3 elements at next resolution level -- todo: make this an array operation
        neighbors = [so3_grid.get_neighbor(quat[i], q_ind[i][0], q_ind[i][1], curr_res) for i in range(len(quat))]
        quat = np.array([x[0] for x in neighbors]) # Bx8x4
        q_ind = np.array([x[1] for x in neighbors]) # Bx8x2
        rot = lie_tools.quaternions_to_SO3(torch.tensor(quat).view(-1,4))
        # get neighboring shifts at next resolution level -- todo: make this an array operation
        neighbors = [shift_grid.get_neighbor(xx,yy, curr_res-1, self.t_extent) for xx,yy in t_ind]
        trans = torch.tensor(np.array([x[0] for x in neighbors]).reshape(-1,2))
        t_ind = np.array([x[1] for x in neighbors]) # Bx4x2
        quat = np.tile(quat,(1,4,1)) # Bx8x4 -> Bx32x4
        q_ind = np.tile(q_ind,(1,4,1)) # Bx8x2 -> Bx32x2
        t_ind = np.repeat(t_ind,8,1) # Bx4x2 -> Bx32x2
        return quat, q_ind, rot, trans, t_ind
    
    def subdivide_rotonly(self, quat, q_ind, curr_res):
        # get neighboring SO3 elements at next resolution level -- todo: make this an array operation
        neighbors = [so3_grid.get_neighbor(quat[i], q_ind[i][0], q_ind[i][1], curr_res) for i in range(len(quat))]
        quat = np.array([x[0] for x in neighbors]) # Bx8x4
        q_ind = np.array([x[1] for x in neighbors]) # Bx8x2
        rot = lie_tools.quaternions_to_SO3(torch.tensor(quat).view(-1,4))
        return quat, q_ind, rot
    
    # TODO: docstring this
    def keep_matrix(self, bound, Lstar, max_poses, probabilistic):
        if probabilistic: raise NotImplementedError
        B = bound.shape[0]
        bound = bound.view(B,-1)
        keep = bound <= Lstar
        Np = keep.sum(1)
        if (Np > max_poses).any():
            w = bound.argsort(1)[:,max_poses:].cpu().contiguous().view(-1)
            maxQ = keep.shape[1]
            keep[np.arange(B).repeat(maxQ-max_poses), w] *= 0
        if (Np == 0).any():
            assert probabilistic
            w = bound.argmin(1)
            keep[np.arange(B),w] = 1
        return keep
            
    def opt_theta(self, images, z, images_tilt=None, niter=5, L=None, probabilistic=False):
        assert L is None # ew: hack for bnnb API consistency
        if probabilistic: # bound += B1 + B2
            raise NotImplementedError
        B = images.size(0)
        assert not self.model.training
        with torch.no_grad():
            if probabilistic:
                self.max_slice = self.estimate_max_slice()
            # expand the base grid B times since each image has a different z
            base_rot = self.base_rot.expand(B,*self.base_rot.shape) # B x 576 x 3 x 3
            bound = self.eval_grid(self.mask_images(images, self.Lmin).unsqueeze(1), 
                                            base_rot, z, self.nbase, self.Lmin, 
                                            images_tilt=self.mask_images(images_tilt, self.Lmin).unsqueeze(1) if images_tilt is not None else None) # Bx1xQ
            # TODO: if probabilistic: bound += B1 + B2
            qi = torch.argmin(bound.view(B,self.nbase),1)
            min_rot = self.base_rot[qi]
            Lstar = self.eval_grid(self.mask_images(images, self.Lmax).unsqueeze(1),
                                        min_rot, z, 1, self.Lmax,
                                        images_tilt = self.mask_images(images_tilt, self.Lmax).unsqueeze(1) if images_tilt is not None else None)
            keep = self.keep_matrix(bound, Lstar.view(B,1), 24, probabilistic) # BxQ
            Np = keep.sum(-1) # per image # todo: filter keep by Np to max_poses
            w = keep.nonzero().cpu() # sum(Np) x 2
            quat = self.base_quat[w[:,1]]
            s2i, s1i = so3_grid.get_base_ind(w[:,1])
            q_ind = np.stack((s2i,s1i),1) # Np x 2
            batch_ind = w[:,0]

            k = int((self.Lmax-self.Lmin)/(niter-1))
            for iter_ in range(1,niter+1):
                vlog(iter_); vlog(Np)
                L = min(self.Lmin +k*iter_, self.Lmax)
                quat, q_ind, rot, = self.subdivide_rotonly(quat, q_ind, iter_)
                batch_ind8 = batch_ind.unsqueeze(1).repeat(1,8).view(-1) # repeat each element 8 times
                bound = self.eval_grid(self.mask_images(images[batch_ind], L).unsqueeze(1), rot, z[batch_ind8], 8, L,
                                            images_tilt=self.mask_images(images_tilt[batch_ind], L).unsqueeze(1) if images_tilt is not None else None) # sum(NP)x1x8
                # TODO: if probabilistic: bound += B1 + B2
                bound2 = self.tile(bound.view(-1,8), Np) # Bxmax(Np)x8?
                min_i = bound2.view(B,-1).argmin(1) 
                min_i[1:] += 8*Np.cumsum(0)[0:-1]
                min_rot = rot[min_i] # CHECK THIS IDK
                Lstar = self.eval_grid(self.mask_images(images, self.Lmax).unsqueeze(1), min_rot, z, 1, self.Lmax,
                                        images_tilt=self.mask_images(images_tilt, self.Lmax).unsqueeze(1) if images_tilt is not None else None) # Bx1x1
                keep = self.keep_matrix(bound2, Lstar.view(B,1), bound2.shape[1], probabilistic) # Bx(max(Np)*8)
                w = keep.nonzero() # sum(Np) x 2
                batch_ind = w[:,0]
                # CHECK THIS IDK
                tmp = 8*Np.cumsum(0)
                tmp[-1] = 0
                w[:,1] += tmp[(w[:,0]-1)]
                w=w.cpu()
                quat = quat.reshape(-1,4)[w[:,1]]
                q_ind = q_ind.reshape(-1,2)[w[:,1]]
                Np = keep.sum(1)
                assert Np.sum() == len(q_ind)
        return min_rot

    def opt_theta_trans(self, images, z, images_tilt=None, niter=5, L=None, probabilistic=False):
        assert L is None # ew: hack for bnnb API consistency
        if probabilistic: # bound += B1 + B2
            raise NotImplementedError
        B = images.size(0)
        assert not self.model.training
        with torch.no_grad():
            if probabilistic:
                self.max_slice = self.estimate_max_slice()
            # expand the base grid B times since each image has a different z
            base_rot = self.base_rot.expand(B,*self.base_rot.shape) # B x 576 x 3 x 3
            bound = self.eval_grid(self.shift_images(images, self.base_shifts, self.Lmin), 
                                            base_rot, z, self.nbase, self.Lmin, 
                                            images_tilt=self.shift_images(images_tilt, self.base_shifts, self.Lmin) if images_tilt is not None else None) # BxTxQ
            # TODO: if probabilistic: bound += B1 + B2
            mini = torch.argmin(bound.view(B,-1),1)
            qi = mini % self.nbase
            ti = mini // self.nbase
            min_trans = self.base_shifts[ti] # unsqueeze by 1 to get right dim for translate_ht func
            min_rot = self.base_rot[qi]
            Lstar = self.eval_grid(self.shift_images(images, min_trans.unsqueeze(1), self.Lmax),
                                        min_rot, z, 1, self.Lmax,
                                        images_tilt = self.shift_images(images_tilt, min_trans.unsqueeze(1), self.Lmax) if images_tilt is not None else None)
            keep = self.keep_matrix(bound, Lstar.view(B,1), 24, probabilistic) # Bx-1
            keep = keep.view(B,len(self.base_shifts), self.nbase) # BxTxQ
            Np = keep.sum((-1,-2)) # per image # todo: filter keep by Np to max_poses
            w = keep.nonzero().cpu() # sum(Np) x 3
            quat = self.base_quat[w[:,2]]
            s2i, s1i = so3_grid.get_base_ind(w[:,2])
            q_ind = np.stack((s2i,s1i),1) # Np x 2
            trans = self.base_shifts[w[:,1]]
            xi, yi = shift_grid.get_base_ind(w[:,1])
            t_ind = np.stack((xi,yi),1) #Np x 2
            batch_ind = w[:,0]

            k = int((self.Lmax-self.Lmin)/(niter-1))
            for iter_ in range(1,niter+1):
                vlog(iter_); vlog(Np)
                L = min(self.Lmin +k*iter_, self.Lmax)
                quat, q_ind, rot, trans, t_ind = self.subdivide(quat, q_ind, t_ind, iter_)
                batch_ind4 = batch_ind.unsqueeze(1).repeat(1,4).view(-1) # repeat each element 4 times
                batch_ind8 = batch_ind.unsqueeze(1).repeat(1,8).view(-1) # repeat each element 8 times
                bound = self.eval_grid(self.shift_images(images[batch_ind4], trans.unsqueeze(1), L).view(len(batch_ind),4,-1), rot, z[batch_ind8], 8, L,
                                            images_tilt=self.shift_images(images_tilt[batch_ind4],trans.unsqueeze(1), L).view(len(batch_ind),4,-1) if images_tilt is not None else None) # sum(NP),4x8
                # TODO: if probabilistic: bound += B1 + B2
                bound2 = self.tile(bound, Np) # Bxmax(Np)x4x8
                min_i = bound2.view(B,-1).argmin(1) 
                min_i[1:] += 32*Np.cumsum(0)[0:-1]
                min_rot = rot[min_i//32*8+min_i%8]
                min_trans = trans[min_i/8]
                Lstar = self.eval_grid(self.shift_images(images, min_trans.unsqueeze(1), self.Lmax), min_rot, z, 1, self.Lmax,
                                        images_tilt=self.shift_images(images_tilt, min_trans.unsqueeze(1), self.Lmax) if images_tilt is not None else None) # Bx1x1
                keep = self.keep_matrix(bound2, Lstar.view(B,1), bound2.shape[1], probabilistic) # Bx(max(Np)*32)
                w = keep.nonzero() # sum(Np) x 2
                batch_ind = w[:,0]
                tmp = 32*Np.cumsum(0)
                tmp[-1] = 0
                w[:,1] += tmp[(w[:,0]-1)]
                w=w.cpu()
                quat = quat.reshape(-1,4)[w[:,1]].reshape(-1,4) # final reshape needed to keep 2D array if len(w) == 1
                q_ind = q_ind.reshape(-1,2)[w[:,1]].reshape(-1,2)
                t_ind = t_ind.reshape(-1,2)[w[:,1]].reshape(-1,2)
                Np = keep.sum(1)
                assert Np.sum() == len(t_ind)
        return min_rot, min_trans


