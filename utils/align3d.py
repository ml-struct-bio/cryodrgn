'''
Hacky script to perform global 3D alignment of two volumes
'''

import argparse
import numpy as np
import sys, os
import time
import pickle

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.utils.data as data

from cryodrgn import utils
from cryodrgn import mrc
from cryodrgn import fft
from cryodrgn import lie_tools
from cryodrgn import so3_grid
from cryodrgn import shift_grid3

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('ref', help='Input volume to align on')
    parser.add_argument('vol', help='Input volume to align')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Aligned mrc')
    parser.add_argument('--niter', type=int, default=20, help='Number of iterations of grid search (default: %(default)s)')
    parser.add_argument('--max-D', type=int, default=40, help='Box size to align')
    parser.add_argument('-rb', type=int, default=300, help='Rotation iterator batch size (default: %(default)s)')
    parser.add_argument('-tb', type=int, default=5, help='Translation iterator batch size (default: %(default)s)')
    parser.add_argument('--t-extent', type=float, help='Extent of volume translation (default: +/-D/16 pixels)')
    parser.add_argument('--t-grid', type=int, default=8, help='Gridsize per dimension')
    parser.add_argument('--r-resol', type=int, default=1, help='Starting resolution for SO3 grid (default: %(default)s')
    parser.add_argument('--keep-r', type=int, default=30, help='Rotations to keep per iteration')
    parser.add_argument('--keep-t', type=int, default=30, help='Translations to keep per iteration')
    parser.add_argument('--flip', action='store_true', help='Flip handedness of volume')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')
    return parser

class VolumeAligner:
    def __init__(self, vol, vol_ref=None, maxD=128, flip=False):
        if flip: vol = vol[::-1]
        nz, ny, nx = vol.shape
        assert nz==ny==nx, 'Volume must be cubic'
        D = nz
        self.vol_real = torch.from_numpy(vol.astype(np.float32)).view(1,1,D,D,D)

        vol = fft.fftn_center(vol)
        # fourier clip the volume
        if D > maxD: 
            assert maxD % 2 == 0
            b, e = D//2 - maxD//2, D//2 + maxD//2
            vol = vol[b:e,b:e,b:e]
            D = maxD

        x2, x1, x0 = np.meshgrid(np.linspace(-1, 1, D, endpoint=True), 
                             np.linspace(-1, 1, D, endpoint=True),
                             np.linspace(-1, 1, D, endpoint=True),
                             indexing='ij')
        lattice = np.stack([x0.ravel(), x1.ravel(), x2.ravel()],1).astype(np.float32)
        self.lattice = torch.from_numpy(lattice)

        x2, x1, x0 = np.meshgrid(np.linspace(-.5, .5, D, endpoint=False), 
                             np.linspace(-.5, .5, D, endpoint=False),
                             np.linspace(-.5, .5, D, endpoint=False),
                             indexing='ij')
        #tcoords = np.stack([x0.ravel(), x1.ravel(), x2.ravel()],1).astype(np.float32)
        tcoords = np.stack([x0,x1,x2],-1).astype(np.float32)
        self.tcoords = torch.from_numpy(tcoords)

        self.volr = torch.from_numpy(vol.real.astype(np.float32)).view(1,1,D,D,D)
        self.voli = torch.from_numpy(vol.imag.astype(np.float32)).view(1,1,D,D,D)

        if vol_ref is not None:
            vol_ref = fft.fftn_center(vol_ref)
            if vol_ref.shape[0] > D:
                b, e = vol_ref.shape[0]//2 - D//2, vol_ref.shape[0]//2 + D//2
                vol_ref = vol_ref[b:e,b:e,b:e]
            self.vol_refr = torch.from_numpy(vol_ref.real.astype(np.float32)).view(1,1,D,D,D)
            self.vol_refi = torch.from_numpy(vol_ref.imag.astype(np.float32)).view(1,1,D,D,D)

        # FT is not symmetric around origin
        c = 2/(D-1)*(D/2) -1 
        self.center = torch.tensor([c,c,c]) # pixel coordinate for vol[D/2,D/2,D/2]
        self.D = D
        self.D2 = D // 2

    def get_mask(self, r):
        raise NotImplementedError
        
    def use_cuda(self):
        self.lattice = self.lattice.cuda()
        self.volr = self.volr.cuda()
        self.voli = self.voli.cuda()
        self.vol_refr = self.vol_refr.cuda()
        self.vol_refi = self.vol_refi.cuda()
        self.vol_real = self.vol_real.cuda()
        self.tcoords = self.tcoords.cuda()

    def real_tform(self, rot, trans):
        if rot.shape[-1] == 4:
            rot = lie_tools.quaternions_to_SO3(rot)
        B = rot.size(0)
        D = self.D
        grid = self.lattice @ rot # B x D^3 x 3 
        grid = grid.view(-1, D, D, D, 3)
        offset = self.center - grid[:,self.D2, self.D2, self.D2]
        grid += offset[:,None,None,None,:]
        grid -= trans 
        grid = grid.view(1, -1, D, D, 3)
        vol = F.grid_sample(self.vol_real, grid)
        return vol.view(B,D,D,D)

    def rotate(self, rot_or_quat):
        rot = rot_or_quat
        if rot_or_quat.shape[-1] == 4:
            rot = lie_tools.quaternions_to_SO3(rot_or_quat)
        B = rot.size(0)
        D = self.D
        grid = self.lattice @ rot # B x D^3 x 3 
        grid = grid.view(-1, D, D, D, 3)
        offset = self.center - grid[:,self.D2, self.D2, self.D2]
        grid += offset[:,None,None,None,:]
        grid = grid.view(1, -1, D, D, 3)
        volr = F.grid_sample(self.volr, grid)
        volr = volr.view(B, D, D, D)
        voli = F.grid_sample(self.voli, grid)
        voli = voli.view(B, D, D, D)
        return volr, voli

    def translate(self, volr, voli, t):
        '''
        Translate a volume by phase shifting its Fourier transform
        
        Inputs:
            #coords: wavevectors between [-.5,.5] (img_dims x 2)
            
            img: FT of image (B x img_dims x 2)
            t: shift in pixels (B x T x 3)

        Returns:
            Shifted images (B x T x img_dims x 2) 

        img_dims can either be 2D or 1D (unraveled image) 
        '''
        # F'(k) = exp(-2*pi*k*x0)*F(k)
        D = self.D
        B = t.size(0)
        T = t.size(1)
        t = t.view(B,T,1,1,3,1) # To be able to do bmm
        tfilt = self.tcoords @ t * -2 * np.pi # BxTxNx1
        tfilt = tfilt.squeeze(-1) # BxTxN
        c = torch.cos(tfilt) # BxTxN
        s = torch.sin(tfilt) # BxTxN
        volr = volr.unsqueeze(1) # Bx1xN
        voli = voli.unsqueeze(1) # Bx1xN
        vr = volr*c-voli*s
        vi = volr*s+voli*c 
        return vr.view(B,T,D,D,D), vi.view(B,T,D,D,D)

    def compute_err_old(self, vr, vi):
        err_r = (self.vol_refr - vr).pow(2).sum((-1,-2,-3))
        err_i = (self.vol_refi - vi).pow(2).sum((-1,-2,-3))
        return err_r + err_i

    def compute_err(self, vr, vi):
        r, i = self.vol_refr, self.vol_refi
        corr = (r*vr+i*vi).sum((-1,-2,-3)) / ((vr*vr+vi*vi).sum((-1,-2,-3))*(r*r+i*i).sum((-1,-2,-3)))**.5
        return -corr

class GridPose(data.Dataset):
    def __init__(self, pose, pose_id):
        assert len(pose) == len(pose_id)
        self.pose = pose
        self.pose_id = pose_id
        self.N = len(pose)
    def __len__(self):
        return self.N
    def __getitem__(self, index):
        return self.pose[index], self.pose_id[index]

class MinPoseTracker():
    def __init__(self, max_keep, pose_dim, pose_id_dim):
        self.min_errs = torch.ones(max_keep)*np.inf
        self.min_poses = torch.empty((max_keep, pose_dim))
        self.min_pose_ids = torch.empty((max_keep, pose_id_dim), device='cpu', dtype=torch.int64)
        self.max_keep = max_keep
        self.pose_dim = pose_dim
        self.pose_id_dim = pose_id_dim

    def add(self, errs, poses, pose_ids):
        e = torch.cat((self.min_errs, errs))
        a = torch.argsort(e)[:self.max_keep]
        self.min_errs = e[a]
        self.min_poses = torch.cat((self.min_poses,poses))[a]
        self.min_pose_ids = torch.cat((self.min_pose_ids, pose_ids))[a]

    def clear(self):
        self.min_errs = torch.ones(self.max_keep)*np.inf
        self.min_poses = torch.empty((self.max_keep, self.pose_dim))
        self.min_pose_ids = torch.empty((self.max_keep, self.pose_id_dim), device='cpu', dtype=torch.int64)

def subdivide_r(quat, grid_id, curr_res):
    neighbors = [so3_grid.get_neighbor(quat[i], grid_id[i][0], grid_id[i][1], curr_res) for i in range(len(quat))]
    quat_next = np.array([x[0] for x in neighbors]).reshape(-1,4) # Qx8x4->8Qx4
    grid_id_next = np.array([x[1] for x in neighbors]).reshape(-1,2) # Qx8x2->8Qx2
    return quat_next, grid_id_next

def subdivide_t(t_ids, resol, extent, ngrid):
    next_ = [shift_grid3.get_neighbor(xi, yi, zi, resol, extent, ngrid) for xi,yi,zi in t_ids]
    trans_next = np.asarray([x[0] for x in next_]).reshape(-1,3)
    t_ids_next = np.asarray([x[1] for x in next_]).reshape(-1,3)
    return trans_next, t_ids_next

def main(args):
    log(args)
    torch.set_grad_enabled(False) 
    use_cuda = torch.cuda.is_available()
    log('Use cuda {}'.format(use_cuda))
    if use_cuda:
        torch.set_default_tensor_type(torch.cuda.FloatTensor)

    t1 = time.time()    
    ref, _ = mrc.parse_mrc(args.ref)
    log('Loaded {} volume'.format(ref.shape))
    vol, _ = mrc.parse_mrc(args.vol)
    log('Loaded {} volume'.format(vol.shape))

    projector = VolumeAligner(vol, vol_ref=ref, maxD=args.max_D, flip=args.flip)
    if use_cuda:
        projector.use_cuda()

    r_resol = args.r_resol
    quats = so3_grid.grid_SO3(r_resol)
    q_id = np.arange(len(quats))
    q_id = np.stack([q_id // (6*2**r_resol), q_id % (6*2**r_resol)], -1)
    rots = GridPose(quats, q_id)

    t_resol = 0
    T_EXTENT = vol.shape[0]/16 if args.t_extent is None else args.t_extent
    T_NGRID = args.t_grid
    trans = shift_grid3.base_shift_grid(T_EXTENT, T_NGRID)
    t_id = np.stack(shift_grid3.get_base_id(np.arange(len(trans)), T_NGRID),-1)
    trans = GridPose(trans, t_id)
    
    max_keep_r = args.keep_r
    max_keep_t = args.keep_t
    #rot_tracker = MinPoseTracker(max_keep_r, 4, 2)
    #tr_tracker = MinPoseTracker(max_keep_t, 3, 3)
    for it in range(args.niter):
        log('Iteration {}'.format(it))
        log('Generating {} rotations'.format(len(rots)))
        log('Generating {} translations'.format(len(trans)))
        pose_err = np.empty((len(rots),len(trans)),dtype=np.float32)
        #rot_tracker.clear()
        #tr_tracker.clear()
        r_iterator = data.DataLoader(rots, batch_size=args.rb, shuffle=False)
        t_iterator = data.DataLoader(trans, batch_size=args.tb, shuffle=False)
        r_it = 0
        for rot, r_id in r_iterator:
            if use_cuda: rot = rot.cuda()
            vr, vi = projector.rotate(rot)
            t_it = 0
            for tr, t_id in t_iterator:
                if use_cuda: tr = tr.cuda()
                vtr, vti = projector.translate(vr, vi, tr.expand(rot.size(0), *tr.shape)) 
                # todo: check volume
                err = projector.compute_err(vtr, vti) # R x T
                pose_err[r_it:r_it+len(rot), t_it:t_it+len(tr)] = err.cpu().numpy()
                #r_err = err.min(1)[0]
                #min_r_err, min_r_i = r_err.sort()
                #rot_tracker.add(min_r_err[:max_keep_r], rot[min_r_i][:max_keep_r], r_id[min_r_i][:max_keep_r])
                #t_err= err.min(0)[0]
                #min_t_err, min_t_i = t_err.sort()
                #tr_tracker.add(min_t_err[:max_keep_t], tr[min_t_i][:max_keep_t], t_id[min_t_i][:max_keep_t])
                t_it += len(tr)
            r_it += len(rot)
        
        r_err = pose_err.min(1)
        r_err_argmin = r_err.argsort()[:max_keep_r]
        t_err = pose_err.min(0)
        t_err_argmin = t_err.argsort()[:max_keep_t]

        # lstart
        #r = rots.pose[r_err_argmin[0]]
        #t = trans.pose[t_err_argmin[0]]
        #log('Best rot: {}'.format(r))
        #log('Best trans: {}'.format(t))
        #vr, vi = projector_full.rotate(torch.tensor(r).unsqueeze(0))
        #vr, vi = projector_full.translate(vr, vi, torch.tensor(t).view(1,1,3))
        #err = projector_full.compute_err(vr,vi)

        #w = np.where(r_err[r_err_argmin] > err.item())[0]
        rots, rots_id = subdivide_r(rots.pose[r_err_argmin], rots.pose_id[r_err_argmin], r_resol)
        rots = GridPose(rots, rots_id)

        t_err = pose_err.min(0)
        t_err_argmin = t_err.argsort()[:max_keep_t]
        trans, trans_id = subdivide_t(trans.pose_id[t_err_argmin], t_resol, T_EXTENT, T_NGRID)
        trans = GridPose(trans, trans_id)
        r_resol += 1
        t_resol += 1
        vlog(r_err[r_err_argmin])
        vlog(t_err[t_err_argmin])
        #log(rot_tracker.min_errs)
        #log(tr_tracker.min_errs)
    r = rots.pose[r_err_argmin[0]]
    t = trans.pose[t_err_argmin[0]]*vol.shape[0]/args.max_D
    log('Best rot: {}'.format(r))
    log('Best trans: {}'.format(t))
    t *= 2/vol.shape[0]
    projector = VolumeAligner(vol, vol_ref=ref, maxD=vol.shape[0], flip=args.flip)
    if use_cuda: projector.use_cuda()
    vr = projector.real_tform(torch.tensor(r).unsqueeze(0),torch.tensor(t).view(1,1,3))
    v = vr.squeeze().cpu().numpy()
    log('Saving {}'.format(args.o))
    mrc.write(args.o, v.astype(np.float32))
    
    td = time.time()-t1
    log('Finished in {}s'.format(td))

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)
