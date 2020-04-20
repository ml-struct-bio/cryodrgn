'''
Generate projections of a 3D volume
'''

import argparse
import numpy as np
import sys, os
import time
import pickle
from scipy.ndimage.fourier import fourier_shift

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.utils.data as data

from cryodrgn import utils
from cryodrgn import mrc
from cryodrgn import fft
from cryodrgn import lie_tools
from cryodrgn import so3_grid

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

log = utils.log
vlog = utils.vlog

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mrc', help='Input volume')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output projection stack (.mrcs)')
    parser.add_argument('--out-pose', type=os.path.abspath, required=True, help='Output poses (.pkl)')
    parser.add_argument('--out-png', type=os.path.abspath, help='Montage of first 9 projections')
    parser.add_argument('-N', type=int, help='Number of random projections')
    parser.add_argument('-b', type=int, default=100, help='Minibatch size (default: %(default)s)')
    parser.add_argument('--t-extent', type=float, default=5, help='Extent of image translation in pixels (default: +/-%(default)s)')
    parser.add_argument('--grid', type=int, help='Generate projections on a uniform deterministic grid on SO3. Specify resolution level')
    parser.add_argument('--tilt', type=float, help='Right-handed x-axis tilt offset in degrees')
    parser.add_argument('--seed', type=int, help='Random seed')
    parser.add_argument('-v','--verbose',action='store_true',help='Increaes verbosity')
    return parser

class Projector:
    def __init__(self, vol, tilt=None):
        nz, ny, nx = vol.shape
        assert nz==ny==nx, 'Volume must be cubic'
        x2, x1, x0 = np.meshgrid(np.linspace(-1, 1, nz, endpoint=True), 
                             np.linspace(-1, 1, ny, endpoint=True),
                             np.linspace(-1, 1, nx, endpoint=True),
                             indexing='ij')

        lattice = np.stack([x0.ravel(), x1.ravel(), x2.ravel()],1).astype(np.float32)
        self.lattice = torch.from_numpy(lattice)

        self.vol = torch.from_numpy(vol.astype(np.float32))
        self.vol = self.vol.unsqueeze(0)
        self.vol = self.vol.unsqueeze(0)

        self.nz = nz
        self.ny = ny
        self.nx = nx

        # FT is not symmetric around origin
        D = nz
        c = 2/(D-1)*(D/2) -1 
        self.center = torch.tensor([c,c,c]) # pixel coordinate for vol[D/2,D/2,D/2]

        if tilt is not None:
            assert tilt.shape == (3,3)
            tilt = torch.tensor(tilt)
        self.tilt = tilt

    def rotate(self, rot):
        B = rot.size(0)
        if self.tilt is not None:
            rot = self.tilt @ rot
        grid = self.lattice @ rot # B x D^3 x 3 
        grid = grid.view(-1, self.nz, self.ny, self.nx, 3)
        offset = self.center - grid[:,int(self.nz/2),int(self.ny/2),int(self.nx/2)]
        grid += offset[:,None,None,None,:]
        grid = grid.view(1, -1, self.ny, self.nx, 3)
        vol = F.grid_sample(self.vol, grid)
        vol = vol.view(B,self.nz,self.ny,self.nx)
        return vol

    def project(self, rot):
        return self.rotate(rot).sum(dim=1)
   
class RandomRot(data.Dataset):
    def __init__(self, N):
        self.N = N
        self.rots = lie_tools.random_SO3(N)
    def __len__(self):
        return self.N
    def __getitem__(self, index):
        return self.rots[index]

class GridRot(data.Dataset):
    def __init__(self, resol):
        quats = so3_grid.grid_SO3(resol)
        self.rots = lie_tools.quaternions_to_SO3(torch.tensor(quats))
        self.N = len(self.rots)
    def __len__(self):
        return self.N
    def __getitem__(self, index):
        return self.rots[index]

def plot_projections(out_png, imgs):
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(10,10))
    axes = axes.ravel()
    for i in range(min(len(imgs),9)):
        axes[i].imshow(imgs[i])
    plt.savefig(out_png)

def mkbasedir(out):
    if not os.path.exists(os.path.dirname(out)):
        os.makedirs(os.path.dirname(out))

def warnexists(out):
    if os.path.exists(out):
        log('Warning: {} already exists. Overwriting.'.format(out))

def translate_img(img, t):
    '''
    img: BxYxX real space image
    t: Bx2 shift in pixels
    '''
    ff = np.fft.fft2(np.fft.fftshift(img))
    ff = fourier_shift(ff, t)
    return np.fft.fftshift(np.fft.ifft2(ff)).real

def main(args):
    for out in (args.o, args.out_png, args.out_pose):
        if not out: continue
        mkbasedir(out)
        warnexists(out)

    if args.t_extent == 0.:
        log('Not shifting images')
    else:
        assert args.t_extent > 0

    if args.seed is not None:
        np.random.seed(args.seed)
        torch.manual_seed(args.seed)

    use_cuda = torch.cuda.is_available()
    log('Use cuda {}'.format(use_cuda))
    if use_cuda:
        torch.set_default_tensor_type(torch.cuda.FloatTensor)

    t1 = time.time()    
    vol, _ = mrc.parse_mrc(args.mrc)
    log('Loaded {} volume'.format(vol.shape))

    if args.tilt:
        theta = args.tilt*np.pi/180
        args.tilt = np.array([[1.,0.,0.],
                        [0, np.cos(theta), -np.sin(theta)],
                        [0, np.sin(theta), np.cos(theta)]]).astype(np.float32)

    projector = Projector(vol, args.tilt)
    if use_cuda:
        projector.lattice = projector.lattice.cuda()
        projector.vol = projector.vol.cuda()

    if args.grid is not None:
        rots = GridRot(args.grid)
        log('Generating {} rotations at resolution level {}'.format(len(rots), args.grid))
    else:
        log('Generating {} random rotations'.format(args.N))
        rots = RandomRot(args.N)
    
    log('Projecting...')
    imgs = []
    iterator = data.DataLoader(rots, batch_size=args.b)
    for i, rot in enumerate(iterator):
        vlog('Projecting {}/{}'.format((i+1)*len(rot), args.N))
        projections = projector.project(rot)
        projections = projections.cpu().numpy()
        imgs.append(projections)

    rots = rots.rots.cpu().numpy()
    imgs = np.vstack(imgs)
    td = time.time()-t1
    log('Projected {} images in {}s ({}s per image)'.format(args.N, td, td/args.N ))

    if args.t_extent:
        log('Shifting images between +/- {} pixels'.format(args.t_extent))
        trans = np.random.rand(args.N,2)*2*args.t_extent - args.t_extent
        imgs = np.asarray([translate_img(img, t) for img,t in zip(imgs,trans)])
        # convention: we want the first column to be x shift and second column to be y shift
        # reverse columns since current implementation of translate_img uses scipy's 
        # fourier_shift, which is flipped the other way
        # convention: save the translation that centers the image
        trans = -trans[:,::-1]
        # convert translation from pixel to fraction
        D = imgs.shape[-1]
        assert D % 2 == 0
        trans /= D

    log('Saving {}'.format(args.o))
    mrc.write(args.o,imgs.astype(np.float32))
    log('Saving {}'.format(args.out_pose))
    with open(args.out_pose,'wb') as f:
        if args.t_extent:
            pickle.dump((rots,trans),f)
        else:
            pickle.dump(rots, f)
    if args.out_png:
        log('Saving {}'.format(args.out_png))
        plot_projections(args.out_png, imgs[:9])

if __name__ == '__main__':
    args = parse_args().parse_args()
    utils._verbose = args.verbose
    main(args)
