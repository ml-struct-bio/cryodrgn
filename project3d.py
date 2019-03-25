'''
Generate projections of a 3D volume
'''

import argparse
import numpy as np
import sys, os
import time
import pickle

import torch
import torch.nn as nn
import torch.nn.functional as F

sys.path.insert(0,'{}/lib-python'.format(os.path.dirname(os.path.abspath(__file__))))
import utils
import mrc
import fft
import lie_tools

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

log = utils.log

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mrc', help='Input volume')
    parser.add_argument('-o', type=os.path.abspath, required=True, help='Output projection stack (.mrcs)')
    parser.add_argument('--out-rot', type=os.path.abspath, required=True, help='Output rotations (.pkl)')
    parser.add_argument('--out-png', type=os.path.abspath, help='Montage of first 9 projections')
    parser.add_argument('-N', type=int, required=True, help='Number of random projections')
    parser.add_argument('-b', type=int, default=100, help='Minibatch size (default: %(default)s)')
    parser.add_argument('--grid', type=int, help='Generate projections on a uniform deterministic grid on SO3. Specify resolution level')
    parser.add_argument('--tilt', type=float, help='Right-handed x-axis tilt offset in degrees')
    parser.add_argument('--seed', type=int, help='Random seed')
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

    def project(self, rot):
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
        return vol.sum(dim=1)
   
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

def main(args):
    for out in (args.o, args.out_rot, args.out_png):
        if not out: continue
        mkbasedir(out)
        warnexists(out)

    if args.seed is not None:
        torch.manual_seed(args.seed)

    use_cuda = torch.cuda.is_available()
    log('Use cuda {}'.format(use_cuda))
    if use_cuda:
        torch.set_default_tensor_type(torch.cuda.FloatTensor)

    t1 = time.time()    
    vol, _ , _ = mrc.parse_mrc(args.mrc)
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

    if args.grid:
        raise NotImplementedError
    
    rots = []
    imgs = []
    for mb in range(int(args.N/args.b)):
        log('Projecting {}/{}'.format((mb+1)*args.b, args.N))
        rot = lie_tools.random_SO3(args.b)
        projections = projector.project(rot)
        projections = projections.cpu().numpy()
        rots.append(rot)
        imgs.append(projections)
    if args.N % args.b:
        log('Projecting {}/{}'.format(args.N, args.N))
        rot = lie_tools.random_SO3(args.N % args.b)
        projections = projector.project(rot)
        projections = projections.cpu().numpy()
        rots.append(rot)
        imgs.append(projections)

    rots = np.vstack(rots)
    imgs = np.vstack(imgs)
    td = time.time()-t1
    log('Projected {} images in {}s ({}s per image)'.format(args.N, td, td/args.N ))

    log('Saving {}'.format(args.o))
    mrc.write(args.o,imgs.astype(np.float32))
    log('Saving {}'.format(args.out_rot))
    with open(args.out_rot,'wb') as f:
        pickle.dump(rots, f)
    if args.out_png:
        log('Saving {}'.format(args.out_png))
        plot_projections(args.out_png, imgs[:9])

if __name__ == '__main__':
    main(parse_args().parse_args())
