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
    parser.add_argument('--seed', type=int, help='Random seed')
    return parser

class ProjectorFT:
    def __init__(self, vol):
        nz, ny, nx = vol.shape
        assert nz==ny==nx, 'Volume must be cubic'
        x0, x1 = np.meshgrid(np.linspace(-1, 1, nx, endpoint=True), 
                             np.linspace(-1, 1, ny, endpoint=True))
        x2 = np.zeros(ny*nx)
        lattice = np.stack([x0.ravel(), x1.ravel(), x2],1).astype(np.float32)
        self.lattice = torch.from_numpy(lattice)

        self.vol_ft = torch.from_numpy(fft.htn_center(vol).astype(np.float32))
        self.vol_ft = self.vol_ft.unsqueeze(0)
        self.vol_ft = self.vol_ft.unsqueeze(0)

        self.ny = ny
        self.nx = nx

        # FT is not symmetric around origin
        D = nz
        c = 2/(D-1)*(D/2) -1 
        self.center = torch.tensor([c,c,c]) # pixel coordinate for vol[D/2,D/2,D/2]

    def project(self, rot):
        grid = self.lattice @ rot # rot.T * coord
        grid = grid.view(1, -1, self.ny, self.nx, 3)
        offset = self.center - grid[0,:,int(self.ny/2),int(self.nx/2)]
        grid += offset[None,:,None,None,:]
        return F.grid_sample(self.vol_ft, grid)
    
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

    t1 = time.time()    
    vol, _ , _ = mrc.parse_mrc(args.mrc)
    log('Loaded {} volume'.format(vol.shape))

    projector = ProjectorFT(vol)

    if args.grid:
        raise NotImplementedError
    
    rots = []
    imgs = []
    for mb in range(int(args.N/args.b)):
        log('Projecting {}/{}'.format((mb+1)*args.b, args.N))
        rot = lie_tools.random_SO3(args.b)
        test = projector.project(rot)
        test = test.squeeze().numpy()
        rots.append(rot)
        imgs.append(test)
    if args.N % args.b:
        log('Projecting {}/{}'.format(args.N, args.N))
        rot = lie_tools.random_SO3(args.N % args.b)
        test = projector.project(rot)
        test = test.squeeze().numpy()
        rots.append(rot)
        imgs.append(test)

    rots = np.vstack(rots)
    imgs = np.vstack(imgs)
    imgs = np.asarray([fft.ihtn_center(x) for x in imgs])
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
