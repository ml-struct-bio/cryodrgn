'''
'''

import argparse
import numpy as np
import sys, os
import pickle
import subprocess

import matplotlib.pyplot as plt
import seaborn as sns

import cryodrgn
from cryodrgn import analysis
from cryodrgn import utils

log = utils.log

def add_args(parser):
    parser.add_argument('workdir', type=os.path.abspath, help='')
    parser.add_argument('epoch', type=int, help='')
    parser.add_argument('--device', type=int, help='CUDA device')

    group = parser.add_argument_group('Extra arguments for volume generation')
    group.add_argument('--Apix', type=float, default=1, help='Pixel size to add to .mrc header (default: %(default)s A/pix)')
    group.add_argument('--flip', action='store_true', help='Flip handedness of output volume')
    group.add_argument('-d','--downsample', type=int, help='Downsample volumes to this box size (pixels)')
    return parser

def analyze_z1(z, outdir, vg):
    '''Plotting and volume generation for 1D z'''
    assert z.shape[1] == 1
    z = z.reshape(-1)
    N = len(z)

    plt.figure(1)
    plt.scatter(np.arange(N), z, alpha=.1, s=2)
    plt.xlabel('particle')
    plt.ylabel('z')
    plt.savefig(f'{outdir}/z.png')

    plt.figure(2)
    sns.distplot(z)
    plt.xlabel('z')
    plt.savefig(f'{outdir}/z_hist.png')

    ztraj = np.linspace(*np.percentile(z,(5,95)), 10) # or np.percentile(z, np.linspace(5,95,10)) ?
    vg.gen_volumes(outdir, ztraj)

def analyze_zN(z, outdir, vg):
    zdim = z.shape[1]

    # Principal component analysis
    pc, pca = analysis.run_pca(z)
    start, end = np.percentile(pc[:,0],(5,95))
    z_pc1 = analysis.get_pc_traj(pca, z.shape[1], 10, 1, start, end)
    start, end = np.percentile(pc[:,1],(5,95))
    z_pc2 = analysis.get_pc_traj(pca, z.shape[1], 10, 2, start, end)

    # kmeans clustering
    kmeans_labels, centers = analysis.cluster_kmeans(z, 20)
    centers, centers_ind = analysis.get_nearest_point(z, centers)

    vg.gen_volumes('f{outdir}/pc1', z_pc1)
    vg.gen_volumes('f{outdir}/pc2', z_pc2)
    vg.gen_volumes('f{outdir}/kmeans20', centers)

    # UMAP
    if zdim > 2:
        umap_emb = analysis.run_umap(z)

    # Make some plots -- todo

    
 
class VolumeGenerator:
    '''Helper class to call analysis.gen_volumes'''
    def __init__(self, weights, config, vol_args={}):
        self.weights = weights
        self.config = config
        self.vol_args = vol_args

    def gen_volumes(self, outdir, z_values):
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        zfile = np.savetxt(f'{outdir}/z_values.txt', z_values)
        analysis.gen_volumes(self.weights, self.config, zfile, outdir, **self.vol_args)

def main(args):
    E = args.epoch
    workdir = args.workdir
    zfile = f'{workdir}/z.{E}.pkl'
    weights = f'{workdir}/weights.{E}.pkl'
    config = f'{workdir}/config.pkl'
    
    outdir = f'{workdir}/analyze.{E}'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    z = utils.load_pkl(zfile)
    zdim = z.shape[1]

    vol_args = dict(Apix=args.Apix, downsample=args.downsample, flip=args.flip, cuda=args.cuda)
    vg = VolumeGenerator(weights, config, vol_args)

    if zdim == 1:
        analyze_z1(z, outdir, vg)
    else:
        analyze_zN(z, outdir, vg)
       
    # copy over template if file doesn't exist
    if not os.path.exists(f'{outdir}/cryoDRGN_viz_template.ipynb'):
        ipynb = f'{cryodrgn._ROOT}/templates/cryoDRGN_viz_template.ipynb'
        cmd = f'cp {ipynb} {outdir}'
        subprocess.check_call(cmd, shell=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
