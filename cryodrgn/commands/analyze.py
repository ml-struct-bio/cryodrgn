'''
Visualize latent space and generate volumes
'''

import argparse
import numpy as np
import sys, os
import pickle
import subprocess
from datetime import datetime as dt

import matplotlib
matplotlib.use('Agg') # non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns

import cryodrgn
from cryodrgn import analysis
from cryodrgn import utils

log = utils.log

def add_args(parser):
    parser.add_argument('workdir', type=os.path.abspath, help='Directory with cryoDRGN results')
    parser.add_argument('epoch', type=int, help='Epoch number N to analyze (0-based indexing, corresponding to z.N.pkl, weights.N.pkl)')
    parser.add_argument('--device', type=int, help='Optionally specify CUDA device')
    parser.add_argument('-o','--outdir', help='Output directory for analysis results (default: [workdir]/analyze.[epoch])')
    parser.add_argument('--skip-vol', action='store_true', help='Skip generation of volumes')
    parser.add_argument('--skip-umap', action='store_true', help='Skip running UMAP')

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

def analyze_zN(z, outdir, vg, skip_umap=False):
    zdim = z.shape[1]

    # Principal component analysis
    log('Perfoming principal component analysis...')
    pc, pca = analysis.run_pca(z)
    start, end = np.percentile(pc[:,0],(5,95))
    z_pc1 = analysis.get_pc_traj(pca, z.shape[1], 10, 1, start, end)
    start, end = np.percentile(pc[:,1],(5,95))
    z_pc2 = analysis.get_pc_traj(pca, z.shape[1], 10, 2, start, end)

    # kmeans clustering
    log('K-means clustering...')
    K = 20
    kmeans_labels, centers = analysis.cluster_kmeans(z, K)
    centers, centers_ind = analysis.get_nearest_point(z, centers)
    if not os.path.exists(f'{outdir}/kmeans20'): 
        os.mkdir(f'{outdir}/kmeans20')
    utils.save_pkl(kmeans_labels, f'{outdir}/kmeans20/labels.pkl')
    np.savetxt(f'{outdir}/kmeans20/centers.txt', centers)
    np.savetxt(f'{outdir}/kmeans20/centers_ind.txt', centers_ind, fmt='%d')

    # Generate volumes
    log('Generating volumes...')
    vg.gen_volumes(f'{outdir}/pc1', z_pc1)
    vg.gen_volumes(f'{outdir}/pc2', z_pc2)
    vg.gen_volumes(f'{outdir}/kmeans20', centers)

    # UMAP -- slow step
    if zdim > 2 and not skip_umap:
        log('Running UMAP...')
        umap_emb = analysis.run_umap(z)
        utils.save_pkl(umap_emb, f'{outdir}/umap.pkl')

    # Make some plots
    log('Generating plots...')
    plt.figure(1)
    plt.scatter(pc[:,0], pc[:,1], alpha=.1, s=2)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.savefig(f'{outdir}/z_pca.png')
    
    if zdim > 2 and not skip_umap:
        plt.figure(2)
        plt.scatter(umap_emb[:,0], umap_emb[:,1], alpha=.1, s=2)
        plt.xlabel('UMAP1')
        plt.ylabel('UMAP2')
        plt.savefig(f'{outdir}/umap.png')

    analysis.plot_by_cluster(pc[:,0], pc[:,1], K, kmeans_labels, centers_ind=centers_ind, annotate=True)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.savefig(f'{outdir}/kmeans20/z_pca.png')

    if zdim > 2 and not skip_umap:
        analysis.plot_by_cluster(umap_emb[:,0], umap_emb[:,1], K, kmeans_labels, centers_ind=centers_ind, annotate=True)
        plt.xlabel('UMAP1')
        plt.ylabel('UMAP2')
        plt.savefig(f'{outdir}/kmeans20/umap.png')
 
class VolumeGenerator:
    '''Helper class to call analysis.gen_volumes'''
    def __init__(self, weights, config, vol_args={}, skip_vol=False):
        self.weights = weights
        self.config = config
        self.vol_args = vol_args
        self.skip_vol = skip_vol

    def gen_volumes(self, outdir, z_values):
        if self.skip_vol: return 
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        zfile = f'{outdir}/z_values.txt'
        np.savetxt(zfile, z_values)
        analysis.gen_volumes(self.weights, self.config, zfile, outdir, **self.vol_args)

def main(args):
    t1 = dt.now()
    E = args.epoch
    workdir = args.workdir
    zfile = f'{workdir}/z.{E}.pkl'
    weights = f'{workdir}/weights.{E}.pkl'
    config = f'{workdir}/config.pkl'
    outdir = f'{workdir}/analyze.{E}'
    if E == -1:
        zfile = f'{workdir}/z.pkl'
        weights = f'{workdir}/weights.pkl'
        outdir = f'{workdir}/analyze'
    
    if args.outdir:
        outdir = args.outdir
    log(f'Saving results to {outdir}')
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    z = utils.load_pkl(zfile)
    zdim = z.shape[1]

    vol_args = dict(Apix=args.Apix, downsample=args.downsample, flip=args.flip, cuda=args.device)
    vg = VolumeGenerator(weights, config, vol_args, skip_vol=args.skip_vol)

    if zdim == 1:
        analyze_z1(z, outdir, vg)
    else:
        analyze_zN(z, outdir, vg, args.skip_umap)
       
    # copy over template if file doesn't exist
    out_ipynb = f'{outdir}/cryoDRGN_viz.ipynb'
    if not os.path.exists(out_ipynb):
        log(f'Creating jupyter notebook...')
        ipynb = f'{cryodrgn._ROOT}/templates/cryoDRGN_viz_template.ipynb'
        cmd = f'cp {ipynb} {out_ipynb}'
        subprocess.check_call(cmd, shell=True)
    else:
        log(f'{out_ipynb} already exists. Skipping')
    log(out_ipynb)
    
    log(f'Finished in {dt.now()-t1}')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
