'''
Pipeline to analyze cryoDRGN volume distribution
'''

import argparse
import numpy as np
import sys, os
import pickle
import shutil
from datetime import datetime as dt
from collections import Counter

from cryodrgn import mrc
from sklearn.decomposition import PCA

from scipy.ndimage.morphology import binary_dilation
from cryodrgn import analysis 
from cryodrgn import utils

import matplotlib
matplotlib.use('Agg')
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.cluster import AgglomerativeClustering

log = utils.log 

def add_args(parser):
    parser.add_argument('workdir', type=os.path.abspath, help='Directory with cryoDRGN results')
    parser.add_argument('epoch', type=int, help='Epoch number N to analyze (0-based indexing, corresponding to z.N.pkl, weights.N.pkl)')
    parser.add_argument('--device', type=int, help='Optionally specify CUDA device')
    parser.add_argument('-o','--outdir', type=os.path.abspath, help='Output directory for landscape analysis results (default: [workdir]/landscape.[epoch])')
    parser.add_argument('--skip-umap', action='store_true', help='Skip running UMAP')
    parser.add_argument('--vol-ind', type=os.path.abspath, help='Index .pkl for filtering volumes')

    group = parser.add_argument_group('Extra arguments for volume generation')
    group.add_argument('-N','--sketch-size', type=int, default=1000, help='Number of volumes to generate for analysis (default: %(default)s)')
    group.add_argument('--Apix', type=float, default=1, help='Pixel size to add to .mrc header (default: %(default)s A/pix)')
    group.add_argument('--flip', action='store_true', help='Flip handedness of output volume')
    group.add_argument('-d','--downsample', type=int, default=128, help='Downsample volumes to this box size (pixels) (default: %(default)s)')
    group.add_argument('--skip-vol', action='store_true', help='Skip generation of volumes')

    group = parser.add_argument_group('Extra arguments for mask generation')
    group.add_argument('--thresh', type=float, help='Density value to threshold for masking (default: half of max density value)')
    group.add_argument('--dilate', type=int, default=5, help='Dilate initial mask by this amount (default: %(default)s pixels)')
    group.add_argument('--mask', metavar='MRC', type=os.path.abspath, help='Path to a custom mask. Must be same box size as generated volumes.')

    group = parser.add_argument_group('Extra arguments for clustering')
    group.add_argument('--linkage', default='average', help='Linkage for agglomerative clustering (e.g. average, ward) (default: %(default)s)')
    group.add_argument('-M', type=int, default=10, help='Number of clusters (default: %(default)s)')

    group = parser.add_argument_group('Extra arguments for landscape visualization')
    group.add_argument('--pc-dim', default=20, help='PCA dimensionality reduction (default: %(default)s)')
    group.add_argument('--plot-dim', type=int, default=5, help='Number of dimensions to plot (default: %(default)s)')

    return parser

def generate_volumes(z, outdir, vg, K):
    # kmeans clustering
    log('Sketching distribution...')
    kmeans_labels, centers = analysis.cluster_kmeans(z, K, on_data=True, reorder=True)
    centers, centers_ind = analysis.get_nearest_point(z, centers)
    if not os.path.exists(f'{outdir}/kmeans{K}'): 
        os.mkdir(f'{outdir}/kmeans{K}')
    utils.save_pkl(kmeans_labels, f'{outdir}/kmeans{K}/labels.pkl')
    np.savetxt(f'{outdir}/kmeans{K}/centers.txt', centers)
    np.savetxt(f'{outdir}/kmeans{K}/centers_ind.txt', centers_ind, fmt='%d')
    log('Generating volumes...')
    vg.gen_volumes(f'{outdir}/kmeans{K}', centers)

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


def make_mask(outdir, K, dilate, thresh, in_mrc=None):
    if in_mrc is None:
        if thresh is None:
            thresh = []
            for i in range(K):
                vol = mrc.parse_mrc(f'{outdir}/kmeans{K}/vol_{i:03d}.mrc')[0]
                thresh.append(np.percentile(vol, 99.99)/2)
            thresh = np.mean(thresh)
        log(f'Threshold: {thresh}')
        log(f'Dilating mask by: {dilate}')
    
        def binary_mask(vol):
            x = (vol>=thresh).astype(bool)
            x = binary_dilation(x, iterations=dilate)
            return x
        
        # combine all masks by taking their union
        vol = mrc.parse_mrc(f'{outdir}/kmeans{K}/vol_000.mrc')[0]
        mask = ~binary_mask(vol)
        for i in range(1,K):
            vol = mrc.parse_mrc(f'{outdir}/kmeans{K}/vol_{i:03d}.mrc')[0]
            mask *= ~binary_mask(vol)
        mask = ~mask
    else:
        # Load provided mrc and convert to a boolean mask
        mask, _ = mrc.parse_mrc(in_mrc)
        mask = mask.astype(bool)

    # save mask
    out_mrc = f'{outdir}/mask.mrc'
    log(f'Saving {out_mrc}')
    mrc.write(out_mrc, mask.astype(np.float32))

    # view slices
    out_png = f'{outdir}/mask_slices.png'
    D = mask.shape[0]
    fig,ax = plt.subplots(1,3, figsize=(10,8))
    ax[0].imshow(mask[D//2,:,:])
    ax[1].imshow(mask[:,D//2,:])
    ax[2].imshow(mask[:,:,D//2])
    plt.savefig(out_png)

def choose_cmap(M):
    if M <= 10:
        cmap = 'tab10'
    elif M <= 20:
        cmap = 'tab20'
    else:
        cmap = ListedColormap(sns.color_palette('husl').as_hex())
    return cmap

def get_colors_for_cmap(cmap, M):
    if M <= 20: 
        colors = plt.cm.get_cmap(cmap)(np.arange(M)/(np.ceil(M/10)*10))
    else:
        colors = plt.cm.get_cmap(cmap)(np.linspace(0,1,M))
    return colors

def analyze_volumes(outdir, K, dim, M, linkage, vol_ind=None, plot_dim=5, particle_ind_orig=None):
    cmap = choose_cmap(M)

    # load mean volume, compute it if it does not exist
    if not os.path.exists(f'{outdir}/kmeans{K}/vol_mean.mrc'):
        volm = np.array([mrc.parse_mrc(f'{outdir}/kmeans{K}/vol_{i:03d}.mrc')[0] for i in range(K)]).mean(axis=0)
        mrc.write(f'{outdir}/kmeans{K}/vol_mean.mrc', volm)
    else:
        volm = mrc.parse_mrc(f'{outdir}/kmeans{K}/vol_mean.mrc')[0]

    # load mask
    mask = mrc.parse_mrc(f'{outdir}/mask.mrc')[0].astype(bool)
    log(f'{mask.sum()} voxels in mask')

    # load volumes
    vols = np.array([mrc.parse_mrc(f'{outdir}/kmeans{K}/vol_{i:03d}.mrc')[0][mask] for i in range(K)])
    vols[vols<0] = 0

    # load umap
    umap = utils.load_pkl(f'{outdir}/umap.pkl')
    ind = np.loadtxt(f'{outdir}/kmeans{K}/centers_ind.txt').astype(int)

    if vol_ind is not None:
        log(f'Filtering to {len(vol_ind)} volumes')
        vols = vols[vol_ind]
        ind = ind[vol_ind]

    # compute PCA
    pca = PCA(dim)
    pca.fit(vols)
    pc = pca.transform(vols)
    utils.save_pkl(pc, f'{outdir}/vol_pca_{K}.pkl')
    utils.save_pkl(pca, f'{outdir}/vol_pca_obj.pkl')
    log('Explained variance ratio:')
    log(pca.explained_variance_ratio_)

    # save rxn coordinates
    for i in range(plot_dim):
        subdir = f'{outdir}/vol_pcs/pc{i+1}'
        if not os.path.exists(subdir):
            os.makedirs(subdir)
        min_, max_ = pc[:,i].min(), pc[:,i].max()
        log((min_, max_))
        for j, val in enumerate(np.linspace(min_,max_,10,endpoint=True)):
            v = volm.copy()
            v[mask] += pca.components_[i]*val
            mrc.write(f'{subdir}/{j}.mrc', v)

    # which plots to show???
    def plot(i,j):
        plt.figure()
        plt.scatter(pc[:,i],pc[:,j])
        plt.xlabel(f'Volume PC{i+1} (EV: {pca.explained_variance_ratio_[i]:03f})')
        plt.ylabel(f'Volume PC{j+1} (EV: {pca.explained_variance_ratio_[j]:03f})')
        plt.savefig(f'{outdir}/vol_pca_{K}_{i+1}_{j+1}.png')
    for i in range(plot_dim-1):
        plot(i,i+1)

    # clustering
    subdir = f'{outdir}/clustering_L2_{linkage}_{M}'
    if not os.path.exists(subdir): 
        os.makedirs(subdir)
    cluster = AgglomerativeClustering(n_clusters=M, affinity='euclidean', linkage=linkage)
    labels = cluster.fit_predict(vols)
    utils.save_pkl(labels, f'{subdir}/state_labels.pkl')

    kmeans_labels = utils.load_pkl(f'{outdir}/kmeans{K}/labels.pkl')
    kmeans_counts = Counter(kmeans_labels)
    for i in range(M):
        vol_i = np.where(labels == i)[0]
        log(f'State {i}: {len(vol_i)} volumes')
        if vol_ind is not None:
            vol_i = np.arange(K)[vol_ind][vol_i]
        vol_i_all = np.array([mrc.parse_mrc(f'{outdir}/kmeans{K}/vol_{i:03d}.mrc')[0] for i in vol_i])
        nparticles = np.array([kmeans_counts[i] for i in vol_i]) 
        vol_i_mean = np.average(vol_i_all, axis=0, weights=nparticles)
        vol_i_std = np.average((vol_i_all-vol_i_mean)**2, axis=0, weights=nparticles)**.5
        mrc.write(f'{subdir}/state_{i}_mean.mrc', vol_i_mean.astype(np.float32))
        mrc.write(f'{subdir}/state_{i}_std.mrc', vol_i_std.astype(np.float32))
        if not os.path.exists(f'{subdir}/state_{i}'):
            os.makedirs(f'{subdir}/state_{i}')
        for v in vol_i:
            os.symlink(f'{outdir}/kmeans{K}/vol_{v:03d}.mrc', f'{subdir}/state_{i}/vol_{v:03d}.mrc')
        particle_ind = analysis.get_ind_for_cluster(kmeans_labels, vol_i)
        log(f'State {i}: {len(particle_ind)} particles')
        if particle_ind_orig is not None:
            utils.save_pkl(particle_ind_orig[particle_ind], f'{subdir}/state_{i}_particle_ind.pkl')
        else:
            utils.save_pkl(particle_ind, f'{subdir}/state_{i}_particle_ind.pkl')  


    # plot clustering results
    def hack_barplot(counts_):
        if M <= 20: # HACK TO GET COLORS
            with sns.color_palette(cmap):
                g = sns.barplot(np.arange(M), counts_)
        else: # default is husl
            g = sns.barplot(np.arange(M), counts_)
        return g
        
    plt.figure()
    counts = Counter(labels)
    g = hack_barplot([counts[i] for i in range(M)])
    for i in range(M):
        g.text(i-.1, counts[i]+2, counts[i])
    plt.xlabel('State')
    plt.ylabel('Count')
    plt.savefig(f'{subdir}/state_volume_counts.png')

    plt.figure()
    particle_counts = [np.sum([kmeans_counts[ii] for ii in np.where(labels == i)[0]]) for i in range(M)]
    g = hack_barplot(particle_counts)
    for i in range(M):
        g.text(i-.1, particle_counts[i]+2, particle_counts[i])
    plt.xlabel('State')
    plt.ylabel('Count')
    plt.savefig(f'{subdir}/state_particle_counts.png')


    def plot_w_labels(i,j):
        plt.figure()
        plt.scatter(pc[:,i],pc[:,j],c=labels, cmap=cmap)
        plt.xlabel(f'Volume PC{i+1} (EV: {pca.explained_variance_ratio_[i]:03f})')
        plt.ylabel(f'Volume PC{j+1} (EV: {pca.explained_variance_ratio_[j]:03f})')
        plt.savefig(f'{subdir}/vol_pca_{K}_{i+1}_{j+1}.png')
    for i in range(plot_dim-1):
        plot_w_labels(i,i+1)

    def plot_w_labels_annotated(i,j):
        fig, ax = plt.subplots(figsize=(16,16))
        plt.scatter(pc[:,i],pc[:,j],c=labels, cmap=cmap)
        annots = np.arange(K)
        if vol_ind is not None:
            annots = annots[vol_ind]
        for ii, k in enumerate(annots):
            ax.annotate(str(k), pc[ii,[i,j]]+np.array([.1,.1]))
        plt.xlabel(f'Volume PC{i+1} (EV: {pca.explained_variance_ratio_[i]:03f})')
        plt.ylabel(f'Volume PC{j+1} (EV: {pca.explained_variance_ratio_[j]:03f})')
        plt.savefig(f'{subdir}/vol_pca_{K}_annotated_{i+1}_{j+1}.png')
    for i in range(plot_dim-1):
        plot_w_labels_annotated(i,i+1)

    # plot clusters on UMAP
    umap_i = umap[ind]
    fig, ax = plt.subplots(figsize=(8,8))
    plt.scatter(umap[:,0], umap[:,1], alpha=.1, s=1, rasterized=True, color='lightgrey')
    colors = get_colors_for_cmap(cmap, M)
    for i in range(M):
        c = umap_i[np.where(labels == i)]
        plt.scatter(c[:,0], c[:,1], label=i, color=colors[i])
    plt.legend()
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.savefig(f'{subdir}/umap.png')

    fig, ax = plt.subplots(figsize=(16,16))
    plt.scatter(umap[:,0], umap[:,1], alpha=.1, s=1, rasterized=True, color='lightgrey')
    plt.scatter(umap_i[:,0], umap_i[:,1], c=labels, cmap=cmap)
    annots = np.arange(K)
    if vol_ind is not None:
        annots = annots[vol_ind]
    for i,k in enumerate(annots):
        ax.annotate(str(k),umap_i[i]+np.array([.1,.1]))
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.savefig(f'{subdir}/umap_annotated.png')

def main(args):
    t1 = dt.now()
    log(args)
    E = args.epoch
    workdir = args.workdir
    zfile = f'{workdir}/z.{E}.pkl'
    weights = f'{workdir}/weights.{E}.pkl'
    config = f'{workdir}/config.pkl'
    outdir = f'{workdir}/landscape.{E}'
    
    if args.outdir:
        outdir = args.outdir
    log(f'Saving results to {outdir}')
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    z = utils.load_pkl(zfile)
    zdim = z.shape[1]
    K = args.sketch_size

    vol_args = dict(Apix=args.Apix, downsample=args.downsample, flip=args.flip, cuda=args.device)
    vg = VolumeGenerator(weights, config, vol_args, skip_vol=args.skip_vol)

    if args.vol_ind is not None:
        args.vol_ind = utils.load_pkl(args.vol_ind)

    if not args.skip_vol:
        generate_volumes(z, outdir, vg, K)
    else:
        log('Skipping volume generation')

    if args.skip_umap:
        assert os.path.exists(f'{outdir}/umap.pkl')
        log('Skipping UMAP')
    else:
        log(f'Copying UMAP from {workdir}/analyze.{E}/umap.pkl')
        if os.path.exists(f'{workdir}/analyze.{E}/umap.pkl'):
            from shutil import copyfile
            copyfile(f'{workdir}/analyze.{E}/umap.pkl', f'{outdir}/umap.pkl')
        else:
            raise NotImplementedError
        
    if args.mask:
        log(f'Using custom mask {args.mask}')
    make_mask(outdir, K, args.dilate, args.thresh, args.mask)

    log('Analyzing volumes...')
    # get particle indices if the dataset was originally filtered
    c = utils.load_pkl(config)
    particle_ind = utils.load_pkl(c['dataset_args']['ind']) if c['dataset_args']['ind'] is not None else None
    analyze_volumes(outdir, K, args.pc_dim, args.M, args.linkage, vol_ind=args.vol_ind, plot_dim=args.plot_dim, particle_ind_orig=particle_ind)
    td = dt.now()-t1
    log(f'Finished in {td}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
