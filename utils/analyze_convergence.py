'''
Visualize convergence and training dynamics 

(BETA -- contributed by Barrett Powell bmp@mit.edu)
'''

import argparse
import sys, os
import numpy as np
from matplotlib import pyplot as plt
import random
from datetime import datetime as dt
import umap
from cryodrgn import analysis, starfile, utils, mrc, fft
from scipy.spatial import distance_matrix
from scipy.ndimage.filters import maximum_filter, gaussian_filter
from scipy.ndimage.morphology import distance_transform_edt, binary_dilation
from scipy import stats, ndimage
from string import ascii_uppercase
import multiprocessing
import itertools
try:
    from cuml.manifold.umap import UMAP as cuUMAP
except ImportError:
    pass
flog = utils.flog

def add_args(parser):
    parser.add_argument('workdir', type=os.path.abspath, help='Directory with cryoDRGN results')
    parser.add_argument('epoch', type=int, help='Latest epoch number N to analyze convergence (0-based indexing, corresponding to z.N.pkl, weights.N.pkl')
    parser.add_argument('-o', '--outdir', type=os.path.abspath, help='Output directory for convergence analysis results (default: [workdir]/convergence.[epoch])')
    parser.add_argument('--epoch-interval', type=int, default=5, help='Interval of epochs between calculating most convergence heuristics')

    group = parser.add_argument_group('UMAP calculation arguments')
    group.add_argument('--force-umap-cpu', action='store_true', help='Override default UMAP GPU-bound implementation via cuML to use umap-learn library instead')
    group.add_argument('--subset', default = 50000, help='Max number of particles to be used for UMAP calculations. \'None\' means use all ptcls')
    group.add_argument('--random-seed', default = None, help='Manually specify the seed used for selection of subset particles')
    group.add_argument('--random-state', type = int, default=42, help='Random state seed used by UMAP for reproducibility at slight cost of performance (default 42, None means slightly faster but non-reproducible)')
    group.add_argument('--n-epochs-umap', type=int, default=25000, help='Number of epochs to train the UMAP embedding via cuML for a given z.pkl, as described in the cuml.UMAP documentation')
    group.add_argument('--skip-umap', action = 'store_true', help='Skip UMAP embedding. Requires that UMAP be precomputed for downstream calcs. Useful for tweaking volume generation settings.')

    group = parser.add_argument_group('Sketching UMAP via local maxima arguments')
    group.add_argument('--n-bins', type = int, default = 30, help='the number of bins along UMAP1 and UMAP2')
    group.add_argument('--smooth', type = bool, default = True, help='smooth the 2D histogram before identifying local maxima')
    group.add_argument('--smooth-width', type = float, default = 1.0, help='width of gaussian kernel for smoothing 2D histogram expressed as multiple of one bin\'s width')
    group.add_argument('--pruned-maxima', type = int, default = 12, help='prune poorly-separated maxima until this many maxima remain')
    group.add_argument('--radius', type = float, default = 5.0, help='distance at which two maxima are considered poorly-separated and are candidates for pruning (euclidean distance in bin-space)')
    group.add_argument('--final-maxima', type = int, default = 10, help='select this many local maxima, sorted by highest bin count after pruning, for which to generate volumes')

    group = parser.add_argument_group('Volume generation arguments')
    group.add_argument('--Apix', type = float, default = 1.0, help='A/pix of output volume')
    group.add_argument('--flip', action='store_true', help='Flip handedness of output volume')
    group.add_argument('--invert', action='store_true', help='Invert contrast of output volume')
    group.add_argument('-d','--downsample', type=int, help='Downsample volumes to this box size (pixels). Recommended for boxes > 250-300px')
    group.add_argument('--cuda', type = int, default = None, help='Specify cuda device for volume generation')
    group.add_argument('--skip-volgen', action = 'store_true', help='Skip volume generation. Requires that volumes already exist for downstream CC + FSC calcs')

    group = parser.add_argument_group('Mask generation arguments')
    group.add_argument('--max-threads', type=int, default=8, help='Max number of threads used to parallelize mask generation')
    group.add_argument('--thresh', type = float, default = None, help='Float, isosurface at which to threshold mask; default None uses 50th percentile')
    group.add_argument('--dilate', type = int, default = 3, help='Number of voxels to dilate thresholded isosurface outwards from mask boundary')
    group.add_argument('--dist', type = int, default = 10, help='Number of voxels over which to apply a soft cosine falling edge from dilated mask boundary')

    return parser


def plot_loss(logfile, outdir, E, LOG):
    '''
    Plots the total loss (reconstruction + regularization) per epoch

    Inputs
        logfile: the run.log auto-generated by cryodrgn train_vae
        outdir: path to base directory to save outputs

    Outputs
        png of total loss vs epochs
    '''

    loss = analysis.parse_loss(logfile)
    plt.plot(loss[:E])
    plt.xlabel('epoch')
    plt.ylabel('total loss')
    plt.savefig(outdir + '/plots/00_total_loss.png',
                dpi=300,
                format='png',
                transparent=True,
                bbox_inches='tight')
    flog(f'Saved total loss plot to {outdir}/plots/00_total_loss.png', LOG)


def encoder_latent_umaps(workdir, outdir, epochs, n_particles_total, subset, random_seed, use_umap_gpu, random_state, n_epochs_umap, LOG):
    '''
    Calculates UMAP embeddings of subset of particles' selected epochs' latent encodings

    Inputs
        workdir: path to directory containing cryodrgn training results
        outdir: path to base directory to save outputs
        epochs: array of epochs for which to calculate UMAPs
        n_particles_total: int of total number of particles trained
        subset: int, size of subset on which to calculate umap, None means all
        random_seed: int, seed for random selection of subset particles
        use_umap_gpu: bool, whether to use the cuML library to GPU accelerate UMAP calculations (if available in env)
        random_state: int, random state seed used by UMAP for reproducibility at slight cost of performance (None means faster but non-reproducible)

    Outputs
        pkl of each UMAP embedding stored in outdir/umaps/umap.epoch.pkl
        png of all UMAPs

    # apparently running multiple UMAP embeddings (i.e. for each epoch's z.pkl) in parallel on CPU requires difficult backend setup
    # see https://github.com/lmcinnes/umap/issues/707
    # therefore not implemented currently
    '''

    if subset == 'None':
        n_particles_subset = n_particles_total
        flog('Using full particle stack for UMAP', LOG)
    else:
        if random_seed is None:
            random_seed = random.randint(0, 100000)
            random.seed(random_seed)
        else:
            random.seed(random_seed)
        n_particles_subset = min(n_particles_total, int(subset))
        flog(f'Randomly selecting {n_particles_subset} particle subset on which to run UMAP (with random seed {random_seed})', LOG)
    ind_subset = sorted(random.sample(range(0, n_particles_total), k=n_particles_subset))
    utils.save_pkl(ind_subset, outdir + '/ind_subset.pkl')

    for epoch in epochs:
        flog(f'Now calculating UMAP for epoch {epoch} with random_state {random_state}', LOG)
        z = utils.load_pkl(workdir + f'/z.{epoch}.pkl')[ind_subset, :]
        if use_umap_gpu: #using cuML library GPU-accelerated UMAP
            reducer = cuUMAP(random_state=random_state, n_epochs=n_epochs_umap)
            umap_embedding = reducer.fit_transform(z)
        else: #using umap-learn library CPU-bound UMAP
            reducer = umap.UMAP(random_state=random_state)
            umap_embedding = reducer.fit_transform(z)
        utils.save_pkl(umap_embedding, f'{outdir}/umaps/umap.{epoch}.pkl')


    n_cols = int(np.ceil(len(epochs) ** 0.5))
    n_rows = int(np.ceil(len(epochs) / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(2 * n_cols, 2 * n_rows), sharex='all', sharey='all')
    fig.tight_layout()

    for i, ax in enumerate(axes.flat):
        try:
            umap_embedding = utils.load_pkl(f'{outdir}/umaps/umap.{epochs[i]}.pkl')
            toplot = ax.hexbin(umap_embedding[:, 0], umap_embedding[:, 1], bins='log', mincnt=1)
            ax.set_title(f'epoch {epochs[i]}')
        except IndexError:
            pass
        except FileNotFoundError:
            flog(f'Could not find file {outdir}/umaps/umap.{epoch}.pkl', LOG)
            pass

    if len(axes.shape) == 1:
        axes[0].set_ylabel('UMAP2')
        for a in axes[:]: a.set_xlabel('UMAP1')
    else:
        assert len(axes.shape) == 2 #there are more than one row and column of axes
        for a in axes[:, 0]: a.set_ylabel('UMAP2')
        for a in axes[-1, :]: a.set_xlabel('UMAP1')
    fig.subplots_adjust(right=0.96)
    cbar_ax = fig.add_axes([0.98, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(toplot, cax=cbar_ax)
    cbar.ax.set_ylabel('particle density', rotation=90)

    plt.subplots_adjust(wspace=0.1)
    plt.subplots_adjust(hspace=0.3)
    plt.savefig(f'{outdir}/plots/01_encoder_umaps.png', dpi=300, format='png', transparent=True, bbox_inches='tight')
    flog(f'Saved UMAP distribution plot to {outdir}/plots/01_encoder_umaps.png', LOG)


def encoder_latent_shifts(workdir, outdir, epochs, E, LOG):
    '''
    Calculates and plots various metrics characterizing the per-particle latent vectors between successive epochs.

    Inputs
        workdir: path to directory containing cryodrgn training results
        outdir: path to base directory to save outputs
        E: int of epoch from which to evaluate convergence (0-indexed)
            Note that because three epochs are needed to define the two inter-epoch vectors analyzed "per epoch",
            the output contains metrics for E-2 epochs. Accordingly, plot x-axes are labeled from epoch 2 - epoch E.

    Outputs
        pkl of all statistics of shape(E, n_metrics)
        png of each statistic plotted over training
    '''
    metrics = ['dot product', 'magnitude', 'cosine distance']

    vector_metrics = np.zeros((E-1, len(metrics)))
    for i in np.arange(E-1):
        flog(f'Calculating vector metrics for epochs {i}-{i+1} and {i+1}-{i+2}', LOG)
        if i == 0:
            z1 = utils.load_pkl(f'{workdir}/z.{i}.pkl')
            z2 = utils.load_pkl(f'{workdir}/z.{i+1}.pkl')
            z3 = utils.load_pkl(f'{workdir}/z.{i+2}.pkl')
        else:
            z1 = z2.copy()
            z2 = z3.copy()
            z3 = utils.load_pkl(workdir + f'/z.{i+2}.pkl')

        diff21 = z2 - z1
        diff32 = z3 - z2

        vector_metrics[i, 0] = np.median(np.einsum('ij,ij->i', diff21, diff32), axis=0)  # median vector dot product
        vector_metrics[i, 1] = np.median(np.linalg.norm(diff32, axis=1), axis=0)  # median vector magnitude
        uv = np.sum(diff32 * diff21, axis=1)
        uu = np.sum(diff32 * diff32, axis=1)
        vv = np.sum(diff21 * diff21, axis=1)
        vector_metrics[i, 2] = np.median(1 - uv / (np.sqrt(uu) * np.sqrt(vv))) #median vector cosine distance

    utils.save_pkl(vector_metrics, f'{outdir}/vector_metrics.pkl')

    fig, axes = plt.subplots(1, len(metrics), figsize=(10,3))
    fig.tight_layout()
    for i,ax in enumerate(axes.flat):
        ax.plot(np.arange(2,E+1), vector_metrics[:,i])
        ax.set_xlabel('epoch')
        ax.set_ylabel(metrics[i])
    plt.savefig(f'{outdir}/plots/02_encoder_latent_vector_shifts.png', dpi=300, format='png', transparent=True, bbox_inches='tight')

    flog(f'Saved latent vector shifts plots to {outdir}/plots/02_encoder_latent_vector_shifts.png', LOG)


def sketch_via_umap_local_maxima(outdir, E, LOG, n_bins=30, smooth=True, smooth_width=1, pruned_maxima=12, radius=5, final_maxima=10):
    '''
    Sketch the UMAP embedding of epoch E latent space via local maxima finding

    Inputs:
        E: epoch for which the (subset, see Convergence 1) umap distribution will be sketched for local maxima
        n_bins: the number of bins along UMAP1 and UMAP2
        smooth: whether to smooth the 2D histogram (aids local maxima finding for particulaly continuous distributions)
        smooth_width: scalar multiple of one-bin-width defining sigma for gaussian kernel smoothing
        pruned_maxima: max number of local maxima above which pruning will be attempted
        radius: radius in bin-space (Euclidean distance) below which points are considered poorly-separated and are candidates for pruning
        final_maxima: the count of local maxima with highest associated bin count that will be returned as final to the user

    Outputs
        binned_ptcls_mask: binary mask of shape ((n_particles_total, n_local_maxima)) labeling all particles in the bin and neighboring 8 bins of a local maxima
        labels: a unique letter assigned to each local maxima
    '''
    def make_edges(umap, n_bins):
        '''
        Helper function to create two 1-D arrays defining @nbins bin edges along axes x and y
        '''
        xedges = np.linspace(umap.min(axis=0)[0], umap.max(axis=0)[0], n_bins + 1)
        yedges = np.linspace(umap.min(axis=0)[1], umap.max(axis=0)[1], n_bins + 1)
        return xedges, yedges

    def local_maxima_2D(data):
        '''
        Helper function to find the coordinates and values of local maxima of a 2d hist
        Evaluates local maxima using a footprint equal to 3x3 set of bins
        '''
        size = 3
        footprint = np.ones((size, size))
        footprint[1, 1] = 0

        filtered = maximum_filter(data, footprint=footprint, mode='mirror')
        mask_local_maxima = data > filtered
        coords = np.asarray(np.where(mask_local_maxima)).T
        values = data[mask_local_maxima]

        return coords, values

    def gen_peaks_img(coords, values, edges):
        '''
        Helper function to scatter the values of the local maxima onto a hist with bins defined by the full umap
        '''
        filtered = np.zeros((edges[0].shape[0], edges[1].shape[0]))
        for peak in range(coords.shape[0]):
            filtered[tuple(coords[peak])] = values[peak]
        return filtered

    def prune_local_maxima(coords, values, n_maxima, radius):
        '''
        Helper function to prune "similar" local maxima and preserve UMAP diversity if more local maxima than desired are found
        Construct distance matrix of all coords to all coords in bin-space
        Find all maxima pairs closer than @radius
        While more than @n_maxima local maxima:
            if there are pairs closer than @radius:
                find single smallest distance d between two points
                compare points connected by d, remove lower value point from coords, values, and distance matrix
        Returns
        * coords
        * values
        '''
        dist_matrix = distance_matrix(coords, coords)
        dist_matrix[dist_matrix > radius] = 0  # ignore points separated by > @radius in bin-space

        while len(values) > n_maxima:
            if not np.count_nonzero(dist_matrix) == 0:  # some peaks are too close and need pruning
                indices_to_compare = np.where(dist_matrix == np.min(dist_matrix[np.nonzero(dist_matrix)]))[0]
                if values[indices_to_compare[0]] > values[indices_to_compare[1]]:
                    dist_matrix = np.delete(dist_matrix, indices_to_compare[1], axis=0)
                    dist_matrix = np.delete(dist_matrix, indices_to_compare[1], axis=1)
                    values = np.delete(values, indices_to_compare[1])
                    coords = np.delete(coords, indices_to_compare[1], axis=0)
                else:
                    dist_matrix = np.delete(dist_matrix, indices_to_compare[0], axis=0)
                    dist_matrix = np.delete(dist_matrix, indices_to_compare[0], axis=1)
                    values = np.delete(values, indices_to_compare[0])
                    coords = np.delete(coords, indices_to_compare[0], axis=0)
            else:  # local maxima are already well separated
                return coords, values
        return coords, values

    def coords_to_umap(umap, binned_ptcls_mask, values):
        '''
        Helper function to convert local maxima coords in bin-space to umap-space
        Calculates each local maximum to be the median UMAP1 and UMAP2 value across all particles in each 3x3 set of bins defining a given local maximum
        '''
        umap_median_peaks = np.zeros((len(values), 2))
        for i in range(len(values)):
            umap_median_peaks[i, :] = np.median(umap[binned_ptcls_mask[:, i], :], axis=0)
        return umap_median_peaks

    flog('Using UMAP local maxima sketching', LOG)
    umap = utils.load_pkl(outdir + f'/umaps/umap.{E}.pkl')
    n_particles_sketch = umap.shape[0]

    # create 2d histogram of umap distribution
    edges = make_edges(umap, n_bins=n_bins)
    hist, xedges, yedges, bincount = stats.binned_statistic_2d(umap[:, 0], umap[:, 1], None, 'count', bins=edges, expand_binnumbers=True)
    to_plot = ['umap', 'hist']

    # optionally smooth the histogram to reduce the number of peaks with sigma=width of two bins
    if smooth:
        hist_smooth = gaussian_filter(hist, smooth_width * np.abs(xedges[1] - xedges[0]))
        coords, values = local_maxima_2D(hist_smooth)
        to_plot[-1] = 'hist_smooth'
    else:
        coords, values = local_maxima_2D(hist)
    flog(f'Found {len(values)} local maxima', LOG)

    # prune local maxima that are densely packed and low in value
    coords, values = prune_local_maxima(coords, values, pruned_maxima, radius)
    flog(f'Pruned to {len(values)} local maxima', LOG)

    # find subset of n_peaks highest local maxima
    indices = (-values).argsort()[:final_maxima]
    coords, values = coords[indices], values[indices]
    peaks_img_top = gen_peaks_img(coords, values, edges)
    to_plot.append('peaks_img_top')
    to_plot.append('sketched_umap')
    flog(f'Filtered to top {len(values)} local maxima', LOG)

    # write list of lists containing indices of all particles within maxima bins + all 8 neighboring bins (assumes footprint = (3,3))
    binned_ptcls_mask = np.zeros((n_particles_sketch, len(values)), dtype=bool)
    for i in range(len(values)):
        binned_ptcls_mask[:, i] = (bincount[0, :] >= coords[i, 0] + 0) & \
                                  (bincount[0, :] <= coords[i, 0] + 2) & \
                                  (bincount[1, :] >= coords[i, 1] + 0) & \
                                  (bincount[1, :] <= coords[i, 1] + 2)

    # find median umap coords of each maxima bin for plotting
    coords = coords_to_umap(umap, binned_ptcls_mask, values)

    # plot the original histogram, all peaks, and highest n_peaks
    fig, axes = plt.subplots(1, len(to_plot), figsize=(len(to_plot) * 3.6, 3))
    fig.tight_layout()
    labels = ascii_uppercase[:len(values)]
    for i, ax in enumerate(axes.flat):
        if to_plot[i] == 'umap':
            ax.hexbin(umap[:, 0], umap[:, 1], mincnt=1)
            ax.vlines(x=xedges, ymin=umap.min(axis=0)[1], ymax=umap.max(axis=0)[1], colors='red', linewidth=0.35)
            ax.hlines(y=yedges, xmin=umap.min(axis=0)[0], xmax=umap.max(axis=0)[0], colors='red', linewidth=0.35)
            ax.set_title(f'epoch {E} UMAP')
            ax.set_xlabel('UMAP1')
            ax.set_ylabel('UMAP2')
        elif to_plot[i] == 'hist':
            ax.imshow(np.rot90(hist))
            ax.set_title('UMAP histogram')
        elif to_plot[i] == 'hist_smooth':
            ax.imshow(np.rot90(hist_smooth))
            ax.set_title('UMAP smoothed histogram')
        elif to_plot[i] == 'peaks_img_top':
            ax.imshow(np.rot90(peaks_img_top))
            ax.set_title(f'final {len(labels)} local maxima')
        elif to_plot[i] == 'sketched_umap':
            ax.hexbin(umap[:, 0], umap[:, 1], mincnt=1)
            ax.scatter(*coords.T, c='r')
            ax.set_title(f'sketched epoch {E} UMAP')
            ax.set_xlabel('UMAP1')
            ax.set_ylabel('UMAP2')
            for k in range(len(values)):
                ax.text(x=coords[k, 0] + 0.3,
                        y=coords[k, 1] + 0.3,
                        s=labels[k],
                        fontdict=dict(color='r', size=10))
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.savefig(f'{outdir}/plots/03_decoder_UMAP-sketching.png', dpi=300, format='png', transparent=True, bbox_inches='tight')
    flog(f'Saved latent sketching plot to {outdir}/plots/03_decoder_UMAP-sketching.png', LOG)

    return binned_ptcls_mask, labels


def follow_candidate_particles(workdir, outdir, epochs, n_dim, binned_ptcls_mask, labels, LOG):
    '''
    Monitor how the labeled set of particles migrates within latent space at selected epochs over training

    Inputs:
        workdir: path to directory containing cryodrgn training results
        outdir: path to base directory to save outputs
        epochs: array of epochs for which to calculate UMAPs
        n_dim: latent dimensionality
        binned_ptcls_mask: (n_particles, len(labels)) binary mask of which particles belong to which class
        labels: unique identifier for each class of representative latent encodings

    Outputs
        plot.png tracking representative latent encodings through epochs
        latent.txt of representative latent encodings for each epoch
    '''

    # track sketched points from epoch E through selected previous epochs and plot overtop UMAP embedding
    n_cols = int(np.ceil(len(epochs) ** 0.5))
    n_rows = int(np.ceil(len(epochs) / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(2 * n_cols, 2 * n_rows), sharex='all', sharey='all')
    fig.tight_layout()

    ind_subset = utils.load_pkl(f'{outdir}/ind_subset.pkl')
    for i, ax in enumerate(axes.flat):
        try:
            umap = utils.load_pkl(f'{outdir}/umaps/umap.{epochs[i]}.pkl')
            z = utils.load_pkl(f'{workdir}/z.{epochs[i]}.pkl')[ind_subset,:]
            z_maxima_median = np.zeros((len(labels), n_dim))

            for k in range(len(labels)):
                z_maxima_median[k, :] = np.median(z[binned_ptcls_mask[:, k]], axis=0) # find median latent value of each maximum in a given epoch

            z_maxima_median_ondata, z_maxima_median_ondata_ind = analysis.get_nearest_point(z, z_maxima_median)  # find on-data latent encoding of each median latent value
            umap_maxima_median_ondata = umap[z_maxima_median_ondata_ind] # find on-data UMAP embedding of each median latent encoding

            # Write out the on-data median latent values of each labeled set of particles for each epoch in epochs
            with open(f'{outdir}/repr_particles/latent_representative.{epochs[i]}.txt', 'w') as f:
                np.savetxt(f, z_maxima_median_ondata, delimiter=' ', newline='\n', header='', footer='', comments='# ')
            flog(f'Saved representative latent encodings for epoch {epochs[i]} to {outdir}/repr_particles/latent_representative.{epochs[i]}.txt', LOG)

            for k in range(len(labels)):
                ax.text(x=umap_maxima_median_ondata[k, 0] + 0.3,
                        y=umap_maxima_median_ondata[k, 1] + 0.3,
                        s=labels[k],
                        fontdict=dict(color='r', size=10))
            toplot = ax.hexbin(*umap.T, bins='log', mincnt=1)
            ax.scatter(umap_maxima_median_ondata[:, 0], umap_maxima_median_ondata[:, 1], s=10, linewidth=0, c='r',
                       alpha=1)
            ax.set_title(f'epoch {epochs[i]}')
        except IndexError:
            pass

    if len(axes.shape) == 1:
        axes[0].set_ylabel('UMAP2')
        for a in axes[:]: a.set_xlabel('UMAP1')
    else:
        assert len(axes.shape) == 2 #there are more than one row and column of axes
        for a in axes[:, 0]: a.set_ylabel('UMAP2')
        for a in axes[-1, :]: a.set_xlabel('UMAP1')
    fig.subplots_adjust(right=0.96)
    cbar_ax = fig.add_axes([0.98, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(toplot, cax=cbar_ax)
    cbar.ax.set_ylabel('Particle Density', rotation=90)

    plt.subplots_adjust(wspace=0.1)
    plt.subplots_adjust(hspace=0.25)

    plt.savefig(f'{outdir}/plots/04_decoder_maxima-sketch-consistency.png', dpi=300, format='png', transparent=True, bbox_inches='tight')
    flog(f'Saved plot tracking representative latent encodings through epochs {epochs} to {outdir}/plots/04_decoder_maxima-sketch-consistency.png', LOG)


def generate_volumes(workdir, outdir, epochs, Apix, flip, invert, downsample, cuda, LOG):
    '''
    Helper function to call cryodrgn.analysis.gen_volumes on all representative latent values in selected epochs
    '''
    for epoch in epochs:
        weights = f'{workdir}/weights.{epoch}.pkl'
        config = f'{workdir}/config.pkl'
        zfile = f'{outdir}/repr_particles/latent_representative.{epoch}.txt'
        volsdir = f'{outdir}/vols.{epoch}'

        analysis.gen_volumes(weights, config, zfile, volsdir, Apix=Apix, flip=flip, invert=invert, downsample=downsample, cuda=cuda)


def mask_volume(volpath, outpath, Apix, thresh=None, dilate=3, dist=10):
    '''
    Helper function to generate a loose mask around the input density
    Density is thresholded to 50% maximum intensity, dilated outwards, and a soft cosine edge is applied

    Inputs
        volpath: an absolute path to the volume to be used for masking
        outpath: an absolute path to write out the mask mrc
        thresh: what intensity threshold between [0, 100] to apply
        dilate: how far to dilate the thresholded density outwards
        dist: how far the cosine edge extends from the density

    Outputs
       volume.masked.mrc written to outdir
    '''
    vol = mrc.parse_mrc(volpath)[0]
    thresh = np.percentile(vol, 99.99) / 2 if thresh is None else thresh
    x = (vol >= thresh).astype(bool)
    x = binary_dilation(x, iterations=dilate)
    y = distance_transform_edt(~x.astype(bool))
    y[y > dist] = dist
    z = np.cos(np.pi * y / dist / 2)

    # check that mask is in range [0,1]
    assert np.all(z >= 0)
    assert np.all(z <= 1)

    # used to write out mask separately from masked volume, now apply and save the masked vol to minimize future I/O
    # mrc.write(outpath, z.astype(np.float32))
    vol *= z
    mrc.write(outpath, vol.astype(np.float32), Apix=Apix)


def mask_volumes(outdir, epochs, labels, max_threads, LOG, Apix, thresh=None, dilate=3, dist=10):
    '''
    Generate a loose mask around each volume in outdir/vols.{epochs}

    Inputs:
        outdir: path to base directory to save outputs
        epochs: array of epochs for which to calculate UMAPs
        labels: unique identifier for each class of representative latent encodings
        thresh: isosurface at which to threshold density when generating mask (default: 50th percentile)
        dilate: number of voxels to dilate thresholded isosurface outwards from mask boundary
        dist: number of voxels over which to apply a soft cosine falling edge from dilated mask boundary

    Outputs:
        volume.masked.mrc for each volume
    '''

    volpaths = []
    outpaths = []
    for epoch in epochs:
        flog(f'Generating and applying masks for epoch {epoch}', LOG)
        volsdir = outdir + f'/vols.{epoch}'
        for cluster in range(len(labels)):
            volpath = f'{volsdir}/vol_{cluster:03d}.mrc'
            outpath = f'{volsdir}/vol_{cluster:03d}.masked.mrc'

            volpaths.append(volpath)
            outpaths.append(outpath)

    args = zip(volpaths, outpaths, itertools.repeat(Apix), itertools.repeat(thresh), itertools.repeat(dilate), itertools.repeat(dist))
    with multiprocessing.Pool(max_threads) as p:
        p.starmap(mask_volume, args, 4)


def calculate_CCs(outdir, epochs, labels, chimerax_colors, LOG):
    '''
    Returns the masked map-map correlation between temporally sequential volume pairs outdir/vols.{epochs}, for each class in labels

    Inputs:
        outdir: path to base directory to save outputs
        epochs: array of epochs for which to calculate UMAPs
        labels: unique identifier for each class of representative latent encodings
        chimerax_colors: approximate colors matching ChimeraX palette to facilitate comparison to volume visualization

    Outputs:
        plot.png of sequential volume pairs map-map CC for each class in labels across training epochs
    '''
    def calc_cc(vol1, vol2):
        '''
        Helper function to calculate the zero-mean correlation coefficient as defined in eq 2 in https://journals.iucr.org/d/issues/2018/09/00/kw5139/index.html
        vol1 and vol2 should be maps of the same box size, structured as numpy arrays with ndim=3, i.e. by loading with cryodrgn.mrc.parse_mrc
        '''
        zmean1 = (vol1 - np.mean(vol1))
        zmean2 = (vol2 - np.mean(vol2))
        cc = (np.sum(zmean1 ** 2) ** -0.5) * (np.sum(zmean2 ** 2) ** -0.5) * np.sum(zmean1 * zmean2)
        return cc

    cc_masked = np.zeros((len(labels), len(epochs) - 1))

    for i in range(len(epochs) - 1):
        for cluster in np.arange(len(labels)):
            vol1, _ = mrc.parse_mrc(f'{outdir}/vols.{epochs[i]}/vol_{cluster:03d}.masked.mrc')
            vol2, _ = mrc.parse_mrc(f'{outdir}/vols.{epochs[i+1]}/vol_{cluster:03d}.masked.mrc')

            cc_masked[cluster, i] = calc_cc(vol1, vol2)

    utils.save_pkl(cc_masked, f'{outdir}/cc_masked.pkl')

    fig, ax = plt.subplots(1, 1)

    ax.set_xlabel('epoch')
    ax.set_ylabel('masked CC')
    for i in range(len(labels)):
        ax.plot(epochs[1:], cc_masked[i,:], c=chimerax_colors[i] * 0.75, linewidth=2.5)
    ax.legend(labels, ncol=3, fontsize='x-small')

    plt.savefig(f'{outdir}/plots/05_decoder_CC.png', dpi=300, format='png', transparent=True, bbox_inches='tight')
    flog(f'Saved map-map correlation plot to {outdir}/plots/05_decoder_CC.png', LOG)


def calculate_FSCs(outdir, epochs, labels, img_size, chimerax_colors, LOG):
    '''
    Returns the masked FSC between temporally sequential volume pairs outdir/vols.{epochs}, for each class in labels

    Inputs:
        outdir: path to base directory to save outputs
        epochs: array of epochs for which to calculate UMAPs
        labels: unique identifier for each class of representative latent encodings
        img_size: box size of input images in pixels
        chimerax_colors: approximate colors matching ChimeraX palette to facilitate comparison to volume visualization

    Outputs:
        plot.png of sequential volume pairs map-map FSC for each class in labels across training epochs
        plot.png of sequential volume pairs map-map FSC at Nyquist for each class in labels across training epochs

    TODO: accelerate via multiprocessing (create iterable list of calc_fsc calls?)

    '''
    def calc_fsc(vol1_path, vol2_path):
        '''
        Helper function to calculate the FSC between two (assumed masked) volumes
        vol1 and vol2 should be maps of the same box size, structured as numpy arrays with ndim=3, i.e. by loading with cryodrgn.mrc.parse_mrc
        '''
        # load masked volumes in fourier space
        vol1, _ = mrc.parse_mrc(vol1_path)
        vol2, _ = mrc.parse_mrc(vol2_path)

        vol1_ft = fft.fftn_center(vol1)
        vol2_ft = fft.fftn_center(vol2)

        # define fourier grid and label into shells
        D = vol1.shape[0]
        x = np.arange(-D // 2, D // 2)
        x0, x1, x2 = np.meshgrid(x, x, x, indexing='ij')
        r = np.sqrt(x0 ** 2 + x1 ** 2 + x2 ** 2)
        r_max = D // 2  # sphere inscribed within volume box
        r_step = 1  # int(np.min(r[r>0]))
        bins = np.arange(0, r_max, r_step)
        bin_labels = np.searchsorted(bins, r, side='right')

        # calculate the FSC via labeled shells
        num = ndimage.sum(np.real(vol1_ft * np.conjugate(vol2_ft)), labels=bin_labels, index=bins + 1)
        den1 = ndimage.sum(np.abs(vol1_ft) ** 2, labels=bin_labels, index=bins + 1)
        den2 = ndimage.sum(np.abs(vol2_ft) ** 2, labels=bin_labels, index=bins + 1)
        fsc = num / np.sqrt(den1 * den2)

        x = bins / D  # x axis should be spatial frequency in 1/px
        return x, fsc

    # calculate masked FSCs for all volumes
    fsc_masked = np.zeros((len(labels), len(epochs) - 1, img_size // 2))

    for cluster in range(len(labels)):
        flog(f'Calculating all FSCs for cluster {cluster}', LOG)

        for i in range(len(epochs) - 1):
            vol1_path = f'{outdir}/vols.{epochs[i]}/vol_{cluster:03d}.masked.mrc'
            vol2_path = f'{outdir}/vols.{epochs[i+1]}/vol_{cluster:03d}.masked.mrc'

            x, fsc_masked[cluster, i, :] = calc_fsc(vol1_path, vol2_path)

    utils.save_pkl(fsc_masked, f'{outdir}/fsc_masked.pkl')
    utils.save_pkl(x, f'{outdir}/fsc_xaxis.pkl')

    # plot all fscs
    n_cols = int(np.ceil(len(labels) ** 0.5))
    n_rows = int(np.ceil(len(labels) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(2 * n_cols, 2 * n_rows), sharex='all', sharey='all')
    fig.tight_layout()
    for cluster, ax in enumerate(axes.flat):
        try:
            colors = plt.cm.viridis(np.linspace(0, 1, len(epochs - 1)))
            ax.set_ylim(0, 1.02)
            ax.set_title(f'maximum {labels[cluster]}')
            legend = []
            for i in range(len(epochs) - 1):
                ax.plot(x, fsc_masked[cluster, i, :], color=colors[i])
                legend.append(f'epoch {epochs[i+1]}')
        except IndexError:
            pass

    x_center, y_center = n_cols//2, n_rows//2
    axes[y_center, 0].set_ylabel('FSC')
    axes[-1,x_center].set_xlabel('frequency (1/px)')
    axes[-1, 0].legend(legend, loc='lower left', ncol=2, fontsize=6.5)
    plt.subplots_adjust(hspace=0.3)
    plt.subplots_adjust(wspace=0.1)
    plt.savefig(f'{outdir}/plots/06_decoder_FSC.png', dpi=300, format='png', transparent=True, bbox_inches='tight')
    flog(f'Saved map-map FSC plot to {outdir}/plots/06_decoder_FSC.png', LOG)

    # plot all FSCs at Nyquist only
    fig, ax = plt.subplots(1, 1)

    ax.set_xlabel('epoch')
    ax.set_ylabel('masked FSC at nyquist')
    for i in range(len(labels)):
        ax.plot(epochs[1:], fsc_masked[i, :, -1], c=chimerax_colors[i] * 0.75, linewidth=2.5)
    ax.legend(labels, ncol=3, fontsize='x-small')

    plt.savefig(f'{outdir}/plots/07_decoder_FSC-nyquist.png', dpi=300, format='png', transparent=True, bbox_inches='tight')
    flog(f'Saved map-map FSC (Nyquist) plot to {outdir}/plots/07_decoder_FSC-nyquist.png', LOG)


def main(args):
    t1 = dt.now()

    # Configure paths
    E = args.epoch
    sampling = args.epoch_interval
    epochs = np.arange(4, E+1, sampling)
    if epochs[-1] != E:
        epochs = np.append(epochs, E)
    workdir = args.workdir
    config = f'{workdir}/config.pkl'
    logfile = f'{workdir}/run.log'

    # assert all required files are locatable
    for i in range(E):
        assert os.path.exists(f'{workdir}/z.{i}.pkl'), f'Could not find training file {workdir}/z.{i}.pkl'
    for epoch in epochs:
        assert os.path.exists(f'{workdir}/weights.{epoch}.pkl'), f'Could not find training file {workdir}/weights.{epoch}.pkl'
    assert os.path.exists(config), f'Could not find training file {config}'
    assert os.path.exists(logfile), f'Could not find training file {logfile}'

    # Configure output paths
    if E == -1:
        outdir = f'{workdir}/convergence'
    if args.outdir:
        outdir = args.outdir
    else:
        outdir = f'{workdir}/convergence.{E}'
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(f'{outdir}/plots', exist_ok=True)
    os.makedirs(f'{outdir}/umaps', exist_ok=True)
    os.makedirs(f'{outdir}/repr_particles', exist_ok=True)
    LOG = f'{outdir}/convergence.log'
    flog(args, LOG)
    if len(epochs) < 3:
        flog('WARNING: Too few epochs have been selected for some analyses. Try decreasing --epoch-interval to a shorter interval, or analyzing a later epoch.', LOG)
    if len(epochs) < 2:
        flog('WARNING: Too few epochs have been selected for any analyses. Try decreasing --epoch-interval to a shorter interval, or analyzing a later epoch.', LOG)
        sys.exit()
    flog(f'Saving all results to {outdir}', LOG)

    # Get total number of particles, latent space dimensionality, input image size
    n_particles_total, n_dim = utils.load_pkl(f'{workdir}/z.{E}.pkl').shape
    cfg  = utils.load_pkl(config)
    img_size = cfg['lattice_args']['D'] - 1

    # Commonly used variables
    #plt.rcParams.update({'font.size': 16})
    plt.rcParams.update({'axes.linewidth': 1.5})
    chimerax_colors = np.divide(((192, 192, 192),
                                 (255, 255, 178),
                                 (178, 255, 255),
                                 (178, 178, 255),
                                 (255, 178, 255),
                                 (255, 178, 178),
                                 (178, 255, 178),
                                 (229, 191, 153),
                                 (153, 191, 229),
                                 (204, 204, 153)), 255)

    # Convergence 1: total loss
    flog('Convergence 1: plotting total loss curve ...', LOG)
    plot_loss(logfile, outdir, E, LOG)

    # Convergence 2: UMAP latent embeddings
    if args.skip_umap:
        flog('Skipping UMAP calculation ...', LOG)
    else:
        flog(f'Convergence 2: calculating and plotting UMAP embeddings of epochs {epochs} ...', LOG)
        if 'cuml.manifold.umap' in sys.modules:
            use_umap_gpu = True
        else:
            use_umap_gpu = False
        if args.force_umap_cpu:
            use_umap_gpu = False
        if use_umap_gpu:
            flog('Using GPU-accelerated UMAP via cuML library', LOG)
        else:
            flog('Using CPU-bound UMAP via umap-learn library', LOG)
        subset = args.subset
        random_state = args.random_state
        random_seed = args.random_seed
        n_epochs_umap = args.n_epochs_umap
        encoder_latent_umaps(workdir, outdir, epochs, n_particles_total, subset, random_seed, use_umap_gpu, random_state, n_epochs_umap, LOG)

    # Convergence 3: latent encoding shifts
    flog(f'Convergence 3: calculating and plotting latent encoding vector shifts for all epochs up to epoch {E} ...', LOG)
    encoder_latent_shifts(workdir, outdir, epochs, E, LOG)

    # Convergence 4: correlation of generated volumes
    flog(f'Convergence 4: sketching epoch {E}\'s latent space to find representative and well-supported latent encodings  ...', LOG)
    n_bins = args.n_bins
    smooth = args.smooth
    smooth_width = args.smooth_width
    pruned_maxima = args.pruned_maxima
    radius = args.radius
    final_maxima = args.final_maxima
    binned_ptcls_mask, labels = sketch_via_umap_local_maxima(outdir, E, LOG, n_bins=n_bins, smooth=smooth, smooth_width=smooth_width, pruned_maxima=pruned_maxima, radius=radius, final_maxima=final_maxima)

    follow_candidate_particles(workdir, outdir, epochs, n_dim, binned_ptcls_mask, labels, LOG)

    if args.skip_volgen:
        flog('Skipping volume generation ...', LOG)
    else:
        flog(f'Generating volumes at representative latent encodings for epochs {epochs} ...', LOG)
        Apix = args.Apix
        flip = args.flip
        invert = True if args.invert else None
        downsample = args.downsample
        cuda = args.cuda
        generate_volumes(workdir, outdir, epochs, Apix, flip, invert, downsample, cuda, LOG)

        flog(f'Generating masked volumes at representative latent encodings for epochs {epochs} ...', LOG)
        thresh = args.thresh
        dilate = args.dilate
        dist = args.dist
        max_threads = min(args.max_threads, multiprocessing.cpu_count())
        flog(f'Using {max_threads} threads to parallelize masking', LOG)
        mask_volumes(outdir, epochs, labels, max_threads, LOG, Apix, thresh=thresh, dilate=dilate, dist=dist)

    flog(f'Calculating masked map-map CCs at representative latent encodings for epochs {epochs} ...', LOG)
    calculate_CCs(outdir, epochs, labels, chimerax_colors, LOG)

    flog(f'Calculating masked map-map FSCs at representative latent encodings for epochs {epochs} ...', LOG)
    if args.downsample:
        img_size = args.downsample
    calculate_FSCs(outdir, epochs, labels, img_size, chimerax_colors, LOG)

    flog(f'Finished in {dt.now() - t1}', LOG)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     epilog='Example usage: $ python analyze_convergence.py [workdir] [epoch]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    add_args(parser)
    main(parser.parse_args())
