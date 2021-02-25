import os
import numpy as np
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import subprocess

from scipy.spatial.distance import cdist, pdist
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture

from . import utils
log = utils.log

def parse_loss(f):
    '''Parse loss from run.log'''
    lines = open(f).readlines()
    lines = [x for x in lines if '====' in x]
    try:
        loss = [x.strip().split()[-1] for x in lines]
        loss = np.asarray(loss).astype(np.float32)
    except:
        loss = [x.split()[-4][:-1] for x in lines]
        loss = np.asarray(loss).astype(np.float32)
    return loss

### Dimensionality reduction ###

def run_pca(z):
    pca = PCA(z.shape[1])
    pca.fit(z)
    log('Explained variance ratio:')
    log(pca.explained_variance_ratio_)
    pc = pca.transform(z)
    return pc, pca

def get_pc_traj(pca, zdim, numpoints, dim, start, end, percentiles=None):
    '''
    Create trajectory along specified principle component
    
    Inputs:
        pca: sklearn PCA object from run_pca
        zdim (int)
        numpoints (int): number of points between @start and @end
        dim (int): PC dimension for the trajectory (1-based index)
        start (float): Value of PC{dim} to start trajectory
        end (float): Value of PC{dim} to stop trajectory
        percentiles (np.array or None): Define percentile array instead of np.linspace(start,stop,numpoints)
    
    Returns:
        np.array (numpoints x zdim) of z values along PC
    '''
    if percentiles is not None:
        assert len(percentiles) == numpoints
    traj_pca = np.zeros((numpoints,zdim))
    traj_pca[:,dim-1] = np.linspace(start, end, numpoints) if percentiles is None else percentiles
    ztraj_pca = pca.inverse_transform(traj_pca)
    return ztraj_pca

def run_tsne(z, n_components=2, perplexity=1000):
    if len(z) > 10000:
        log('WARNING: {} datapoints > {}. This may take awhile.'.format(len(z), 10000))
    z_embedded = TSNE(n_components=n_components, perplexity=perplexity).fit_transform(z)
    return z_embedded

def run_umap(z, **kwargs):
    import umap # CAN GET STUCK IN INFINITE IMPORT LOOP
    reducer = umap.UMAP(**kwargs)
    z_embedded = reducer.fit_transform(z)
    return z_embedded

### Clustering ###

def cluster_kmeans(z, K, on_data=True, reorder=True):
    '''
    Cluster z by K means clustering
    Returns cluster labels, cluster centers
    If reorder=True, reorders clusters according to agglomerative clustering of cluster centers
    '''
    kmeans = KMeans(n_clusters=K,
                    random_state=0,
                    max_iter=10)
    labels = kmeans.fit_predict(z)
    centers = kmeans.cluster_centers_

    if on_data:
        centers, centers_ind = get_nearest_point(z, centers)

    if reorder:
        g = sns.clustermap(centers)
        reordered = g.dendrogram_row.reordered_ind
        centers = centers[reordered]
        if on_data: centers_ind = centers_ind[reordered]
        tmp = {k:i for i,k in enumerate(reordered)}
        labels = np.array([tmp[k] for k in labels])
    return labels, centers

def cluster_gmm(z, K, on_data=True, random_state=None, **kwargs):
    '''
    Cluster z by a K-component full covariance Gaussian mixture model
    
    Inputs:
        z (Ndata x zdim np.array): Latent encodings
        K (int): Number of clusters
        on_data (bool): Compute cluster center as nearest point on the data manifold
        random_state (int or None): Random seed used for GMM clustering
        **kwargs: Additional keyword arguments passed to sklearn.mixture.GaussianMixture

    Returns: 
        np.array (Ndata,) of cluster labels
        np.array (K x zdim) of cluster centers
    '''
    clf = GaussianMixture(n_components=K, covariance_type='full', random_state=None, **kwargs)
    labels = clf.fit_predict(z)
    centers = clf.means_
    if on_data:
        centers, centers_ind = get_nearest_point(z, centers)
    return labels, centers

def get_nearest_point(data, query):
    '''
    Find closest point in @data to @query
    Return datapoint, index
    '''
    ind = cdist(query, data).argmin(axis=1)
    return data[ind], ind

### HELPER FUNCTIONS FOR INDEX ARRAY MANIPULATION

def convert_original_indices(ind, N_orig, orig_ind):
    '''
    Convert index array into indices into the original particle stack
    ''' # todo -- finish docstring
    return np.arange(N_orig)[orig_ind][ind]

def combine_ind(N, sel1, sel2, kind='intersection'):
    # todo -- docstring
    if kind == 'intersection':
        ind_selected = set(sel1) & set(sel2)
    elif kind == 'union':
        ind_selected = set(sel1) | set(sel2)
    else:
        raise RuntimeError(f"Mode {kind} not recognized. Choose either 'intersection' or 'union'")
    ind_selected_not = np.array(sorted(set(np.arange(N)) - ind_selected))
    ind_selected = np.array(sorted(ind_selected))
    return ind_selected, ind_selected_not

def get_ind_for_cluster(labels, selected_clusters):
    '''Return index array of the selected clusters
    
    Inputs:
        labels: np.array of cluster labels for each particle
        selected_clusters: list of cluster labels to select

    Return:
        ind_selected: np.array of particle indices with the desired cluster labels

    Example usage:
        ind_keep = get_ind_for_cluster(kmeans_labels, [0,4,6,14])
    '''
    ind_selected = np.array([i for i,label in enumerate(labels) if label in selected_clusters])
    return ind_selected


### PLOTTING ###

def _get_colors(K, cmap=None):
    if cmap is not None:
        cm = plt.get_cmap(cmap)
        colors = [cm(i/float(K)) for i in range(K)]
    else:
        colors = ['C{}'.format(i) for i in range(10)]
        colors = [colors[i%len(colors)] for i in range(K)]
    return colors
   
def scatter_annotate(x, y, centers=None, centers_ind=None, annotate=True, labels=None, alpha=.1, s=1):
    fig, ax = plt.subplots()
    plt.scatter(x, y, alpha=alpha, s=s, rasterized=True)

    # plot cluster centers
    if centers_ind is not None:
        assert centers is None
        centers = np.array([[x[i],y[i]] for i in centers_ind])
    if centers is not None:
        plt.scatter(centers[:,0], centers[:,1], c='k')
    if annotate:
        assert centers is not None
        if labels is None:
            labels = range(len(centers))
        for i in labels:
            ax.annotate(str(i), centers[i,0:2]+np.array([.1,.1]))
    return fig, ax

def scatter_annotate_hex(x, y, centers=None, centers_ind=None, annotate=True, labels=None):
    g = sns.jointplot(x=x, y=y, kind='hex')

    # plot cluster centers
    if centers_ind is not None:
        assert centers is None
        centers = np.array([[x[i],y[i]] for i in centers_ind])
    if centers is not None:
        g.ax_joint.scatter(centers[:,0], centers[:,1], color='k', edgecolor='grey')
    if annotate:
        assert centers is not None
        if labels is None:
            labels = range(len(centers))
        for i in labels:
            g.ax_joint.annotate(str(i), centers[i,0:2]+np.array([.1,.1]), color='black',
                                bbox=dict(boxstyle='square,pad=.1', ec='None', fc='1', alpha=.5))
    return g

def scatter_color(x, y, c, cmap='viridis', s=1, alpha=.1, label=None, figsize=None):
    fig, ax = plt.subplots(figsize=figsize)
    assert len(x) == len(y) == len(c)
    sc = plt.scatter(x, y, s=s, alpha=alpha, rasterized=True, cmap=cmap, c=c)
    cbar = plt.colorbar(sc)
    cbar.set_alpha(1)
    cbar.draw_all()
    if label:
        cbar.set_label(label)
    return fig, ax

def plot_by_cluster(x, y, K, labels, centers=None, centers_ind=None, annotate=False, 
                    s=2, alpha=0.1, colors=None, cmap=None, figsize=None):
    fig, ax = plt.subplots(figsize=figsize)
    if type(K) is int:
        K = list(range(K))

    if colors is None:
        colors = _get_colors(len(K), cmap)

    # scatter by cluster
    for i in K:
        ii = labels == i
        x_sub = x[ii]
        y_sub = y[ii]
        plt.scatter(x_sub, y_sub, s=s, alpha=alpha, label='cluster {}'.format(i), color=colors[i], rasterized=True)

    # plot cluster centers
    if centers_ind is not None:
        assert centers is None
        centers = np.array([[x[i],y[i]] for i in centers_ind])
    if centers is not None:
        plt.scatter(centers[:,0], centers[:,1], c='k')
    if annotate:
        assert centers is not None
        for i in K:
            ax.annotate(str(i), centers[i,0:2])
    return fig, ax

def plot_by_cluster_subplot(x, y, K, labels, 
                            s=2, alpha=.1, colors=None, cmap=None, figsize=None):
    if type(K) is int:
        K = list(range(K))
    ncol = int(np.ceil(len(K)**.5))
    nrow = int(np.ceil(len(K)/ncol))
    fig, ax = plt.subplots(ncol, nrow, sharex=True, sharey=True, figsize=(10,10))
    if colors is None:
        colors = _get_colors(len(K), cmap)
    for i in K:
        ii = labels == i
        x_sub = x[ii]
        y_sub = y[ii]
        a = ax.ravel()[i]
        a.scatter(x_sub, y_sub, s=s, alpha=alpha, rasterized=True, color=colors[i])
        a.set_title(i)
    return fig, ax

def plot_euler(theta,phi,psi,plot_psi=True):
    sns.jointplot(x=theta,y=phi,kind='hex',
              xlim=(-180,180),
              ylim=(0,180)).set_axis_labels("theta", "phi")
    if plot_psi:
        plt.figure()
        plt.hist(psi)
        plt.xlabel('psi')

def ipy_plot_interactive_annotate(df, ind, opacity=.3):
    '''Interactive plotly widget for a cryoDRGN pandas dataframe with annotated points'''
    import plotly.graph_objs as go
    from ipywidgets import interactive
    if 'labels' in df.columns:
        text = [f'Class {k}: index {i}' for i,k in zip(df.index, df.labels)] # hovertext
    else:
        text = [f'index {i}' for i in df.index] # hovertext
    xaxis, yaxis = df.columns[0], df.columns[1]
    scatter = go.Scattergl(x=df[xaxis], 
                           y=df[yaxis], 
                           mode='markers',
                           text=text,
                           marker=dict(size=2,
                                       opacity=opacity,
                                       ))
    sub = df.loc[ind]
    text = [f'{k}){i}' for i,k in zip(sub.index, sub.labels)]
    scatter2 = go.Scatter(x=sub[xaxis],
                            y=sub[yaxis],
                            mode='markers+text',
                            text=text,
                            textposition="top center",
                            textfont=dict(size=9,color='black'),
                            marker=dict(size=5,color='black'))
    f = go.FigureWidget([scatter,scatter2])
    f.update_layout(xaxis_title=xaxis, yaxis_title=yaxis)
    
    def update_axes(xaxis, yaxis, color_by, colorscale):
        scatter = f.data[0]
        scatter.x = df[xaxis]
        scatter.y = df[yaxis]
        
        scatter.marker.colorscale = colorscale
        if colorscale is None:
            scatter.marker.color = None
        else:
            scatter.marker.color = df[color_by] if color_by != 'index' else df.index
    
        scatter2 = f.data[1]
        scatter2.x = sub[xaxis]
        scatter2.y = sub[yaxis]
        with f.batch_update(): # what is this for??
            f.layout.xaxis.title = xaxis
            f.layout.yaxis.title = yaxis
        
    widget = interactive(update_axes, 
                    yaxis = df.select_dtypes('number').columns, 
                    xaxis = df.select_dtypes('number').columns,
                    color_by = df.columns,
                    colorscale = [None,'hsv','plotly3','deep','portland','picnic','armyrose'])
    return widget, f

def ipy_plot_interactive(df, opacity=.3):
    '''Interactive plotly widget for a cryoDRGN pandas dataframe'''
    import plotly.graph_objs as go
    from ipywidgets import interactive
    if 'labels' in df.columns:
        text = [f'Class {k}: index {i}' for i,k in zip(df.index, df.labels)] # hovertext
    else:
        text = [f'index {i}' for i in df.index] # hovertext
    
    xaxis, yaxis = df.columns[0], df.columns[1]
    f = go.FigureWidget([go.Scattergl(x=df[xaxis],
                                  y=df[yaxis],
                                  mode='markers',
                                  text=text,
                                  marker=dict(size=2,
                                              opacity=opacity,
                                              color=np.arange(len(df)),
                                              colorscale='hsv'
                                             ))])
    scatter = f.data[0]
    N = len(df)
    f.update_layout(xaxis_title=xaxis, yaxis_title=yaxis)
    f.layout.dragmode = 'lasso'

    def update_axes(xaxis, yaxis, color_by, colorscale):
        scatter = f.data[0]
        scatter.x = df[xaxis]
        scatter.y = df[yaxis]
        
        scatter.marker.colorscale = colorscale
        if colorscale is None:
            scatter.marker.color = None
        else:
            scatter.marker.color = df[color_by] if color_by != 'index' else df.index
        with f.batch_update(): # what is this for??
            f.layout.xaxis.title = xaxis
            f.layout.yaxis.title = yaxis
 
    widget = interactive(update_axes, 
                         yaxis=df.select_dtypes('number').columns, 
                         xaxis=df.select_dtypes('number').columns,
                         color_by = df.columns,
                         colorscale = [None,'hsv','plotly3','deep','portland','picnic','armyrose'])

    t = go.FigureWidget([go.Table(
                        header=dict(values=['index']),
                        cells=dict(values=[df.index]),
                        )])

    def selection_fn(trace, points, selector):
        t.data[0].cells.values = [df.loc[points.point_inds].index]

    scatter.on_selection(selection_fn)
    return widget, f, t

def plot_projections(imgs, labels=None):
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(10,10))
    axes = axes.ravel()
    for i in range(min(len(imgs),9)):
        axes[i].imshow(imgs[i], cmap='Greys_r') 
        axes[i].axis('off')
        if labels is not None:
            axes[i].set_title(labels[i])
    return fig, axes

def gen_volumes(weights, config, zfile, outdir, cuda=None, 
                Apix=None, flip=False, downsample=None):
    '''Call cryodrgn eval_vol to generate volumes at specified z values
    Input:
        weights (str): Path to model weights .pkl
        config (str): Path to config.pkl
        zfile (str): Path to .txt file of z values
        outdir (str): Path to output directory for volumes,
        cuda (int or None): Specify cuda device
        Apix (float or None): Apix of output volume
        flip (bool): Flag to flip chirality of output volumes
        downsample (int or None): Generate volumes at this box size
    '''
    cmd = f'cryodrgn eval_vol {weights} --config {config} --zfile {zfile} -o {outdir}'
    if Apix is not None:
        cmd += f' --Apix {Apix}'
    if flip:
        cmd += f' --flip'
    if downsample is not None:
        cmd += f' -d {downsample}'
    if cuda is not None:
        cmd = f'CUDA_VISIBLE_DEVICES={cuda} {cmd}'
    log(f'Running command:\n{cmd}')
    return subprocess.check_call(cmd, shell=True)
    
def load_dataframe(z=None, pc=None, euler=None, trans=None, labels=None, tsne=None, umap=None, **kwargs):
    '''Load results into a pandas dataframe for downstream analysis'''
    data = {}
    if umap is not None:
        data['UMAP1'] = umap[:,0]
        data['UMAP2'] = umap[:,1]
    if tsne is not None:
        data['TSNE1'] = tsne[:,0]
        data['TSNE2'] = tsne[:,1]
    if pc is not None:
        zD = pc.shape[1]
        for i in range(zD):
            data[f'PC{i+1}'] = pc[:,i]
    if labels is not None:
        data['labels'] = labels
    if euler is not None:
        data['theta'] = euler[:,0]
        data['phi'] = euler[:,1]
        data['psi'] = euler[:,2]
    if trans is not None:
        data['tx'] = trans[:,0]
        data['ty'] = trans[:,1]
    if z is not None:
        zD = z.shape[1]
        for i in range(zD):
            data[f'z{i}'] = z[:,i]
    for kk,vv in kwargs.items():
        data[kk] = vv
    df = pd.DataFrame(data=data)
    df['index'] = df.index
    return df


