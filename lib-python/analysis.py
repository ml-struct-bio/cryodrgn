import numpy as np
import pickle
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.spatial.distance import cdist, pdist
import umap
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

import utils
log = utils.log

def run_pca(z):
    pca = PCA(z.shape[1])
    pca.fit(z)
    log('Explained variance ratio:')
    log(pca.explained_variance_ratio_)
    pc = pca.transform(z)
    return pc, pca

def run_tsne(z, n_components=2, perplexity=50):
    if len(z) > 10000:
        log('WARNING: {} datapoints > {}. This may take awhile.'.format(len(z), 10000))
    z_embedded = TSNE(n_components=n_components, perplexity=perplexity).fit_transform(z)
    return z_embedded

def run_umap(z):
    reducer = umap.UMAP()
    z_embedded = reducer.fit_transform(z)
    return z_embedded

def cluster_kmeans(z, K):
    '''
    Cluster z by K means clustering
    Returns cluster labels, cluster centers
    '''
    kmeans = KMeans(n_clusters=K,
                    random_state=0,
                    max_iter=10).fit(z)
    centers = kmeans.cluster_centers_
    labels = kmeans.predict(z)
    return labels, centers

def get_nearest_point(data, query):
    '''
    Find closest point in @data to @query
    Return datapoint, index
    '''
    ind = cdist(query, data).argmin(axis=1)
    return data[ind], ind

def _get_colors(K, cmap=None):
    if cmap is not None:
        cm = plt.get_cmap(cmap)
        colors = [cm(i/float(K)) for i in range(K)]
    else:
        colors = ['C{}'.format(i) for i in range(10)]
        colors = [colors[i%len(colors)] for i in range(K)]
    return colors
   
def scatter_by_cluster(x, y, K, labels, centers=None, centers_i=None, annotate=False, s=2, alpha=0.1, cmap=None):
    fig, ax = plt.subplots()
    colors = _get_colors(K, cmap)

    # scatter by cluster
    for i in range(K):
        ii = labels == i
        x_sub = x[ii]
        y_sub = y[ii]
        plt.scatter(x_sub, y_sub, s=s, alpha=alpha, label='cluster {}'.format(i), color=colors[i], rasterized=True)

    # plot cluster centers
    if centers_i is not None:
        assert centers is None
        centers = np.array([[x[i],y[i]] for i in centers_i])
    if centers is not None:
        plt.scatter(centers[:,0], centers[:,1], c='k')
    if annotate:
        assert centers is not None
        for i in range(K):
            ax.annotate(str(i), centers[i,0:2])
    return fig, ax

def scatter_by_cluster_subplot(x, y, K, labels, s=2, alpha=.1, cmap=None):
    ncol = int(np.ceil(K**.5))
    nrow = int(np.ceil(K/ncol))
    fig, ax = plt.subplots(ncol, nrow, sharex=True, sharey=True)
    colors = _get_colors(K, cmap)
    for i in range(K):
        ii = labels == i
        x_sub = x[ii]
        y_sub = y[ii]
        a = ax.ravel()[i]
        a.scatter(x_sub, y_sub, s=s, alpha=alpha, rasterized=True, color=colors[i])
        a.set_title(i)
    return fig, ax

def load_workdir_results(workdir, e):
    '''Load results into a pandas dataframe for downstream analysis'''
    pass

def load_results(zfile, config, pose=None, labels=None, pca=None, tsne=None, umap=None):
    pass
