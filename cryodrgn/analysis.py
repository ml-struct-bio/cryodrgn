import argparse
import re
import logging
import matplotlib.pyplot as plt
from matplotlib.figure import Figure, Axes
import numpy as np
import numpy.typing as npt
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import cdist
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.mixture import GaussianMixture
from typing import Optional, Union, Tuple, List
from cryodrgn.commands import eval_vol

logger = logging.getLogger(__name__)


def parse_loss(f: str) -> np.ndarray:
    """Parse loss from run.log"""
    lines = open(f).readlines()
    lines = [x for x in lines if "====" in x]
    regex = "total\sloss\s=\s(\d.\d+)"  # type: ignore  # noqa: W605
    matches = [re.search(regex, x) for x in lines]
    loss = []
    for m in matches:
        assert m is not None
        loss.append(m.group(1))
    loss = np.asarray(loss).astype(np.float32)

    return loss


# Dimensionality reduction


def run_pca(z: np.ndarray) -> Tuple[np.ndarray, PCA]:
    pca = PCA(z.shape[1])
    pca.fit(z)
    logger.info("Explained variance ratio:")
    logger.info(pca.explained_variance_ratio_)
    pc = pca.transform(z)
    return pc, pca


def get_pc_traj(
    pca: PCA,
    zdim: int,
    numpoints: int,
    dim: int,
    start: Optional[float],
    end: Optional[float],
    percentiles: Optional[np.ndarray] = None,
) -> npt.NDArray[np.float32]:
    """
    Create trajectory along specified principal component

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
    """
    if percentiles is not None:
        assert len(percentiles) == numpoints
    traj_pca = np.zeros((numpoints, zdim))
    if percentiles is not None:
        traj_pca[:, dim - 1] = percentiles
    else:
        assert start is not None
        assert end is not None
        traj_pca[:, dim - 1] = np.linspace(start, end, numpoints)
    ztraj_pca = pca.inverse_transform(traj_pca)
    return ztraj_pca


def run_tsne(
    z: np.ndarray, n_components: int = 2, perplexity: float = 1000
) -> np.ndarray:
    if len(z) > 10000:
        logger.warning(
            "WARNING: {} datapoints > {}. This may take awhile.".format(len(z), 10000)
        )
    z_embedded = TSNE(n_components=n_components, perplexity=perplexity).fit_transform(z)
    return z_embedded


def run_umap(z: np.ndarray, **kwargs) -> np.ndarray:
    import umap  # CAN GET STUCK IN INFINITE IMPORT LOOP

    reducer = umap.UMAP(**kwargs)
    z_embedded = reducer.fit_transform(z)
    return z_embedded


# Clustering


def cluster_kmeans(
    z: np.ndarray, K: int, on_data: bool = True, reorder: bool = True
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Cluster z by K means clustering
    Returns cluster labels, cluster centers
    If reorder=True, reorders clusters according to agglomerative clustering of cluster centers
    """
    kmeans = KMeans(n_clusters=K, random_state=0, max_iter=10)
    labels = kmeans.fit_predict(z)
    centers = kmeans.cluster_centers_

    centers_ind = None
    if on_data:
        centers, centers_ind = get_nearest_point(z, centers)

    if reorder:
        g = sns.clustermap(centers)
        reordered = g.dendrogram_row.reordered_ind
        centers = centers[reordered]
        if centers_ind is not None:
            centers_ind = centers_ind[reordered]
        tmp = {k: i for i, k in enumerate(reordered)}
        labels = np.array([tmp[k] for k in labels])
    return labels, centers


def cluster_gmm(
    z,
    K: int,
    on_data: bool = True,
    random_state: Union[int, np.random.RandomState, None] = None,
    **kwargs,
) -> Tuple[np.ndarray, np.ndarray]:
    """
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
    """
    clf = GaussianMixture(
        n_components=K, covariance_type="full", random_state=random_state, **kwargs
    )
    labels = clf.fit_predict(z)
    centers = clf.means_
    if on_data:
        centers, centers_ind = get_nearest_point(z, centers)
    return labels, centers


def get_nearest_point(
    data: np.ndarray, query: np.ndarray
) -> Tuple[npt.NDArray[np.float32], np.ndarray]:
    """
    Find closest point in @data to @query
    Return datapoint, index
    """
    ind = cdist(query, data).argmin(axis=1)
    return data[ind], ind


# HELPER FUNCTIONS FOR INDEX ARRAY MANIPULATION


def convert_original_indices(
    ind: np.ndarray, N_orig: int, orig_ind: np.ndarray
) -> np.ndarray:
    """
    Convert index array into indices into the original particle stack
    """  # todo -- finish docstring
    return np.arange(N_orig)[orig_ind][ind]


def combine_ind(
    N: int, sel1: np.ndarray, sel2: np.ndarray, kind: str = "intersection"
) -> Tuple[np.ndarray, np.ndarray]:
    # todo -- docstring
    if kind == "intersection":
        ind_selected = set(sel1) & set(sel2)
    elif kind == "union":
        ind_selected = set(sel1) | set(sel2)
    else:
        raise RuntimeError(
            f"Mode {kind} not recognized. Choose either 'intersection' or 'union'"
        )
    ind_selected_not = np.array(sorted(set(np.arange(N)) - ind_selected))
    ind_selected = np.array(sorted(ind_selected))
    return ind_selected, ind_selected_not


def get_ind_for_cluster(
    labels: np.ndarray, selected_clusters: np.ndarray
) -> np.ndarray:
    """Return index array of the selected clusters

    Inputs:
        labels: np.array of cluster labels for each particle
        selected_clusters: list of cluster labels to select

    Return:
        ind_selected: np.array of particle indices with the desired cluster labels

    Example usage:
        ind_keep = get_ind_for_cluster(kmeans_labels, [0,4,6,14])
    """
    ind_selected = np.array(
        [i for i, label in enumerate(labels) if label in selected_clusters]
    )
    return ind_selected


# PLOTTING


def _get_chimerax_colors(K: int) -> List:
    colors = [
        "#b2b2b2",
        "#ffffb2",
        "#b2ffff",
        "#b2b2ff",
        "#ffb2ff",
        "#ffb2b2",
        "#b2ffb2",
        "#e5bf99",
        "#99bfe5",
        "#cccc99",
    ]
    colors = [colors[i % len(colors)] for i in range(K)]
    return colors


def _get_colors(K: int, cmap: Optional[str] = None) -> List:
    if cmap is not None:
        cm = plt.get_cmap(cmap)
        colors = [cm(i / float(K)) for i in range(K)]
    else:
        colors = ["C{}".format(i) for i in range(10)]
        colors = [colors[i % len(colors)] for i in range(K)]
    return colors


def scatter_annotate(
    x: np.ndarray,
    y: np.ndarray,
    centers: Optional[np.ndarray] = None,
    centers_ind: Optional[np.ndarray] = None,
    annotate: bool = True,
    labels: Optional[np.ndarray] = None,
    alpha: Union[float, np.ndarray, None] = 0.1,
    s: Union[float, np.ndarray, None] = 1,
    colors: Union[list, str, None] = None,
) -> Tuple[Figure, Axes]:
    fig, ax = plt.subplots(figsize=(4, 4))
    plt.scatter(x, y, alpha=alpha, s=s, rasterized=True)

    # plot cluster centers
    if centers_ind is not None:
        assert centers is None
        centers = np.array([[x[i], y[i]] for i in centers_ind])
    if centers is not None:
        if colors is None:
            colors = "k"
        plt.scatter(centers[:, 0], centers[:, 1], c=colors, edgecolor="black")
    if annotate:
        assert centers is not None
        if labels is None:
            labels = np.arange(len(centers))
        assert labels is not None
        for i in labels:
            ax.annotate(str(i), centers[i, 0:2] + np.array([0.1, 0.1]))
    return fig, ax


def scatter_annotate_hex(
    x: np.ndarray,
    y: np.ndarray,
    centers: Optional[np.ndarray] = None,
    centers_ind: Optional[np.ndarray] = None,
    annotate: bool = True,
    labels: Optional[np.ndarray] = None,
    colors: Union[list, str, None] = None,
) -> sns.JointGrid:
    g = sns.jointplot(x=x, y=y, kind="hex", height=4)

    # plot cluster centers
    if centers_ind is not None:
        assert centers is None
        centers = np.array([[x[i], y[i]] for i in centers_ind])
    if centers is not None:
        if colors is None:
            colors = "k"
        g.ax_joint.scatter(centers[:, 0], centers[:, 1], c=colors, edgecolor="black")
    if annotate:
        assert centers is not None
        if labels is None:
            labels = np.arange(len(centers))
        assert labels is not None
        for i in labels:
            g.ax_joint.annotate(
                str(i),
                centers[i, 0:2] + np.array([0.1, 0.1]),
                color="black",
                bbox=dict(boxstyle="square,pad=.1", ec="None", fc="1", alpha=0.5),
            )
    return g


def scatter_color(
    x: np.ndarray,
    y: np.ndarray,
    c: np.ndarray,
    cmap: str = "viridis",
    s=1,
    alpha: float = 0.1,
    label: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
) -> Tuple[Figure, Axes]:
    fig, ax = plt.subplots(figsize=figsize)
    assert len(x) == len(y) == len(c)
    sc = plt.scatter(x, y, s=s, alpha=alpha, rasterized=True, cmap=cmap, c=c)
    cbar = plt.colorbar(sc)
    cbar.set_alpha(1)
    cbar.draw_all()
    if label:
        cbar.set_label(label)
    return fig, ax


def plot_by_cluster(
    x,
    y,
    K,
    labels,
    centers=None,
    centers_ind=None,
    annotate=False,
    s=2,
    alpha=0.1,
    colors=None,
    cmap=None,
    figsize=None,
):
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
        plt.scatter(
            x_sub,
            y_sub,
            s=s,
            alpha=alpha,
            label="cluster {}".format(i),
            color=colors[i],
            rasterized=True,
        )

    # plot cluster centers
    if centers_ind is not None:
        assert centers is None
        centers = np.array([[x[i], y[i]] for i in centers_ind])
    if centers is not None:
        plt.scatter(centers[:, 0], centers[:, 1], c="k")
    if annotate:
        assert centers is not None
        for i in K:
            ax.annotate(str(i), centers[i, 0:2])
    return fig, ax


def plot_by_cluster_subplot(
    x, y, K, labels, s=2, alpha=0.1, colors=None, cmap=None, figsize=None
):
    if type(K) is int:
        K = list(range(K))
    ncol = int(np.ceil(len(K) ** 0.5))
    nrow = int(np.ceil(len(K) / ncol))
    fig, ax = plt.subplots(ncol, nrow, sharex=True, sharey=True, figsize=(10, 10))
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


def plot_euler(theta, phi, psi, plot_psi=True):
    sns.jointplot(
        x=theta, y=phi, kind="hex", xlim=(-180, 180), ylim=(0, 180)
    ).set_axis_labels("theta", "phi")
    if plot_psi:
        plt.figure()
        plt.hist(psi)
        plt.xlabel("psi")


def ipy_plot_interactive_annotate(df, ind, opacity=0.3):
    """Interactive plotly widget for a cryoDRGN pandas dataframe with annotated points"""
    import plotly.graph_objs as go
    from ipywidgets import interactive

    if "labels" in df.columns:
        text = [
            f"Class {k}: index {i}" for i, k in zip(df.index, df.labels)
        ]  # hovertext
    else:
        text = [f"index {i}" for i in df.index]  # hovertext
    xaxis, yaxis = df.columns[0], df.columns[1]
    scatter = go.Scattergl(
        x=df[xaxis],
        y=df[yaxis],
        mode="markers",
        text=text,
        marker=dict(
            size=2,
            opacity=opacity,
        ),
    )
    sub = df.loc[ind]
    text = [f"{k}){i}" for i, k in zip(sub.index, sub.labels)]
    scatter2 = go.Scatter(
        x=sub[xaxis],
        y=sub[yaxis],
        mode="markers+text",
        text=text,
        textposition="top center",
        textfont=dict(size=9, color="black"),
        marker=dict(size=5, color="black"),
    )
    f = go.FigureWidget([scatter, scatter2])
    f.update_layout(xaxis_title=xaxis, yaxis_title=yaxis)

    def update_axes(xaxis, yaxis, color_by, colorscale):
        scatter = f.data[0]
        scatter.x = df[xaxis]
        scatter.y = df[yaxis]

        scatter.marker.colorscale = colorscale
        if colorscale is None:
            scatter.marker.color = None
        else:
            scatter.marker.color = df[color_by] if color_by != "index" else df.index

        scatter2 = f.data[1]
        scatter2.x = sub[xaxis]
        scatter2.y = sub[yaxis]
        with f.batch_update():  # what is this for??
            f.layout.xaxis.title = xaxis
            f.layout.yaxis.title = yaxis

    widget = interactive(
        update_axes,
        yaxis=df.select_dtypes("number").columns,
        xaxis=df.select_dtypes("number").columns,
        color_by=df.columns,
        colorscale=[None, "hsv", "plotly3", "deep", "portland", "picnic", "armyrose"],
    )
    return widget, f


def ipy_plot_interactive(df, opacity=0.3):
    """Interactive plotly widget for a cryoDRGN pandas dataframe"""
    import plotly.graph_objs as go
    from ipywidgets import interactive

    if "labels" in df.columns:
        text = [
            f"Class {k}: index {i}" for i, k in zip(df.index, df.labels)
        ]  # hovertext
    else:
        text = [f"index {i}" for i in df.index]  # hovertext

    xaxis, yaxis = df.columns[0], df.columns[1]
    f = go.FigureWidget(
        [
            go.Scattergl(
                x=df[xaxis],
                y=df[yaxis],
                mode="markers",
                text=text,
                marker=dict(
                    size=2, opacity=opacity, color=np.arange(len(df)), colorscale="hsv"
                ),
            )
        ]
    )
    scatter = f.data[0]
    f.update_layout(xaxis_title=xaxis, yaxis_title=yaxis)
    f.layout.dragmode = "lasso"

    def update_axes(xaxis, yaxis, color_by, colorscale):
        scatter = f.data[0]
        scatter.x = df[xaxis]
        scatter.y = df[yaxis]

        scatter.marker.colorscale = colorscale
        if colorscale is None:
            scatter.marker.color = None
        else:
            scatter.marker.color = df[color_by] if color_by != "index" else df.index
        with f.batch_update():  # what is this for??
            f.layout.xaxis.title = xaxis
            f.layout.yaxis.title = yaxis

    widget = interactive(
        update_axes,
        yaxis=df.select_dtypes("number").columns,
        xaxis=df.select_dtypes("number").columns,
        color_by=df.columns,
        colorscale=[None, "hsv", "plotly3", "deep", "portland", "picnic", "armyrose"],
    )

    t = go.FigureWidget(
        [
            go.Table(
                header=dict(values=["index"]),
                cells=dict(values=[df.index]),
            )
        ]
    )

    def selection_fn(trace, points, selector):
        t.data[0].cells.values = [df.loc[points.point_inds].index]

    scatter.on_selection(selection_fn)
    return widget, f, t


def plot_projections(imgs, labels=None, max_imgs=25):
    if len(imgs) > max_imgs:
        imgs = imgs[:max_imgs]
    N = len(imgs)
    nrows = int(np.floor(N**0.5))
    ncols = int(np.ceil(N**0.5))
    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols, figsize=(1.5 * ncols, 1.5 * nrows)
    )
    axes = axes.ravel()
    for i in range(N):
        axes[i].imshow(imgs[i], cmap="Greys_r")
        if labels is not None:
            axes[i].set_title(labels[i])
    for i in range(nrows * ncols):
        axes[i].axis("off")
    plt.tight_layout()
    return fig, axes


def gen_volumes(
    weights,
    config,
    zfile,
    outdir,
    device=None,
    Apix=None,
    flip=False,
    downsample=None,
    invert=None,
    vol_start_index=0,
):
    """Call cryodrgn eval_vol to generate volumes at specified z values
    Input:
        weights (str): Path to model weights .pkl
        config (str): Path to config.yaml
        zfile (str): Path to .txt file of z values
        outdir (str): Path to output directory for volumes,
        device (int or None): Specify cuda device
        Apix (float or None): Apix of output volume
        flip (bool): Flag to flip chirality of output volumes
        downsample (int or None): Generate volumes at this box size
        invert (bool): Invert contrast of output volumes
        vol_start_index (int): Start index for generated volumes
    """
    args = [weights, "--config", config, "--zfile", zfile, "-o", outdir]
    if Apix is not None:
        args += ["--Apix", f"{Apix}"]
    if flip:
        args += ["--flip"]
    if downsample is not None:
        args += ["-d", f"{downsample}"]
    if invert:
        args += ["--invert"]
    if device is not None:
        args += ["--device", f"{device}"]
    if vol_start_index is not None:
        args += ["--vol-start-index", f"{vol_start_index}"]

    args = eval_vol.add_args(argparse.ArgumentParser()).parse_args(args)
    return eval_vol.main(args)


def load_dataframe(
    z=None, pc=None, euler=None, trans=None, labels=None, tsne=None, umap=None, **kwargs
):
    """Load results into a pandas dataframe for downstream analysis"""
    data = {}
    if umap is not None:
        data["UMAP1"] = umap[:, 0]
        data["UMAP2"] = umap[:, 1]
    if tsne is not None:
        data["TSNE1"] = tsne[:, 0]
        data["TSNE2"] = tsne[:, 1]
    if pc is not None:
        zD = pc.shape[1]
        for i in range(zD):
            data[f"PC{i+1}"] = pc[:, i]
    if labels is not None:
        data["labels"] = labels
    if euler is not None:
        data["theta"] = euler[:, 0]
        data["phi"] = euler[:, 1]
        data["psi"] = euler[:, 2]
    if trans is not None:
        data["tx"] = trans[:, 0]
        data["ty"] = trans[:, 1]
    if z is not None:
        zD = z.shape[1]
        for i in range(zD):
            data[f"z{i}"] = z[:, i]
    for kk, vv in kwargs.items():
        data[kk] = vv
    df = pd.DataFrame(data=data)
    df["index"] = df.index
    return df
