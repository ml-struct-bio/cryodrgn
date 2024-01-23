"""Summaries produced by models in the course of training."""

import matplotlib.pyplot as plt
import numpy as np
import torch
import logging
from matplotlib import colors
import seaborn as sns

from cryodrgn import analysis
from cryodrgn import lie_tools
from cryodrgn import utils
from cryodrgn import fft

logger = logging.getLogger(__name__)


def make_scalar_summary(writer, scalars, step):
    for key in scalars:
        writer.add_scalar(key, scalars[key], step)


def make_conf_summary(
    writer, conf, steps, labels, prefix="", pca=None, logvar=None, palette_type=None
):
    if conf.shape[1] > 1:
        if pca is None:
            pc, pca = analysis.run_pca(conf)
        else:
            pc = pca.transform(conf)

        fig = plt.figure(dpi=96, figsize=(5, 5))
        plt.hexbin(pc[:, 0], pc[:, 1], cmap="Blues", gridsize=35)
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        writer.add_figure(prefix + "Conformation - PCA", fig, global_step=steps)

        if palette_type is not None:
            n_colors = 50
            if palette_type == "rainbow":
                palette = sns.color_palette("hls", n_colors)
            else:
                # linear
                palette = sns.color_palette("blend:#0098FF,#FF009C", n_colors=n_colors)

            rainbow_labels = np.repeat(
                np.arange(n_colors)[:, None], int(len(pc) / n_colors) + 1, axis=-1
            ).reshape(-1)[: len(pc)]
            f = plt.figure(figsize=(5, 5), dpi=96)
            sns.scatterplot(
                x=pc[:, 0], y=pc[:, 1], hue=rainbow_labels, palette=palette, alpha=0.01
            )
            plt.xlabel("PC1")
            plt.ylabel("PC2")
            plt.legend([], [], frameon=False)
            writer.add_figure(
                prefix + "Conformation - Rainbow PCA", f, global_step=steps
            )

        if labels is not None:
            # colored PCA
            centers = []
            n_classes = int(np.max(labels) - np.min(labels)) + 1
            for i in range(n_classes):
                z_sub = conf[labels == i]
                centers.append(np.mean(z_sub, axis=0))
            centers = np.array(centers)
            centers2, centers_i = analysis.get_nearest_point(conf, centers)

            fig, _ = analysis.plot_by_cluster(
                pc[:, 0],
                pc[:, 1],
                n_classes,
                labels,
                centers_ind=centers_i,
                annotate=True,
                dpi=96,
                figsize=(5, 5),
                s=5,
            )
            writer.add_figure(
                prefix + "Conformation - Colored PCA", fig, global_step=steps
            )

    if conf.shape[1] == 1:
        if labels is not None:
            fig = plt.figure(dpi=96, figsize=(6, 4))
            i_min = int(np.min(labels))
            i_max = int(np.max(labels))
            for i in range(i_min, i_max + 1):
                indices_at_i = (labels > (i - 0.5)) & (labels < (i + 0.5))
                plt.hist(
                    conf[indices_at_i].flatten(), bins=50, histtype="step", label=str(i)
                )
            plt.legend(loc="best")
            plt.grid(True)
            writer.add_figure(
                prefix + "Conformation - Colored Histogram", fig, global_step=steps
            )
        else:
            fig = plt.figure(dpi=96, figsize=(6, 4))
            plt.plot(conf[:, 0], "k.")
            writer.add_figure(
                prefix + "Conformation - Histogram", fig, global_step=steps
            )

    if logvar is not None:
        std_max = np.max(np.exp(0.5 * logvar), -1)
        std_min = np.min(np.exp(0.5 * logvar), -1)

        fig = plt.figure(dpi=96, figsize=(7, 5))
        plt.hist(std_max.flatten(), bins=100, alpha=0.5, color="red", label="max")
        plt.hist(std_min.flatten(), bins=100, alpha=0.5, color="blue", label="min")
        plt.xlabel("Std. Dev.")
        plt.legend(loc="best")
        writer.add_figure(
            prefix + "Conformation - Histogram of Std. Dev.", fig, global_step=steps
        )

    return pca


def make_img_summary(writer, in_dict, recon_y, output_mask, steps, prefix=""):
    """
    writer: Writer
    in_dict: dict
        y: [batch_size(, n_tilts), D, D]
        y_real: [batch_size(, n_tilts), low_res, low_res]
    recon_y: [batch_size(, n_tilts), n_pts]
    output_mask: [D**2]
    steps: int
    prefix: str
    """
    if in_dict["y"].dim() == 3:
        y = in_dict["y"][0]
        recon_y = recon_y[0]
        y_real = in_dict["y_real"][0]
    else:
        y = in_dict["y"][0, 0]
        recon_y = recon_y[0, 0]
        y_real = in_dict["y_real"][0, 0]

    # HT
    fig = plt.figure(dpi=96, figsize=(10, 4))
    # GT
    plt.subplot(121)
    y_np = utils.to_numpy(y)
    plt.imshow(y_np, cmap="plasma")
    plt.colorbar()
    plt.axis("off")
    plt.title(f"Ground Truth ({y_np.shape[-2]}x{y_np.shape[-1]})")

    # Pred
    plt.subplot(122)
    y_pred = torch.zeros_like(y).reshape(-1)
    y_pred[output_mask.binary_mask] = recon_y.to(dtype=y_pred.dtype)
    y_pred = utils.to_numpy(y_pred)
    resolution = int(np.sqrt(y_pred.shape[0]))
    y_pred = y_pred.reshape(resolution, resolution)
    plt.imshow(y_pred, cmap="plasma")
    plt.colorbar()
    plt.axis("off")
    plt.title(f"Prediction ({resolution}x{resolution})")
    writer.add_figure(prefix + "Hartley Transform", fig, global_step=steps)

    # HT (log)
    fig = plt.figure(dpi=96, figsize=(10, 4))
    # GT
    plt.subplot(121)
    plt.imshow(np.abs(y_np) + 1e-9, norm=colors.LogNorm(), cmap="RdPu")
    plt.colorbar()
    plt.axis("off")
    plt.title(f"Ground Truth ({y_np.shape[-2]}x{y_np.shape[-1]})")

    # Pred
    plt.subplot(122)
    plt.imshow(np.abs(y_pred) + 1e-9, norm=colors.LogNorm(), cmap="RdPu")
    plt.colorbar()
    plt.axis("off")
    plt.title(f"Prediction ({resolution}x{resolution})")
    writer.add_figure(prefix + "Log Hartley Transform", fig, global_step=steps)

    # Image
    fig = plt.figure(dpi=96, figsize=(13, 4))
    # GT
    plt.subplot(131)
    y_real_np = utils.to_numpy(y_real)
    plt.imshow(y_real_np)
    plt.colorbar()
    plt.axis("off")
    plt.title(f"Input Encoder ({y_real_np.shape[-2]}x{y_real_np.shape[-1]})")

    plt.subplot(132)
    y_real_np = fft.ihtn_center(y_np)
    plt.imshow(y_real_np)
    plt.colorbar()
    plt.axis("off")
    plt.title(f"Ground Truth ({y_real_np.shape[-2]}x{y_real_np.shape[-1]})")

    # Pred
    plt.subplot(133)
    y_real_pred = fft.ihtn_center(y_pred)
    plt.imshow(y_real_pred.reshape(resolution, resolution))
    plt.colorbar()
    plt.axis("off")
    plt.title(f"Prediction ({resolution}x{resolution})")
    writer.add_figure(prefix + "Image", fig, global_step=steps)


def single_loss_summary(writer, per_image_loss, angles, steps, suffix):
    """
    writer: Writer
    per_image_loss: [n_img, 1]
    angles: [n_img]
    steps: int
    suffix: str
    """
    losses = per_image_loss.flatten()
    x_min = np.percentile(losses, 1)
    x_max = np.percentile(losses, 99)
    x_90 = np.percentile(losses, 90)
    bins = np.linspace(x_min, x_max, 50)

    fig = plt.figure(dpi=96, figsize=(7, 5))
    plt.hist(losses, bins=bins, histtype="step")
    plt.yscale("log")
    plt.xlabel("L2 Loss " + suffix)
    plt.ylabel("Number of Images")
    plt.grid(which="both")
    plt.xlim(x_min, x_max)
    writer.add_figure("Log Loss Histogram " + suffix, fig, global_step=steps)

    if angles is not None:
        fig = plt.figure(dpi=96, figsize=(5, 5))
        plt.hexbin(losses, angles, cmap="Blues", gridsize=35)
        plt.xlabel("L2 Loss " + suffix)
        plt.ylabel("Out-of-Plane Error")
        writer.add_figure(
            "Out-of-Plane Error vs. Loss" + suffix, fig, global_step=steps
        )

        fig = plt.figure(dpi=96, figsize=(7, 5))
        plt.hist(angles[losses > x_90], bins=100, histtype="step")
        plt.xlabel("Out-of-Plane Error (deg)")
        writer.add_figure(
            "Out-of-Plane Errors 10% Junk" + suffix, fig, global_step=steps
        )


def make_loss_summary(writer, per_image_losses, angles, steps):
    """
    writer: Writer
    per_image_losses: dict
        full: [n_img, 1]
        lf: [n_img, 1]
        hf: [n_img, 1]
    angles: [n_img] (numpy)
    steps: int
    """
    for key in per_image_losses.keys():
        suffix = "(" + key + ")"
        single_loss_summary(writer, per_image_losses[key], angles, steps, suffix)


def make_pose_summary(
    writer,
    rots,
    trans,
    rotmat_gt,
    trans_gt,
    steps,
    shift=False,
    prefix="",
    view_direction_only=False,
):
    """
    writer: Writer
    rots: [n_img, 3, 3] (numpy)
    trans: [n_img, 2] (numpy)
    rotmat_gt: [n_img, 3, 3] (tensor)
    trans_gt: [n_img, 2] (tensor)
    steps: int
    shift: bool
    prefix: str
    view_direction_only: bool

    output: [n_img] (numpy)
    """
    rots = torch.tensor(rots).float()

    unitvec_gt = torch.tensor([0, 0, 1], dtype=torch.float32).reshape(3, 1)
    # plot out-of-plane angles (not aligned)
    out_of_planes = torch.sum(rots * unitvec_gt, dim=-2)
    out_of_planes = out_of_planes.numpy()
    out_of_planes /= np.linalg.norm(out_of_planes, axis=-1, keepdims=True)
    azimuth, elevation = lie_tools.direction_to_azimuth_elevation(out_of_planes)

    fig = plt.figure(dpi=96, figsize=(10, 4))
    plt.subplot(111)
    plt.hist2d(
        azimuth.flatten(),
        elevation.flatten(),
        cmap="Reds",
        bins=100,
        range=[[-np.pi, np.pi], [-np.pi / 2.0, np.pi / 2.0]],
    )

    plt.grid(True)
    plt.xlabel("Azimuth (rad.)")
    plt.ylabel("Elevation (rad.)")
    writer.add_figure(prefix + "Out-of-Plane Distribution", fig, global_step=steps)

    if rotmat_gt is not None:
        if not view_direction_only:
            rot_pred_aligned, rot_gt_aligned, median_frob = align_rot_best(
                rotmat_gt, rots
            )
            rot_pred_aligned = torch.tensor(rot_pred_aligned).float()
            rot_gt_aligned = torch.tensor(rot_gt_aligned).float()
            # writer.add_scalar(prefix + 'Median Frobenius Rotations', median_frob, steps)
        else:
            rots, rot_pred_aligned, rot_gt_aligned = align_view_dir(rotmat_gt, rots)

        angles, mean_oop_error, median_oop_error = get_angular_error(
            rotmat_gt, rot_pred_aligned
        )
        writer.add_scalar(
            prefix + "Out-of-Plane Angle - Mean Error (deg.)", mean_oop_error, steps
        )
        writer.add_scalar(
            prefix + "Out-of-Plane Angle - Median Error (deg.)", median_oop_error, steps
        )

        fig = plt.figure(dpi=96, figsize=(7, 5))
        plt.hist(angles.flatten(), bins=100, histtype="step")
        plt.xlabel("Error (deg.)")
        writer.add_figure(
            prefix + "Out-of-Plane Angle - Histogram of Errors", fig, global_step=steps
        )

        euler_gt = lie_tools.rotmat_to_euler(rotmat_gt.numpy())
        euler_pred = lie_tools.rotmat_to_euler(torch.clone(rots).numpy())
        euler_pred_aligned = lie_tools.rotmat_to_euler(rot_pred_aligned.numpy())
        euler_gt_aligned = lie_tools.rotmat_to_euler(rot_gt_aligned.numpy())

        fig = plt.figure(dpi=96, figsize=(15, 5))
        plt.subplot(131)
        plt.hexbin(
            euler_gt[:, 0], euler_pred_aligned[:, 0], cmap="Oranges", gridsize=35
        )
        plt.xlabel("alpha (gt)")
        plt.ylabel("alpha (pred)")

        plt.subplot(132)
        plt.hexbin(euler_gt[:, 1], euler_pred_aligned[:, 1], cmap="Greens", gridsize=35)
        plt.xlabel("beta (gt)")
        plt.ylabel("beta (pred)")

        plt.subplot(133)
        plt.hexbin(
            euler_gt[:, 2], euler_pred_aligned[:, 2], cmap="Purples", gridsize=35
        )
        plt.xlabel("gamma (gt)")
        plt.ylabel("gamma (pred)")

        writer.add_figure(
            prefix + "Euler Angles - Prediction vs. Ground Truth",
            fig,
            global_step=steps,
        )

        fig = plt.figure(dpi=96, figsize=(15, 5))
        plt.subplot(131)
        plt.hexbin(
            euler_pred[:, 0], euler_gt_aligned[:, 0], cmap="Oranges", gridsize=35
        )
        plt.xlabel("alpha (pred)")
        plt.ylabel("alpha (gt)")
        plt.subplot(132)
        plt.hexbin(euler_pred[:, 1], euler_gt_aligned[:, 1], cmap="Greens", gridsize=35)
        plt.xlabel("beta (pred)")
        plt.ylabel("beta (gt)")
        plt.subplot(133)
        plt.hexbin(
            euler_pred[:, 2], euler_gt_aligned[:, 2], cmap="Purples", gridsize=35
        )
        plt.xlabel("gamma (pred)")
        plt.ylabel("gamma (gt)")
        writer.add_figure(
            prefix + "Euler Angles - Ground Truth vs. Prediction",
            fig,
            global_step=steps,
        )

        ip_angular_errors = np.min(
            np.concatenate(
                [
                    np.abs(euler_pred_aligned[:, 0] - euler_gt[:, 0])[None],
                    np.abs(euler_pred_aligned[:, 0] - (euler_gt[:, 0] - 2.0 * np.pi))[
                        None
                    ],
                    np.abs(euler_pred_aligned[:, 0] - (euler_gt[:, 0] + 2.0 * np.pi))[
                        None
                    ],
                ],
                0,
            ),
            0,
        )
        ip_angular_errors = ip_angular_errors * 180.0 / np.pi
        mean_ip_error = np.mean(ip_angular_errors)
        median_ip_error = np.median(ip_angular_errors)
        writer.add_scalar(
            prefix + "In-Plane Angle - Mean Error (deg.)", mean_ip_error, steps
        )
        writer.add_scalar(
            prefix + "In-Plane Angle - Median Error (deg.)", median_ip_error, steps
        )

        fig = plt.figure(dpi=96, figsize=(7, 5))
        plt.hist(ip_angular_errors.flatten(), bins=100, histtype="step")
        plt.xlabel("Error (deg.)")
        writer.add_figure(
            prefix + "In-Plane Angle - Histogram of Errors", fig, global_step=steps
        )

        if shift:
            trans = torch.tensor(trans).float()
            (
                trans_pred_corr,
                trans_gt_corr,
                trans_metrics,
                dist,
            ) = get_trans_metrics(trans_gt, trans, rots, correct_global_trans=True)
            writer.add_scalar(
                prefix + "Translation - RMSE (pix.)", trans_metrics["rmse"], steps
            )
            writer.add_scalar(
                prefix + "Translation - RMedSE (pix.)", trans_metrics["rmedse"], steps
            )
            writer.add_scalar(
                prefix + "Translation - Mean Absolute Error (pix.)",
                trans_metrics["mae"],
                steps,
            )
            writer.add_scalar(
                prefix + "Translation - Median Absolute Error (pix.)",
                trans_metrics["medae"],
                steps,
            )

            fig = plt.figure(dpi=96, figsize=(7, 5))
            plt.hist(np.sqrt(dist).flatten(), bins=100, histtype="step")
            plt.xlabel("Error (pix.)")
            writer.add_figure(
                prefix + "Translation - Histogram of Errors", fig, global_step=steps
            )

            fig = plt.figure(dpi=96, figsize=(10, 4))
            plt.subplot(121)
            plt.hexbin(
                trans_gt_corr[:, 0],
                trans_pred_corr[:, 0],
                cmap="Blues",
                gridsize=35,
                extent=(-5, 5, -5, 5),
            )
            plt.xlabel("x gt (pix.)")
            plt.ylabel("x pred (pix.)")
            plt.subplot(122)
            plt.hexbin(
                trans_gt_corr[:, 1],
                trans_pred_corr[:, 1],
                cmap="Blues",
                gridsize=35,
                extent=(-5, 5, -5, 5),
            )
            plt.xlabel("y gt (pix.)")
            plt.ylabel("y pred (pix.)")
            writer.add_figure(
                prefix + "Translation - Prediction vs. Ground Truth",
                fig,
                global_step=steps,
            )


def get_ref_matrix(r1, r2, i, flip=False):
    if flip:
        return np.matmul(r2[i].T, _flip(r1[i]))
    else:
        return np.matmul(r2[i].T, r1[i])


def _flip(rot):
    x = np.diag([1, 1, -1]).astype(rot.dtype)
    return np.matmul(x, rot)


def align_rot(r1, r2, i, flip=False):
    if flip:
        return np.matmul(_flip(r2), get_ref_matrix(r1, r2, i, flip=True))
    else:
        return np.matmul(r2, get_ref_matrix(r1, r2, i, flip=False))


def rigid_transform_3d(a, b):
    """
    a: [batch_size, 3] (numpy)
    b: [batch_size, 3] (numpy)

    output: [3, 3] (numpy)
    """
    a = np.transpose(a)
    b = np.transpose(b)

    H = a @ np.transpose(b)

    # find rotation
    u, s, vt = np.linalg.svd(H)
    r = vt.T @ u.T

    return r


def align_view_dir(rot_gt_tensor, rot_pred_tensor):
    """
    rot_gt_tensor: [n_rots, 3, 3]
    rot_pred_tensor: [n_rots, 3, 3]

    output: [n_rots, 3, 3], [n_rots, 3, 3], [n_rots, 3, 3]
    """
    rot_pred = rot_pred_tensor.clone().numpy()
    rot_gt = rot_gt_tensor.clone().numpy()

    alignment = torch.tensor(
        rigid_transform_3d(
            lie_tools.rotmat_to_r3(rot_pred), lie_tools.rotmat_to_r3(rot_gt)
        )
    ).float()

    rotmat_pred_aligned = alignment @ rot_pred_tensor
    rotmat_gt_aligned = alignment.permute(0, 1) @ rot_gt_tensor.float()

    euler_gt = lie_tools.rotmat_to_euler(rot_gt)
    euler_pred_aligned = lie_tools.rotmat_to_euler(rotmat_pred_aligned.clone().numpy())
    euler_pred_aligned[..., 0] = euler_gt[..., 0]
    rotmat_pred_aligned = torch.tensor(
        lie_tools.euler_to_rotmat(euler_pred_aligned)
    ).float()

    rotmat_pred = alignment.permute(0, 1) @ rotmat_pred_aligned

    return rotmat_pred, rotmat_pred_aligned, rotmat_gt_aligned


def align_rot_best(rot_gt_tensor, rot_pred_tensor, n_tries=100):
    """
    rot_gt_tensor: [n_rots, 3, 3]
    rot_pred_tensor: [n_rots, 3, 3]
    n_tries: int

    output: [n_rots, 3, 3] (numpy), [n_rots, 3, 3] (numpy), float
    """
    n_tries = min(n_tries, rot_gt_tensor.shape[0])
    rot_gt = rot_gt_tensor.clone().numpy()
    rot_pred = rot_pred_tensor.clone().numpy()

    median = []
    for i in range(n_tries):
        rot_pred_aligned = align_rot(rot_gt, rot_pred, i, flip=False)
        dists = frob_norm(rot_gt, rot_pred_aligned)
        median.append(np.median(dists))

    median_flip = []
    for i in range(n_tries):
        rot_pred_aligned = align_rot(rot_gt, rot_pred, i, flip=True)
        dists = frob_norm(rot_gt, rot_pred_aligned)
        median_flip.append(np.median(dists))

    if np.min(median) < np.min(median_flip):
        logger.info("Correct Handedness")
        i_best = np.argmin(median)
        alignment_matrix = get_ref_matrix(rot_gt, rot_pred, i_best, flip=False)
        rot_pred_aligned = np.matmul(rot_pred, alignment_matrix)
        rot_gt_aligned = np.matmul(rot_gt, alignment_matrix.T)
        median_frob = np.min(median)

    else:
        logger.info("Flipped Handedness")
        i_best = np.argmin(median_flip)
        alignment_matrix = get_ref_matrix(rot_gt, rot_pred, i_best, flip=True)
        rot_pred_aligned = np.matmul(_flip(rot_pred), alignment_matrix)
        rot_gt_aligned = _flip(np.matmul(rot_gt, alignment_matrix.T))
        median_frob = np.min(median_flip)

    return rot_pred_aligned, rot_gt_aligned, median_frob


def frob_norm(r1, r2):
    """
    r1: [n_rots, 3, 3]
    r2: [n_rots, 3, 3]

    output: float
    """
    return np.sum((r1 - r2) ** 2, axis=(1, 2))


def get_angular_error(rot_gt, rot_pred):
    """
    rot_gt: [n_rots, 3, 3]
    rot_pred: [n_rots, 3, 3]

    output: [n_rots] (numpy), float, float
    """
    unitvec_gt = torch.tensor([0, 0, 1], dtype=torch.float32).reshape(3, 1)

    out_of_planes_gt = torch.sum(rot_gt * unitvec_gt, dim=-2)
    out_of_planes_gt = out_of_planes_gt.numpy()
    out_of_planes_gt /= np.linalg.norm(out_of_planes_gt, axis=-1, keepdims=True)

    out_of_planes_pred = torch.sum(rot_pred * unitvec_gt, dim=-2)
    out_of_planes_pred = out_of_planes_pred.numpy()
    out_of_planes_pred /= np.linalg.norm(out_of_planes_pred, axis=-1, keepdims=True)

    angles = (
        np.arccos(np.clip(np.sum(out_of_planes_gt * out_of_planes_pred, -1), -1.0, 1.0))
        * 180.0
        / np.pi
    )

    return angles, np.mean(angles), np.median(angles)


def get_trans_metrics(trans_gt, trans_pred, rotmat, correct_global_trans=False):
    if correct_global_trans:
        b = torch.cat([(trans_pred - trans_gt)[:, 0], (trans_pred - trans_gt)[:, 1]], 0)

        matrix_a = torch.cat([rotmat[:, 0, :], rotmat[:, 1, :]], 0)
        u = torch.tensor(np.linalg.lstsq(matrix_a, b, rcond=-1)[0]).float()
        matrix_n = torch.tensor([[1, 0, 0], [0, 1, 0]]).reshape(1, 2, 3).float()
        batch_size = rotmat.shape[0]

        trans_pred_corr = trans_pred - torch.bmm(
            matrix_n.repeat(batch_size, 1, 1),
            (u @ rotmat.permute(0, 2, 1)).reshape(-1, 3, 1),
        ).reshape(-1, 2)

    else:
        trans_pred_corr = trans_pred

    dist = np.sum((trans_gt.numpy() - trans_pred_corr.numpy()) ** 2, axis=1)

    trans_metrics = {
        "rmse": np.sqrt(np.mean(dist)),
        "rmedse": np.sqrt(np.median(dist)),
        "mae": np.mean(np.sqrt(dist)),
        "medae": np.median(np.sqrt(dist)),
    }

    return trans_pred_corr, trans_gt, trans_metrics, dist
