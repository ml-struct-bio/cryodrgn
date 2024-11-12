"""Homogeneous neural net ab initio reconstruction with hierarchical pose optimization.

Example usages
--------------
$ cryodrgn abinit_homo particles.256.txt --ctf ctf.pkl --ind chosen-particles.pkl \
                                         -o cryodrn-out/256_abinit-homo

"""
import argparse
import os
from datetime import datetime as dt
import logging
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from cryodrgn import ctf, dataset, lie_tools, models, utils
from cryodrgn.mrc import MRCFile
import cryodrgn.config

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument("outdir", help="experiment output location")

    parser.add_argument(
        "--no-analysis",
        action="store_true",
        help="Lazy loading if full dataset is too large to fit in memory",
    )
    parser.add_argument(
        "--shuffler-size",
        type=int,
        default=0,
        help="If non-zero, will use a data shuffler for faster lazy data loading.",
    )

    group = parser.add_argument_group("Tilt series")
    group.add_argument("--tilt", help="Particle stack file (.mrcs)")
    group.add_argument(
        "--tilt-deg",
        type=float,
        default=45,
        help="X-axis tilt offset in degrees (default: %(default)s)",
    )

    group = parser.add_argument_group("Training parameters")
    group.add_argument(
        "--t-extent",
        type=float,
        default=10,
        help="+/- pixels to search over translations (default: %(default)s)",
    )
    group.add_argument(
        "--t-ngrid",
        type=float,
        default=7,
        help="Initial grid size for translations (default: %(default)s)",
    )
    group.add_argument(
        "--t-xshift",
        type=float,
        default=0,
        help="X-axis translation shift (default: %(default)s)",
    )
    group.add_argument(
        "--t-yshift",
        type=float,
        default=0,
        help="Y-axis translation shift (default: %(default)s)",
    )
    group.add_argument(
        "--no-trans", action="store_true", help="Don't search over translations"
    )
    group.add_argument(
        "--pretrain",
        type=int,
        default=10000,
        help="Number of initial iterations with random poses (default: %(default)s)",
    )
    group.add_argument(
        "--ps-freq",
        type=int,
        default=5,
        help="Frequency of pose inference (default: every %(default)s epochs)",
    )
    group.add_argument(
        "-n",
        "--num-epochs",
        type=int,
        default=30,
        help="Number of training epochs (default: %(default)s)",
    )
    group.add_argument(
        "-b",
        "--batch-size",
        type=int,
        default=8,
        help="Minibatch size (default: %(default)s)",
    )
    group.add_argument(
        "--wd",
        type=float,
        default=0,
        help="Weight decay in Adam optimizer (default: %(default)s)",
    )
    group.add_argument(
        "--lr",
        type=float,
        default=1e-4,
        help="Learning rate in Adam optimizer (default: %(default)s)",
    )
    group.add_argument(
        "--reset-model-every", type=int, help="If set, reset the model every N epochs"
    )
    group.add_argument(
        "--reset-optim-every",
        type=int,
        help="If set, reset the optimizer every N epochs",
    )
    group.add_argument(
        "--reset-optim-after-pretrain",
        type=int,
        help="If set, reset the optimizer every N epochs",
    )

    group = parser.add_argument_group("Pose search parameters")
    group.add_argument(
        "--l-start",
        type=int,
        default=12,
        help="Starting L radius (default: %(default)s)",
    )
    group.add_argument(
        "--l-end", type=int, default=32, help="End L radius (default: %(default)s)"
    )
    group.add_argument(
        "--niter",
        type=int,
        default=4,
        help="Number of iterations of grid subdivision (default: %(default)s)",
    )
    group.add_argument(
        "--l-ramp-epochs",
        type=int,
        default=25,
        help="Number of epochs to ramp up to --l-end (default: %(default)s)",
    )
    group.add_argument(
        "--probabilistic", action="store_true", help="Use probabilistic bound"
    )
    group.add_argument(
        "--nkeptposes",
        type=int,
        default=8,
        help="Number of poses to keep at each refinement interation "
        "during branch and bound (default: %(default)s)",
    )
    group.add_argument(
        "--base-healpy",
        type=int,
        default=2,
        help="Base healpy grid for pose search. Higher means exponentially higher "
        "resolution (default: %(default)s)",
    )
    group.add_argument(
        "--pose-model-update-freq",
        type=int,
        help="If set, only update the model used for pose search every N examples",
    )

    group = parser.add_argument_group("Network Architecture")
    group.add_argument(
        "--layers",
        type=int,
        default=3,
        help="Number of hidden layers (default: %(default)s)",
    )
    group.add_argument(
        "--dim",
        type=int,
        default=256,
        help="Number of nodes in hidden layers (default: %(default)s)",
    )
    group.add_argument(
        "--l-extent",
        type=float,
        default=0.5,
        help="Coordinate lattice size (if not using positional encoding) "
        "(default: %(default)s)",
    )
    group.add_argument(
        "--pe-type",
        choices=(
            "geom_ft",
            "geom_full",
            "geom_lowf",
            "geom_nohighf",
            "linear_lowf",
            "gaussian",
            "none",
        ),
        default="gaussian",
        help="Type of positional encoding (default: %(default)s)",
    )
    group.add_argument(
        "--pe-dim",
        type=int,
        help="Num frequencies in positional encoding (default: D/2)",
    )
    group.add_argument(
        "--domain",
        choices=("hartley", "fourier"),
        default="hartley",
        help="Volume decoder representation (default: %(default)s)",
    )
    group.add_argument(
        "--activation",
        choices=("relu", "leaky_relu"),
        default="relu",
        help="Activation (default: %(default)s)",
    )
    group.add_argument(
        "--feat-sigma",
        type=float,
        default=0.5,
        help="Scale for random Gaussian features (default: %(default)s)",
    )


def main(args: argparse.Namespace, configs: Optional[dict[str, Any]] = None) -> None:
    if configs is None:
        configs = SetupHelper(args.outdir, update_existing=False).create_configs()

    if "z_dim" in configs and configs["z_dim"] > 0:
        configs["z_dim"] = 0

        print(
            "WARNING: given configurations specify heterogeneous reconstruction, "
            "updating them to specify homogeneous reconstruction!"
        )

    # pose inference
    if poses is not None:
        rot = poses[0]
        if not no_trans:
            trans = poses[1]
    else:  # BNB
        model.eval()
        with torch.no_grad():
            rot, trans, base_pose = ps.opt_theta_trans(
                y, images_tilt=yt, init_poses=base_pose, ctf_i=ctf_i
            )
            base_pose = base_pose.detach().cpu().numpy()

    # reconstruct circle of pixels instead of whole image
    mask = lattice.get_circular_mask(lattice.D // 2)
    # mask = lattice.get_circular_mask(ps.Lmax)

    def gen_slice(R):
        slice_ = model(lattice.coords[mask] @ R).view(B, -1)
        if ctf_i is not None:
            slice_ *= ctf_i.view(B, -1)[:, mask]
        return slice_

    def translate(img):
        img = lattice.translate_ht(img, trans.unsqueeze(1), mask)
        return img.view(B, -1)

    # Train model
    model.train()
    optim.zero_grad()

    y = y.view(B, -1)[:, mask]
    if tilt_rot is not None:
        yt = yt.view(B, -1)[:, mask]
    if not no_trans:
        y = translate(y)
        if tilt_rot is not None:
            yt = translate(yt)

    if tilt_rot is not None:
        loss = 0.5 * F.mse_loss(gen_slice(rot), y) + 0.5 * F.mse_loss(
            gen_slice(tilt_rot @ rot), yt
        )
    else:
        loss = F.mse_loss(gen_slice(rot), y)
    loss.backward()
    optim.step()
    save_pose = [rot.detach().cpu().numpy()]
    if not no_trans:
        save_pose.append(trans.detach().cpu().numpy())
    return loss.item(), save_pose, base_pose


def get_latest(args):
    # assumes args.num_epochs > latest checkpoint
    logger.info("Detecting latest checkpoint...")
    weights = [f"{args.outdir}/weights.{i}.pkl" for i in range(args.num_epochs)]
    weights = [f for f in weights if os.path.exists(f)]
    args.load = weights[-1]
    logger.info(f"Loading {args.load}")
    i = args.load.split(".")[-2]
    args.load_poses = f"{args.outdir}/pose.{i}.pkl"
    assert os.path.exists(args.load_poses)  # Might need to relax this assert
    logger.info(f"Loading {args.load_poses}")
    return args


def make_model(args, D: int):
    activation = {"relu": nn.ReLU, "leaky_relu": nn.LeakyReLU}[args.activation]
    return get_decoder(
        3,
        D,
        args.layers,
        args.dim,
        args.domain,
        args.pe_type,
        enc_dim=args.pe_dim,
        activation=activation,
        feat_sigma=args.feat_sigma,
    )


def save_config(args, dataset, lattice, out_config):
    dataset_args = dict(
        particles=args.particles,
        norm=dataset.norm,
        invert_data=args.invert_data,
        ind=args.ind,
        window=args.window,
        window_r=args.window_r,
        datadir=args.datadir,
        ctf=args.ctf,
    )
    if args.tilt is not None:
        dataset_args["particles_tilt"] = args.tilt
    lattice_args = dict(D=lattice.D, extent=lattice.extent, ignore_DC=lattice.ignore_DC)
    model_args = dict(
        qlayers=args.layers,
        qdim=args.dim,
        pe_type=args.pe_type,
        feat_sigma=args.feat_sigma,
        pe_dim=args.pe_dim,
        domain=args.domain,
        activation=args.activation,
    )
    config = dict(
        dataset_args=dataset_args, lattice_args=lattice_args, model_args=model_args
    )

    cryodrgn.config.save(config, out_config)


def main(args):

    if args.pose_model_update_freq:
        pose_model = make_model(args, D)
        pose_model.to(device)
        pose_model.eval()
    else:
        pose_model = model

    sorted_poses = []
    if args.load:
        args.pretrain = 0
        logger.info("Loading checkpoint from {}".format(args.load))
        checkpoint = torch.load(args.load)
        model.load_state_dict(checkpoint["model_state_dict"])
        optim.load_state_dict(checkpoint["optimizer_state_dict"])
        start_epoch = checkpoint["epoch"] + 1
        model.train()
        if args.load_poses:
            rot, trans = utils.load_pkl(args.load_poses)
            assert np.all(
                trans <= 1
            ), "ERROR: Old pose format detected. Translations must be in units of fraction of box."
            # Convert translations to pixel units to feed back to the model
            sorted_poses = (rot, trans * model.D)

            # sorted_base_poses = None   # TODO: need to save base_poses if we are going to use it
    else:
        start_epoch = 0

    data_iterator = dataset.make_dataloader(
        data, batch_size=args.batch_size, shuffler_size=args.shuffler_size
    )

    # pretrain decoder with random poses
    global_it = 0
    logger.info("Using random poses for {} iterations".format(args.pretrain))
    while global_it < args.pretrain:
        for batch in data_iterator:
            global_it += len(batch[0])
            batch = (
                (batch[0].to(device), None)
                if tilt is None
                else (batch[0].to(device), batch[1].to(device))
            )
            loss = pretrain(model, lattice, optim, batch, tilt=ps.tilt)
            if global_it % args.log_interval == 0:
                logger.info(f"[Pretrain Iteration {global_it}] loss={loss:4f}")
            if global_it > args.pretrain:
                break
    out_mrc = "{}/pretrain.reconstruct.mrc".format(args.outdir)
    model.eval()
    vol = model.eval_volume(lattice.coords, lattice.D, lattice.extent, tuple(data.norm))
    MRCFile.write(out_mrc, vol)

    # reset model after pretraining
    if args.reset_optim_after_pretrain:
        logger.info(">> Resetting optim after pretrain")
        optim = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.wd)

    # training loop
    cc = 0
    if args.pose_model_update_freq:
        pose_model.load_state_dict(model.state_dict())

    epoch = None
    for epoch in range(start_epoch, args.num_epochs):
        t2 = dt.now()
        batch_it = 0
        loss_accum = 0
        poses, base_poses = [], []

        if args.l_ramp_epochs > 0:
            Lramp = args.l_start + int(
                epoch / args.l_ramp_epochs * (args.l_end - args.l_start)
            )
            ps.Lmin = min(Lramp, args.l_start)
            ps.Lmax = min(Lramp, args.l_end)

        if args.reset_model_every and (epoch - 1) % args.reset_model_every == 0:
            logger.info(">> Resetting model")
            model = make_model(args, D)

        if args.reset_optim_every and (epoch - 1) % args.reset_optim_every == 0:
            logger.info(">> Resetting optim")
            optim = torch.optim.Adam(
                model.parameters(), lr=args.lr, weight_decay=args.wd
            )

        if epoch % args.ps_freq != 0:
            logger.info("Using previous iteration poses")
        for batch in data_iterator:
            ind = batch[-1]
            ind_np = ind.cpu().numpy()
            batch = (
                (batch[0].to(device), None)
                if tilt is None
                else (batch[0].to(device), batch[1].to(device))
            )
            batch_it += len(batch[0])

            # train the model
            if epoch % args.ps_freq != 0:
                p = [torch.tensor(x[ind_np], device=device) for x in sorted_poses]  # type: ignore
                # bp = sorted_base_poses[ind_np]
                bp = None
            else:
                p, bp = None, None

            cc += len(batch[0])
            if args.pose_model_update_freq and cc > args.pose_model_update_freq:
                pose_model.load_state_dict(model.state_dict())
                cc = 0

            c = ctf_params[ind] if ctf_params is not None else None
            loss_item, pose, base_pose = train(
                model,
                lattice,
                ps,
                optim,
                batch,
                tilt,
                args.no_trans,
                poses=p,
                base_pose=bp,
                ctf_params=c,
            )
            poses.append((ind.cpu().numpy(), pose))
            base_poses.append((ind_np, base_pose))
            # logging
            loss_accum += loss_item * len(batch[0])
            if batch_it % args.log_interval == 0:
                logger.info(
                    "# [Train Epoch: {}/{}] [{}/{} images] loss={:.4f}".format(
                        epoch + 1, args.num_epochs, batch_it, Nimg, loss_item
                    )
                )

        logger.info(
            "# =====> Epoch: {} Average loss = {:.4}; Finished in {}".format(
                epoch + 1, loss_accum / Nimg, dt.now() - t2
            )
        )

        # sort pose
        sorted_poses = sort_poses(poses) if poses else None
        # sorted_base_poses = sort_base_poses(base_poses)

        if args.checkpoint and epoch % args.checkpoint == 0:
            out_mrc = "{}/reconstruct.{}.mrc".format(args.outdir, epoch)
            out_weights = "{}/weights.{}.pkl".format(args.outdir, epoch)
            out_poses = "{}/pose.{}.pkl".format(args.outdir, epoch)
            save_checkpoint(
                model,
                lattice,
                sorted_poses,
                optim,
                epoch,
                data.norm,
                out_mrc,
                out_weights,
                out_poses,
            )

    if epoch is not None:
        # save model weights and evaluate the model on 3D lattice
        out_mrc = "{}/reconstruct.mrc".format(args.outdir)
        out_weights = "{}/weights.pkl".format(args.outdir)
        out_poses = "{}/pose.pkl".format(args.outdir)
        save_checkpoint(
            model,
            lattice,
            sorted_poses,
            optim,
            epoch,
            data.norm,
            out_mrc,
            out_weights,
            out_poses,
        )

        td = dt.now() - t1
        logger.info(
            "Finished in {} ({} per epoch)".format(
                td, td / (args.num_epochs - start_epoch)
            )
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_args(parser)
    main(parser.parse_args())
