"""
Heterogeneous NN reconstruction with hierarchical pose optimization
"""
import argparse
import os
import pickle
import sys
import logging
from datetime import datetime as dt
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.parallel import DataParallel
from typing import Union
from cryodrgn import ctf, dataset, lie_tools, utils
from cryodrgn.beta_schedule import LinearSchedule, get_beta_schedule
import cryodrgn.models.config
from cryodrgn.lattice import Lattice
from cryodrgn.losses import EquivarianceLoss
from cryodrgn.models.variational_autoencoder import unparallelize, HetOnlyVAE
from cryodrgn.pose_search import PoseSearch

logger = logging.getLogger(__name__)


def add_args(parser):
    parser.add_argument("outdir", help="experiment output location")

    parser.add_argument(
        "--no-analysis",
        action="store_true",
        help="just do the training stage",
        dest="no_anlz",
    )


def make_model(args, lattice, enc_mask, in_dim) -> HetOnlyVAE:
    return HetOnlyVAE(
        lattice,
        args.qlayers,
        args.qdim,
        args.players,
        args.pdim,
        in_dim,
        args.zdim,
        encode_mode=args.encode_mode,
        enc_mask=enc_mask,
        enc_type=args.pe_type,
        enc_dim=args.pe_dim,
        domain=args.domain,
        activation={"relu": nn.ReLU, "leaky_relu": nn.LeakyReLU}[args.activation],
        feat_sigma=args.feat_sigma,
    )


def pretrain(model, lattice, optim, minibatch, tilt, zdim):
    y, yt = minibatch
    use_tilt = yt is not None
    B = y.size(0)

    model.train()
    optim.zero_grad()

    rot = lie_tools.random_SO3(B, device=y.device)
    z = torch.randn((B, zdim), device=y.device)

    # reconstruct circle of pixels instead of whole image
    mask = lattice.get_circular_mask(lattice.D // 2)

    def gen_slice(R):
        _model = unparallelize(model)
        assert isinstance(_model, HetOnlyVAE)
        return _model.decode(lattice.coords[mask] @ R, z).view(B, -1)

    y = y.view(B, -1)[:, mask]
    if use_tilt:
        yt = yt.view(B, -1)[:, mask]
        gen_loss = 0.5 * F.mse_loss(gen_slice(rot), y) + 0.5 * F.mse_loss(
            gen_slice(tilt @ rot), yt
        )
    else:
        gen_loss = F.mse_loss(gen_slice(rot), y)

    gen_loss.backward()
    optim.step()
    return gen_loss.item()


def train(
    model: Union[DataParallel, HetOnlyVAE],
    lattice,
    ps,
    optim,
    L,
    minibatch,
    beta,
    beta_control=None,
    equivariance=None,
    enc_only=False,
    poses=None,
    ctf_params=None,
):
    y, yt = minibatch
    use_tilt = yt is not None
    use_ctf = ctf_params is not None
    B = y.size(0)
    D = lattice.D

    ctf_i = None
    if use_ctf:
        freqs = lattice.freqs2d.unsqueeze(0).expand(
            B, *lattice.freqs2d.shape
        ) / ctf_params[:, 0].view(B, 1, 1)
        ctf_i = ctf.compute_ctf(freqs, *torch.split(ctf_params[:, 1:], 1, 1)).view(
            B, D, D
        )

    # TODO: Center image?
    # We do this in pose-supervised train_vae

    # VAE inference of z
    model.train()
    optim.zero_grad()
    input_ = (y, yt) if use_tilt else (y,)
    if ctf_i is not None:
        input_ = (x * ctf_i.sign() for x in input_)  # phase flip by the ctf

    _model = unparallelize(model)
    assert isinstance(_model, HetOnlyVAE)
    z_mu, z_logvar = _model.encode(*input_)
    z = _model.reparameterize(z_mu, z_logvar)

    lamb = eq_loss = None
    if equivariance is not None:
        lamb, equivariance_loss = equivariance
        eq_loss = equivariance_loss(y, z_mu)

    # pose inference
    if poses is not None:  # use provided poses
        rot = poses[0]
        trans = poses[1]
    else:  # pose search
        model.eval()
        with torch.no_grad():
            rot, trans, _base_pose = ps.opt_theta_trans(
                y,
                z=z,
                images_tilt=None if enc_only else yt,
                ctf_i=ctf_i,
            )
        model.train()

    # reconstruct circle of pixels instead of whole image
    mask = lattice.get_circular_mask(L)

    def gen_slice(R):
        slice_ = model(lattice.coords[mask] @ R, z).view(B, -1)
        if ctf_i is not None:
            slice_ *= ctf_i.view(B, -1)[:, mask]
        return slice_

    def translate(img):
        img = lattice.translate_ht(img, trans.unsqueeze(1), mask)
        return img.view(B, -1)

    y = y.view(B, -1)[:, mask]
    if use_tilt:
        yt = yt.view(B, -1)[:, mask]
    y = translate(y)
    if use_tilt:
        yt = translate(yt)

    if use_tilt:
        gen_loss = 0.5 * F.mse_loss(gen_slice(rot), y) + 0.5 * F.mse_loss(
            gen_slice(bnb.tilt @ rot), yt  # type: ignore  # noqa: F821
        )
    else:
        gen_loss = F.mse_loss(gen_slice(rot), y)

    # latent loss
    kld = torch.mean(
        -0.5 * torch.sum(1 + z_logvar - z_mu.pow(2) - z_logvar.exp(), dim=1), dim=0
    )
    if torch.isnan(kld):
        logger.info(z_mu[0])
        logger.info(z_logvar[0])
        raise RuntimeError("KLD is nan")

    if beta_control is None:
        loss = gen_loss + beta * kld / mask.sum().float()
    else:
        loss = gen_loss + beta_control * (beta - kld) ** 2 / mask.sum().float()

    if loss is not None and eq_loss is not None:
        loss += lamb * eq_loss

    loss.backward()

    optim.step()
    save_pose = [rot.detach().cpu().numpy()]
    save_pose.append(trans.detach().cpu().numpy())
    return (
        gen_loss.item(),
        kld.item(),
        loss.item(),
        eq_loss.item() if eq_loss else None,
        save_pose,
    )


def eval_z(
    model,
    lattice,
    data,
    batch_size,
    device,
    use_tilt=False,
    ctf_params=None,
    shuffler_size=0,
):
    assert not model.training
    z_mu_all = []
    z_logvar_all = []
    data_generator = dataset.make_dataloader(
        data, batch_size=batch_size, shuffler_size=shuffler_size, shuffle=False
    )

    for minibatch in data_generator:
        ind = minibatch[-1]
        y = minibatch[0].to(device)
        yt = None
        if use_tilt:
            yt = minibatch[1].to(device)
        B = len(ind)
        D = lattice.D
        c = None
        if ctf_params is not None:
            freqs = lattice.freqs2d.unsqueeze(0).expand(
                B, *lattice.freqs2d.shape
            ) / ctf_params[ind, 0].view(B, 1, 1)
            c = ctf.compute_ctf(freqs, *torch.split(ctf_params[ind, 1:], 1, 1)).view(
                B, D, D
            )
        # if trans is not None:
        #    y = lattice.translate_ht(y.view(B,-1), trans[ind].unsqueeze(1)).view(B,D,D)
        #    if yt is not None: yt = lattice.translate_ht(yt.view(B,-1), trans[ind].unsqueeze(1)).view(B,D,D)
        input_ = (y, yt) if yt is not None else (y,)
        if c is not None:
            input_ = (x * c.sign() for x in input_)  # phase flip by the ctf
        _model = unparallelize(model)
        assert isinstance(_model, HetOnlyVAE)
        z_mu, z_logvar = _model.encode(*input_)
        z_mu_all.append(z_mu.detach().cpu().numpy())
        z_logvar_all.append(z_logvar.detach().cpu().numpy())
    z_mu_all = np.vstack(z_mu_all)
    z_logvar_all = np.vstack(z_logvar_all)
    return z_mu_all, z_logvar_all


def save_checkpoint(
    model,
    lattice,
    optim,
    epoch,
    norm,
    search_pose,
    z_mu,
    z_logvar,
    out_mrc_dir,
    out_weights,
    out_z,
    out_poses,
):
    """Save model weights, latent encoding z, and decoder volumes"""
    # save model weights
    torch.save(
        {
            "epoch": epoch,
            "model_state_dict": unparallelize(model).state_dict(),
            "optimizer_state_dict": optim.state_dict(),
            "search_pose": search_pose,
        },
        out_weights,
    )
    # save z
    with open(out_z, "wb") as f:
        pickle.dump(z_mu, f)
        pickle.dump(z_logvar, f)
    with open(out_poses, "wb") as f:
        rot, trans = search_pose
        # When saving translations, save in box units (fractional)
        if isinstance(model, DataParallel):
            _model = model.module
            assert isinstance(_model, HetOnlyVAE)
            D = _model.lattice.D
        else:
            D = model.lattice.D
        pickle.dump((rot, trans / D), f)


def save_config(args, dataset, lattice, model, out_config):
    dataset_args = dict(
        particles=args.particles,
        norm=dataset.norm,
        invert_data=args.invert_data,
        ind=args.ind,
        keepreal=args.use_real,
        window=args.window,
        window_r=args.window_r,
        datadir=args.datadir,
        ctf=args.ctf,
    )
    if args.tilt is not None:
        dataset_args["particles_tilt"] = args.tilt
    lattice_args = dict(D=lattice.D, extent=lattice.extent, ignore_DC=lattice.ignore_DC)
    model_args = dict(
        qlayers=args.qlayers,
        qdim=args.qdim,
        players=args.players,
        pdim=args.pdim,
        zdim=args.zdim,
        encode_mode=args.encode_mode,
        enc_mask=args.enc_mask,
        pe_type=args.pe_type,
        feat_sigma=args.feat_sigma,
        pe_dim=args.pe_dim,
        domain=args.domain,
        activation=args.activation,
    )
    config = dict(
        dataset_args=dataset_args, lattice_args=lattice_args, model_args=model_args
    )

    cryodrgn.models.config.save(config, out_config)


def sort_poses(poses):
    ind = [x[0] for x in poses]
    ind = np.concatenate(ind)
    rot = [x[1][0] for x in poses]
    rot = np.concatenate(rot)
    rot = rot[np.argsort(ind)]
    if len(poses[0][1]) == 2:
        trans = [x[1][1] for x in poses]
        trans = np.concatenate(trans)
        trans = trans[np.argsort(ind)]
        return (rot, trans)
    return (rot,)


def get_latest(args):
    logger.info("Detecting latest checkpoint...")
    weights = [f"{args.outdir}/weights.{i}.pkl" for i in range(args.num_epochs)]
    weights = [f for f in weights if os.path.exists(f)]
    args.load = weights[-1]
    logger.info(f"Loading {args.load}")
    i = args.load.split(".")[-2]
    args.load_poses = f"{args.outdir}/pose.{i}.pkl"
    assert os.path.exists(args.load_poses)
    logger.info(f"Loading {args.load_poses}")
    return args


def main(args):
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    t1 = dt.now()
    if args.outdir is not None and not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    logger.addHandler(logging.FileHandler(f"{args.outdir}/run.log"))

    if args.load == "latest":
        args = get_latest(args)
    logger.info(" ".join(sys.argv))
    logger.info(args)

    # set the random seed
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    # set the device
    use_cuda = torch.cuda.is_available()
    device = torch.device("cuda" if use_cuda else "cpu")
    logger.info("Use cuda {}".format(use_cuda))
    if not use_cuda:
        logger.warning("WARNING: No GPUs detected")

    # set beta schedule
    if args.beta is None:
        args.beta = 1.0 / args.zdim
    try:
        args.beta = float(args.beta)
    except ValueError:
        assert (
            args.beta_control
        ), "Need to set beta control weight for schedule {}".format(args.beta)
    beta_schedule = get_beta_schedule(args.beta)

    # load index filter
    if args.ind is not None:
        logger.info("Filtering image dataset with {}".format(args.ind))
        ind = pickle.load(open(args.ind, "rb"))
    else:
        ind = None

    # load dataset
    logger.info(f"Loading dataset from {args.particles}")
    if args.tilt is None:
        tilt = None
        args.use_real = args.encode_mode == "conv"
    else:
        assert args.encode_mode == "tilt"
        tilt = torch.tensor(utils.xrot(args.tilt_deg).astype(np.float32), device=device)

    data = dataset.ImageDataset(
        mrcfile=args.particles,
        norm=args.norm,
        invert_data=args.invert_data,
        ind=ind,
        window=args.window,
        keepreal=args.use_real,
        datadir=args.datadir,
        window_r=args.window_r,
    )

    Nimg = data.N
    D = data.D

    if args.encode_mode == "conv":
        assert D - 1 == 64, "Image size must be 64x64 for convolutional encoder"

    # load ctf
    if args.ctf is not None:
        logger.info("Loading ctf params from {}".format(args.ctf))
        ctf_params = ctf.load_ctf_for_training(D - 1, args.ctf)
        if args.ind is not None:
            ctf_params = ctf_params[ind]
        assert ctf_params.shape == (Nimg, 8), ctf_params.shape
        ctf_params = torch.tensor(ctf_params, device=device)
    else:
        ctf_params = None

    lattice = Lattice(D, extent=0.5, device=device)
    if args.enc_mask is None:
        args.enc_mask = D // 2
    if args.enc_mask > 0:
        assert args.enc_mask <= D // 2
        enc_mask = lattice.get_circular_mask(args.enc_mask)
        in_dim = enc_mask.sum()
    elif args.enc_mask == -1:
        enc_mask = None
        in_dim = D**2
    else:
        raise RuntimeError(
            "Invalid argument for encoder mask radius {}".format(args.enc_mask)
        )

    model = make_model(args, lattice, enc_mask, in_dim)
    model.to(device)
    logger.info(model)
    logger.info(
        "{} parameters in model".format(
            sum(p.numel() for p in model.parameters() if p.requires_grad)
        )
    )

    equivariance_lambda = equivariance_loss = None

    if args.equivariance:
        assert args.equivariance > 0, "Regularization weight must be positive"
        equivariance_lambda = LinearSchedule(
            0, args.equivariance, args.eq_start_it, args.eq_end_it
        )
        equivariance_loss = EquivarianceLoss(model, D)

    optim = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.wd)

    if args.load == "latest":
        args = get_latest(args)

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
            if isinstance(model, DataParallel):
                _model = model.module
                assert isinstance(_model, HetOnlyVAE)
                D = _model.lattice.D
            else:
                D = model.lattice.D
            sorted_poses = (rot, trans * D)
    else:
        start_epoch = 0

    # parallelize
    if args.multigpu and torch.cuda.device_count() > 1:
        logger.info(f"Using {torch.cuda.device_count()} GPUs!")
        args.batch_size *= torch.cuda.device_count()
        logger.info(f"Increasing batch size to {args.batch_size}")
        model = DataParallel(model)
    elif args.multigpu:
        logger.warning(
            f"WARNING: --multigpu selected, but {torch.cuda.device_count()} GPUs detected"
        )

    if args.pose_model_update_freq:
        assert not args.multigpu, "TODO"
        pose_model = make_model(args, lattice, enc_mask, in_dim)
        pose_model.to(device)
        pose_model.eval()
    else:
        pose_model = model

    # save configuration
    out_config = "{}/config.yaml".format(args.outdir)
    save_config(args, data, lattice, model, out_config)

    ps = PoseSearch(
        pose_model,
        lattice,
        args.l_start,
        args.l_end,
        tilt,
        t_extent=args.t_extent,
        t_ngrid=args.t_ngrid,
        niter=args.niter,
        nkeptposes=args.nkeptposes,
        base_healpy=args.base_healpy,
        t_xshift=args.t_xshift,
        t_yshift=args.t_yshift,
        device=device,
    )

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
            loss = pretrain(model, lattice, optim, batch, tilt=ps.tilt, zdim=args.zdim)
            if global_it % args.log_interval == 0:
                logger.info(f"[Pretrain Iteration {global_it}] loss={loss:4f}")
            if global_it > args.pretrain:
                break

    # reset model after pretraining
    if args.reset_optim_after_pretrain:
        logger.info(">> Resetting optim after pretrain")
        optim = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.wd)

    # training loop
    num_epochs = args.num_epochs
    cc = 0
    if args.pose_model_update_freq:
        pose_model.load_state_dict(model.state_dict())

    epoch = None
    for epoch in range(start_epoch, num_epochs):
        t2 = dt.now()
        kld_accum = 0
        gen_loss_accum = 0
        loss_accum = 0
        eq_loss_accum = 0
        batch_it = 0
        poses = []

        L_model = lattice.D // 2
        if args.l_ramp_epochs > 0:
            Lramp = args.l_start + int(
                epoch / args.l_ramp_epochs * (args.l_end - args.l_start)
            )
            ps.Lmin = min(Lramp, args.l_start)
            ps.Lmax = min(Lramp, args.l_end)
            if epoch < args.l_ramp_epochs and args.l_ramp_model:
                L_model = ps.Lmax

        if args.reset_model_every and (epoch - 1) % args.reset_model_every == 0:
            logger.info(">> Resetting model")
            model = make_model(args, lattice, enc_mask, in_dim)

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
            global_it = Nimg * epoch + batch_it

            lamb = None
            beta = beta_schedule(global_it)
            if equivariance_lambda is not None and equivariance_loss is not None:
                lamb = equivariance_lambda(global_it)
                equivariance_tuple = (lamb, equivariance_loss)
            else:
                equivariance_tuple = None

            # train the model
            p = None
            if epoch % args.ps_freq != 0:
                p = [torch.tensor(x[ind_np], device=device) for x in sorted_poses]  # type: ignore

            cc += len(batch[0])
            if args.pose_model_update_freq and cc > args.pose_model_update_freq:
                pose_model.load_state_dict(model.state_dict())
                cc = 0

            ctf_i = ctf_params[ind] if ctf_params is not None else None
            gen_loss, kld, loss, eq_loss, pose = train(
                model,
                lattice,
                ps,
                optim,
                L_model,
                batch,
                beta,
                args.beta_control,
                equivariance_tuple,
                enc_only=args.enc_only,
                poses=p,
                ctf_params=ctf_i,
            )
            # logging
            poses.append((ind.cpu().numpy(), pose))
            kld_accum += kld * len(ind)
            gen_loss_accum += gen_loss * len(ind)
            if args.equivariance:
                assert eq_loss is not None
                eq_loss_accum += eq_loss * len(ind)

            loss_accum += loss * len(ind)
            if batch_it % args.log_interval == 0:
                eq_log = (
                    f"equivariance={eq_loss:.4f}, lambda={lamb:.4f}, "
                    if eq_loss is not None and lamb is not None
                    else ""
                )
                logger.info(
                    f"# [Train Epoch: {epoch+1}/{num_epochs}] [{batch_it}/{Nimg} images] gen loss={gen_loss:.4f}, "
                    f"kld={kld:.4f}, beta={beta:.4f}, {eq_log}loss={loss:.4f}"
                )

        eq_log = (
            "equivariance = {:.4f}, ".format(eq_loss_accum / Nimg)
            if args.equivariance
            else ""
        )
        logger.info(
            "# =====> Epoch: {} Average gen loss = {:.4}, KLD = {:.4f}, {}total loss = {:.4f}; Finished in {}".format(
                epoch + 1,
                gen_loss_accum / Nimg,
                kld_accum / Nimg,
                eq_log,
                loss_accum / Nimg,
                dt.now() - t2,
            )
        )

        sorted_poses = sort_poses(poses) if poses else None

        # save checkpoint
        if args.checkpoint and epoch % args.checkpoint == 0:
            out_mrc = "{}/reconstruct.{}.mrc".format(args.outdir, epoch)
            out_weights = "{}/weights.{}.pkl".format(args.outdir, epoch)
            out_poses = "{}/pose.{}.pkl".format(args.outdir, epoch)
            out_z = "{}/z.{}.pkl".format(args.outdir, epoch)
            model.eval()
            with torch.no_grad():
                z_mu, z_logvar = eval_z(
                    model,
                    lattice,
                    data,
                    args.batch_size,
                    device,
                    use_tilt=tilt is not None,
                    ctf_params=ctf_params,
                    shuffler_size=args.shuffler_size,
                )
                save_checkpoint(
                    model,
                    lattice,
                    optim,
                    epoch,
                    data.norm,
                    sorted_poses,
                    z_mu,
                    z_logvar,
                    out_mrc,
                    out_weights,
                    out_z,
                    out_poses,
                )

    if epoch is not None:
        # save model weights and evaluate the model on 3D lattice
        model.eval()

        out_mrc = "{}/reconstruct".format(args.outdir)
        out_weights = "{}/weights.pkl".format(args.outdir)
        out_poses = "{}/pose.pkl".format(args.outdir)
        out_z = "{}/z.pkl".format(args.outdir)
        with torch.no_grad():
            z_mu, z_logvar = eval_z(
                model,
                lattice,
                data,
                args.batch_size,
                device,
                use_tilt=tilt is not None,
                ctf_params=ctf_params,
            )
            save_checkpoint(
                model,
                lattice,
                optim,
                epoch,
                data.norm,
                sorted_poses,
                z_mu,
                z_logvar,
                out_mrc,
                out_weights,
                out_z,
                out_poses,
            )

        td = dt.now() - t1
        logger.info(
            "Finished in {} ({} per epoch)".format(td, td / (num_epochs - start_epoch))
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = add_args(parser).parse_args()
    main(args)
