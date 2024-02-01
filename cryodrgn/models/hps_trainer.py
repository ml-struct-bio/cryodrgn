import os
import pickle
import logging
from datetime import datetime as dt
from collections import OrderedDict
from typing import Any

import numpy as np
import torch
from torch import nn
from torch.nn.parallel import DataParallel
import torch.nn.functional as F

from cryodrgn import ctf, dataset, lie_tools, utils
from cryodrgn.beta_schedule import LinearSchedule, get_beta_schedule
import cryodrgn.models.config
from cryodrgn.lattice import Lattice
from cryodrgn.losses import EquivarianceLoss
from cryodrgn.models.variational_autoencoder import unparallelize, HetOnlyVAE
from cryodrgn.models.neural_nets import get_decoder
from cryodrgn.models.config import ModelConfigurations
from cryodrgn.pose_search import PoseSearch


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


class HierarchicalSearchTrainer:
    def make_model(self) -> HetOnlyVAE:
        if self.zdim == 0:
            model = get_decoder(
                3,
                self.D,
                self.configs.dec_layers,
                self.configs.dec_dim,
                self.configs.domain,
                self.configs.pe_type,
                enc_dim=self.configs.pe_dim,
                activation=self.activation,
                feat_sigma=self.configs.feat_sigma,
            )
        else:
            model = HetOnlyVAE(
                self.lattice,
                self.configs.enc_layers,
                self.configs.enc_dim,
                self.configs.dec_layers,
                self.configs.dec_dim,
                self.in_dim,
                self.zdim,
                encode_mode=self.configs.encode_mode,
                enc_mask=self.enc_mask,
                enc_type=self.configs.pe_type,
                enc_dim=self.configs.pe_dim,
                domain=self.configs.domain,
                activation=self.activation,
                feat_sigma=self.configs.feat_sigma,
            )

        return model

    def __init__(self, configs: dict) -> None:
        self.logger = logging.getLogger(__name__)
        if configs["verbose"]:
            self.logger.setLevel(logging.DEBUG)

        self.configs = AbInitioConfigurations(configs)
        self.outdir = self.configs.outdir
        os.makedirs(self.outdir, exist_ok=True)
        self.logger.addHandler(
            logging.FileHandler(os.path.join(self.outdir, "training.log"))
        )

        if self.configs.load == "latest":
            configs = self.get_latest()
            self.configs = AbInitioConfigurations(configs)

        self.logger.info(str(self.configs))

        # set the random seed
        np.random.seed(self.configs.seed)
        torch.manual_seed(self.configs.seed)

        # set the device
        use_cuda = torch.cuda.is_available()
        self.device = torch.device("cuda" if use_cuda else "cpu")
        self.logger.info("Use cuda {}".format(use_cuda))
        if not use_cuda:
            self.logger.warning("WARNING: No GPUs detected")

        # set beta schedule
        self.zdim = self.configs.z_dim
        if self.zdim > 0:
            if self.configs.beta is None:
                self.configs.beta = 1.0 / self.zdim

            self.beta_schedule = get_beta_schedule(self.configs.beta)

        # load index filter
        if self.configs.ind is not None:
            self.logger.info(f"Filtering image dataset with {self.configs.ind}")
            self.ind = pickle.load(open(self.configs.ind, "rb"))
        else:
            self.ind = None

        # load dataset
        self.logger.info(f"Loading dataset from {self.configs.particles}")
        if self.configs.tilt is None:
            self.tilt = None
        else:
            self.tilt = torch.tensor(
                utils.xrot(self.configs.tilt_deg).astype(np.float32), device=self.device
            )

        self.data = dataset.ImageDataset(
            mrcfile=self.configs.particles,
            norm=self.configs.data_norm,
            invert_data=self.configs.invert_data,
            ind=self.ind,
            window=self.configs.window,
            keepreal=self.configs.use_real,
            datadir=self.configs.datadir,
            window_r=self.configs.window_r,
        )

        self.Nimg = self.data.N
        self.D = self.data.D

        if self.configs.encode_mode == "conv":
            if self.D - 1 != 64:
                raise ValueError("Image size must be 64x64 for convolutional encoder!")

        # load ctf
        if self.configs.ctf is not None:
            self.logger.info(f"Loading ctf params from {self.configs.ctf}")
            ctf_params = ctf.load_ctf_for_training(self.D - 1, self.configs.ctf)
            if self.ind is not None:
                ctf_params = ctf_params[self.ind]
            assert ctf_params.shape == (self.Nimg, 8), ctf_params.shape
            self.ctf_params = torch.tensor(ctf_params, device=self.device)
        else:
            self.ctf_params = None

        self.lattice = Lattice(self.D, extent=0.5, device=self.device)

        if self.zdim > 0:
            if self.configs.enc_mask is None:
                self.configs.enc_mask = self.D // 2

            if self.configs.enc_mask > 0:
                assert self.configs.enc_mask <= self.D // 2
                self.enc_mask = self.lattice.get_circular_mask(self.configs.enc_mask)
                self.in_dim = self.enc_mask.sum()
            elif self.configs.enc_mask == -1:
                self.enc_mask = None
                self.in_dim = self.D**2
            else:
                raise RuntimeError(
                    f"Invalid argument for encoder mask radius {self.configs.enc_mask}"
                )

        else:
            self.enc_mask = None
            self.in_dim = None

        self.activation = {"relu": nn.ReLU, "leaky_relu": nn.LeakyReLU}[
            self.configs.activation
        ]
        self.model = self.make_model()
        self.model.to(self.device)

        self.logger.info(self.model)
        self.logger.info(
            "{} parameters in model".format(
                sum(p.numel() for p in self.model.parameters() if p.requires_grad)
            )
        )

        self.equivariance_lambda = self.equivariance_loss = None
        if self.zdim > 0:
            if self.configs.equivariance:
                self.equivariance_lambda = LinearSchedule(
                    0,
                    self.configs.equivariance,
                    self.configs.equivariance_start,
                    self.configs.equivariance_stop,
                )
                self.equivariance_loss = EquivarianceLoss(self.model, self.D)

        self.optim = torch.optim.Adam(
            self.model.parameters(),
            lr=self.configs.learning_rate,
            weight_decay=self.configs.weight_decay,
        )

        self.sorted_poses = list()
        if self.configs.load:
            self.configs.pretrain = 0
            self.logger.info("Loading checkpoint from {}".format(self.configs.load))
            checkpoint = torch.load(self.configs.load)
            self.model.load_state_dict(checkpoint["model_state_dict"])
            self.optim.load_state_dict(checkpoint["optimizer_state_dict"])
            self.start_epoch = checkpoint["epoch"] + 1
            self.model.train()

            if self.configs.load_poses:
                rot, trans = utils.load_pkl(self.configs.load_poses)
                assert np.all(
                    trans <= 1
                ), "ERROR: Old pose format detected. Translations must be in units of fraction of box."

                # Convert translations to pixel units to feed back to the model
                if isinstance(self.model, DataParallel):
                    _model = self.model.module
                    assert isinstance(_model, HetOnlyVAE)
                    D = _model.lattice.D
                else:
                    D = self.model.lattice.D

                self.sorted_poses = (rot, trans * D)
        else:
            self.start_epoch = 0

        # parallelize
        if self.zdim > 0:
            if self.configs.multigpu and torch.cuda.device_count() > 1:
                self.logger.info(f"Using {torch.cuda.device_count()} GPUs!")
                self.configs.batch_size *= torch.cuda.device_count()
                self.logger.info(f"Increasing batch size to {self.configs.batch_size}")
                self.model = DataParallel(self.model)
            elif self.configs.multigpu:
                self.logger.warning(
                    f"WARNING: --multigpu selected, but {torch.cuda.device_count()} GPUs detected"
                )
        else:
            if self.configs.multigpu:
                raise NotImplementedError("--multigpu")

        if self.configs.pose_model_update_freq:
            assert not self.configs.multigpu, "TODO"
            self.pose_model = self.make_model()
            self.pose_model.to(self.device)
            self.pose_model.eval()
        else:
            self.pose_model = self.model

        # save configuration
        self.configs.write(os.path.join(self.outdir, "train-configs.yaml"))

        self.pose_search = PoseSearch(
            self.pose_model,
            self.lattice,
            self.configs.l_start,
            self.configs.l_end,
            self.configs.tilt,
            t_extent=self.configs.t_extent,
            t_ngrid=self.configs.t_ngrid,
            niter=self.configs.grid_niter,
            nkeptposes=self.configs.n_kept_poses,
            base_healpy=self.configs.base_healpy,
            t_xshift=self.configs.t_xshift,
            t_yshift=self.configs.t_yshift,
            device=self.device,
        )

        self.data_iterator = dataset.make_dataloader(
            self.data,
            batch_size=self.configs.batch_size,
            shuffler_size=self.configs.shuffler_size,
        )

    def train(self):

        t0 = dt.now()

        # pretrain decoder with random poses
        global_it = 0
        self.logger.info(f"Using random poses for {self.configs.pretrain} iterations")
        while global_it < self.configs.pretrain:
            for batch in self.data_iterator:
                global_it += len(batch[0])
                batch = (
                    (batch[0].to(self.device), None)
                    if self.configs.tilt is None
                    else (batch[0].to(self.device), batch[1].to(self.device))
                )

                loss = self.pretrain(batch)
                if global_it % self.configs.log_interval == 0:
                    self.logger.info(f"[Pretrain Iteration {global_it}] loss={loss:4f}")
                if global_it > self.configs.pretrain:
                    break

        # reset model after pretraining
        if self.configs.reset_optim_after_pretrain:
            self.logger.info(">> Resetting optim after pretrain")
            self.optim = torch.optim.Adam(
                self.model.parameters(),
                lr=self.configs.learning_rate,
                weight_decay=self.configs.weight_decay,
            )

        # training loop
        num_epochs = self.configs.num_epochs
        cc = 0
        if self.configs.pose_model_update_freq:
            self.pose_model.load_state_dict(self.model.state_dict())

        epoch = None
        for epoch in range(self.start_epoch, num_epochs):
            t2 = dt.now()
            kld_accum = 0
            gen_loss_accum = 0
            loss_accum = 0
            eq_loss_accum = 0
            batch_it = 0
            poses = []

            L_model = self.lattice.D // 2
            if self.configs.l_ramp_epochs > 0:
                Lramp = self.configs.l_start + int(
                    epoch
                    / self.configs.l_ramp_epochs
                    * (self.configs.l_end - self.configs.l_start)
                )

                self.pose_search.Lmin = min(Lramp, self.configs.l_start)
                self.pose_search.Lmax = min(Lramp, self.configs.l_end)
                if epoch < self.configs.l_ramp_epochs and self.configs.l_ramp_model:
                    L_model = self.pose_search.Lmax

            if (
                self.configs.reset_model_every
                and (epoch - 1) % self.configs.reset_model_every == 0
            ):
                self.logger.info(">> Resetting model")
                self.model = self.make_model()

            if (
                self.configs.reset_optim_every
                and (epoch - 1) % self.configs.reset_optim_every == 0
            ):
                self.logger.info(">> Resetting optim")
                self.optim = torch.optim.Adam(
                    self.model.parameters(),
                    lr=self.configs.learning_rate,
                    weight_decay=self.configs.weight_decay,
                )

            if epoch % self.configs.ps_freq != 0:
                self.logger.info("Using previous iteration poses")
            for batch in self.data_iterator:
                ind = batch[-1]
                ind_np = ind.cpu().numpy()
                batch = (
                    (batch[0].to(self.device), None)
                    if self.tilt is None
                    else (batch[0].to(self.device), batch[1].to(self.device))
                )
                batch_it += len(batch[0])
                global_it = self.Nimg * epoch + batch_it

                lamb = None
                beta = self.beta_schedule(global_it)
                if (
                    self.equivariance_lambda is not None
                    and self.equivariance_loss is not None
                ):
                    lamb = self.equivariance_lambda(global_it)
                    equivariance_tuple = (lamb, self.equivariance_loss)
                else:
                    equivariance_tuple = None

                # train the model
                p = None
                if epoch % self.configs.ps_freq != 0:
                    p = [torch.tensor(x[ind_np], device=self.device) for x in self.sorted_poses]  # type: ignore

                cc += len(batch[0])
                if (
                    self.configs.pose_model_update_freq
                    and cc > self.configs.pose_model_update_freq
                ):
                    self.pose_model.load_state_dict(self.model.state_dict())
                    cc = 0

                ctf_i = self.ctf_params[ind] if self.ctf_params is not None else None
                gen_loss, kld, loss, eq_loss, pose = self.train_step(
                    batch,
                    L_model,
                    beta,
                    equivariance_tuple,
                    poses=p,
                    ctf_params=ctf_i,
                )

                # logging
                poses.append((ind.cpu().numpy(), pose))
                kld_accum += kld * len(ind)
                gen_loss_accum += gen_loss * len(ind)
                if self.configs.equivariance:
                    assert eq_loss is not None
                    eq_loss_accum += eq_loss * len(ind)

                loss_accum += loss * len(ind)
                if batch_it % self.configs.log_interval == 0:
                    eq_log = (
                        f"equivariance={eq_loss:.4f}, lambda={lamb:.4f}, "
                        if eq_loss is not None and lamb is not None
                        else ""
                    )
                    self.logger.info(
                        f"# [Train Epoch: {epoch + 1}/{num_epochs}] "
                        f"[{batch_it}/{self.Nimg} images] gen loss={gen_loss:.4f}, "
                        f"kld={kld:.4f}, beta={beta:.4f}, {eq_log}loss={loss:.4f}"
                    )

            eq_log = (
                "equivariance = {:.4f}, ".format(eq_loss_accum / self.Nimg)
                if self.configs.equivariance
                else ""
            )
            self.logger.info(
                "# =====> Epoch: {} Average gen loss = {:.4}, KLD = {:.4f}, {}total loss = {:.4f}; Finished in {}".format(
                    epoch + 1,
                    gen_loss_accum / self.Nimg,
                    kld_accum / self.Nimg,
                    eq_log,
                    loss_accum / self.Nimg,
                    dt.now() - t2,
                )
            )

            if poses:
                ind = [x[0] for x in poses]
                ind = np.concatenate(ind)
                rot = [x[1][0] for x in poses]
                rot = np.concatenate(rot)
                rot = rot[np.argsort(ind)]

                if len(poses[0][1]) == 2:
                    trans = [x[1][1] for x in poses]
                    trans = np.concatenate(trans)
                    trans = trans[np.argsort(ind)]
                    self.sorted_poses = (rot, trans)
                else:
                    self.sorted_poses = (rot,)

            else:
                self.sorted_poses = tuple()

            # save checkpoint
            save_epoch = (
                epoch == self.configs.num_epochs - 1
                or self.configs.checkpoint
                and epoch % self.configs.checkpoint == 0
            )

            if save_epoch:
                out_mrc = os.path.join(self.outdir, f"reconstruct.{epoch}.mrc")
                out_weights = os.path.join(self.outdir, f"weights.{epoch}.pkl")
                out_poses = os.path.join(self.outdir, f"pose.{epoch}.pkl")
                out_z = os.path.join(self.outdir, f"z.{epoch}.pkl")

                self.model.eval()
                with torch.no_grad():
                    self.save_checkpoint(
                        epoch,
                        out_mrc,
                        out_weights,
                        out_z,
                        out_poses,
                    )

            td = dt.now() - t0
            self.logger.info(
                f"Finished in {td} ({td / (num_epochs - self.start_epoch)} per epoch)"
            )

    def train_step(
        self,
        minibatch,
        L,
        beta,
        equivariance=None,
        poses=None,
        ctf_params=None,
    ):
        y, yt = minibatch
        use_tilt = yt is not None
        use_ctf = ctf_params is not None
        B = y.size(0)
        D = self.lattice.D

        ctf_i = None
        if use_ctf:
            freqs = self.lattice.freqs2d.unsqueeze(0).expand(
                B, *self.lattice.freqs2d.shape
            ) / ctf_params[:, 0].view(B, 1, 1)
            ctf_i = ctf.compute_ctf(freqs, *torch.split(ctf_params[:, 1:], 1, 1)).view(
                B, D, D
            )

        # TODO: Center image?
        # We do this in pose-supervised train_vae

        # VAE inference of z
        self.model.train()
        self.optim.zero_grad()
        input_ = (y, yt) if use_tilt else (y,)
        if ctf_i is not None:
            input_ = (x * ctf_i.sign() for x in input_)  # phase flip by the ctf

        _model = unparallelize(self.model)
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
            self.model.eval()
            with torch.no_grad():
                rot, trans, _base_pose = self.pose_search.opt_theta_trans(
                    y,
                    z=z,
                    images_tilt=None if self.configs.enc_only else yt,
                    ctf_i=ctf_i,
                )
            self.model.train()

        # reconstruct circle of pixels instead of whole image
        mask = self.lattice.get_circular_mask(L)

        def gen_slice(R):
            slice_ = self.model(self.lattice.coords[mask] @ R, z).view(B, -1)
            if ctf_i is not None:
                slice_ *= ctf_i.view(B, -1)[:, mask]
            return slice_

        def translate(img):
            img = self.lattice.translate_ht(img, trans.unsqueeze(1), mask)
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
            self.logger.info(z_mu[0])
            self.logger.info(z_logvar[0])
            raise RuntimeError("KLD is nan")

        if self.configs.beta_control is None:
            loss = gen_loss + beta * kld / mask.sum().float()
        else:
            loss = self.configs.beta_control * (beta - kld) ** 2 / mask.sum().float()
            loss += gen_loss

        if loss is not None and eq_loss is not None:
            loss += lamb * eq_loss

        loss.backward()

        self.optim.step()
        save_pose = [rot.detach().cpu().numpy()]
        save_pose.append(trans.detach().cpu().numpy())
        return (
            gen_loss.item(),
            kld.item(),
            loss.item(),
            eq_loss.item() if eq_loss else None,
            save_pose,
        )

    def pretrain(self, batch):
        if self.zdim > 0:
            y, yt = batch
            use_tilt = yt is not None
            B = y.size(0)

            self.model.train()
            self.optim.zero_grad()

            rot = lie_tools.random_SO3(B, device=y.device)
            z = torch.randn((B, self.zdim), device=y.device)

            # reconstruct circle of pixels instead of whole image
            mask = self.lattice.get_circular_mask(self.lattice.D // 2)

            def gen_slice(R):
                _model = unparallelize(self.model)
                assert isinstance(_model, HetOnlyVAE)
                return _model.decode(self.lattice.coords[mask] @ R, z).view(B, -1)

            y = y.view(B, -1)[:, mask]
            if use_tilt:
                yt = yt.view(B, -1)[:, mask]
                gen_loss = 0.5 * F.mse_loss(gen_slice(rot), y) + 0.5 * F.mse_loss(
                    gen_slice(self.configs.tilt @ rot), yt
                )
            else:
                gen_loss = F.mse_loss(gen_slice(rot), y)

            gen_loss.backward()
            self.optim.step()

            return gen_loss.item()

        else:
            y, yt = batch
            B = y.size(0)
            self.model.train()
            self.optim.zero_grad()

            mask = self.lattice.get_circular_mask(self.lattice.D // 2)

            def gen_slice(R):
                slice_ = self.model(self.lattice.coords[mask] @ R)
                return slice_.view(B, -1)

            rot = lie_tools.random_SO3(B, device=y.device)

            y = y.view(B, -1)[:, mask]
            if self.tilt is not None:
                yt = yt.view(B, -1)[:, mask]
                loss = 0.5 * F.mse_loss(gen_slice(rot), y) + 0.5 * F.mse_loss(
                    gen_slice(self.tilt @ rot), yt
                )
            else:
                loss = F.mse_loss(gen_slice(rot), y)
            loss.backward()
            self.optim.step()

            return loss.item()

    def save_checkpoint(
        self,
        epoch,
        search_pose,
        out_weights,
        out_z,
        out_poses,
    ):
        """Save model weights, latent encoding z, and decoder volumes"""
        # save model weights
        torch.save(
            {
                "epoch": epoch,
                "model_state_dict": unparallelize(self.model).state_dict(),
                "optimizer_state_dict": self.optim.state_dict(),
                "search_pose": search_pose,
            },
            out_weights,
        )

        if self.zdim > 0:
            assert not self.model.training
            z_mu_all = []
            z_logvar_all = []
            data_generator = dataset.make_dataloader(
                self.data,
                batch_size=self.configs.batch_size,
                shuffler_size=self.configs.shuffler_size,
                shuffle=False,
            )

            for minibatch in data_generator:
                ind = minibatch[-1]
                y = minibatch[0].to(self.device)
                yt = None
                if self.configs.tilt:
                    yt = minibatch[1].to(self.device)
                B = len(ind)
                D = self.lattice.D
                c = None

                if self.ctf_params is not None:
                    freqs = self.lattice.freqs2d.unsqueeze(0).expand(
                        B, *self.lattice.freqs2d.shape
                    ) / self.ctf_params[ind, 0].view(B, 1, 1)
                    c = ctf.compute_ctf(
                        freqs, *torch.split(self.ctf_params[ind, 1:], 1, 1)
                    ).view(B, D, D)

                input_ = (y, yt) if yt is not None else (y,)
                if c is not None:
                    input_ = (x * c.sign() for x in input_)  # phase flip by the ctf

                _model = unparallelize(self.model)
                assert isinstance(_model, HetOnlyVAE)
                z_mu, z_logvar = _model.encode(*input_)
                z_mu_all.append(z_mu.detach().cpu().numpy())
                z_logvar_all.append(z_logvar.detach().cpu().numpy())

            # save z
            z_mu_all, z_logvar_all = np.vstack(z_mu_all), np.vstack(z_logvar_all)
            with open(out_z, "wb") as f:
                pickle.dump(z_mu_all, f)
                pickle.dump(z_logvar_all, f)

        with open(out_poses, "wb") as f:
            rot, trans = search_pose
            # When saving translations, save in box units (fractional)
            if isinstance(self.model, DataParallel):
                _model = self.model.module
                assert isinstance(_model, HetOnlyVAE)
                D = _model.lattice.D
            else:
                D = self.model.lattice.D
            pickle.dump((rot, trans / D), f)

    def get_latest(self) -> None:
        self.logger.info("Detecting latest checkpoint...")

        weights = [
            os.path.join(self.configs.outdir, f"weights.{epoch}.pkl")
            for epoch in range(self.configs.num_epochs)
        ]
        weights = [f for f in weights if os.path.exists(f)]

        self.configs.load = weights[-1]
        self.logger.info(f"Loading {self.configs.load}")
        epoch = self.configs.load.split(".")[-2]
        self.configs.load_poses = os.path.join(self.configs.outdir, f"pose.{epoch}.pkl")
        assert os.path.exists(self.configs.load_poses)
        self.logger.info(f"Loading {self.configs.load_poses}")


class AbInitioConfigurations(ModelConfigurations):

    __slots__ = (
        "particles",
        "ctf",
        "dataset",
        "datadir",
        "ind",
        "log_interval",
        "verbose",
        "load",
        "load_poses",
        "checkpoint",
        "z_dim",
        "invert_data",
        "window",
        "window_r",
        "lazy",
        "shuffler_size",
        "max_threads",
        "tilt",
        "tilt_deg",
        "enc_only",
        "num_epochs",
        "batch_size",
        "weight_decay",
        "learning_rate",
        "beta",
        "beta_control",
        "equivariance",
        "equivariance_start",
        "equivariance_stop",
        "data_norm",
        "l_ramp_epochs",
        "l_ramp_model",
        "reset_model_every",
        "reset_optim_every",
        "reset_optim_after_pretrain",
        "multigpu",
        "l_start",
        "l_end",
        "grid_niter",
        "t_extent",
        "t_ngrid",
        "t_xshift",
        "t_yshift",
        "pretrain",
        "ps_freq",
        "n_kept_poses",
        "base_healpy",
        "pose_model_update_freq",
        "enc_layers",
        "enc_dim",
        "encode_mode",
        "enc_mask",
        "use_real",
        "dec_layers",
        "dec_dim",
        "pe_type",
        "feat_sigma",
        "pe_dim",
        "domain",
        "activation",
        "seed",
    )
    defaults = OrderedDict(
        {
            "particles": None,
            "ctf": None,
            "dataset": None,
            "datadir": None,
            "ind": None,
            "log_interval": 1000,
            "verbose": False,
            "load": None,
            "load_poses": None,
            "checkpoint": 1,
            "z_dim": None,
            "invert_data": True,
            "window": True,
            "window_r": 0.85,
            "lazy": False,
            "shuffler_size": 0,
            "max_threads": 16,
            "tilt": None,
            "tilt_deg": 45,
            "enc_only": False,
            "num_epochs": 30,
            "batch_size": 8,
            "weight_decay": 0,
            "learning_rate": 1e-4,
            "beta": None,
            "beta_control": None,
            "equivariance": None,
            "equivariance_start": 100000,
            "equivariance_stop": 200000,
            "data_norm": None,
            "l_ramp_epochs": 0,
            "l_ramp_model": 0,
            "reset_model_every": None,
            "reset_optim_every": None,
            "reset_optim_after_pretrain": None,
            "multigpu": False,
            "l_start": 12,
            "l_end": 32,
            "grid_niter": 4,
            "t_extent": 10,
            "t_ngrid": 7,
            "t_xshift": 0,
            "t_yshift": 0,
            "pretrain": 10000,
            "ps_freq": 5,
            "n_kept_poses": 8,
            "base_healpy": 2,
            "pose_model_update_freq": None,
            "enc_layers": 3,
            "enc_dim": 256,
            "encode_mode": "resid",
            "enc_mask": None,
            "use_real": False,
            "dec_layers": 3,
            "dec_dim": 256,
            "pe_type": "gaussian",
            "feat_sigma": 0.5,
            "pe_dim": None,
            "domain": "hartley",
            "activation": "relu",
            "seed": None,
        }
    )

    def __init__(self, config_vals: dict[str, Any]) -> None:
        assert config_vals["model"] == "hps"
        super().__init__(config_vals)

        if self.dataset is None:
            if self.particles is None:
                raise ValueError(
                    "As dataset was not specified, please " "specify particles!"
                )
            if self.ctf is None:
                raise ValueError("As dataset was not specified, please " "specify ctf!")

        if self.beta is not None:
            if self.z_dim == 0:
                raise ValueError("Cannot use beta with homogeneous reconstruction!.")

            if not isinstance(self.beta, (int, float)) and not self.beta_control:
                raise ValueError(
                    f"Need to set beta control weight for schedule {self.beta}"
                )

        if self.tilt is None:
            if self.z_dim > 0:
                if self.use_real != (self.encode_mode == "conv"):
                    raise ValueError(
                        "Using real space image is only available "
                        "for convolutional encoder in SPA heterogeneous reconstruction!"
                    )
        else:
            if self.z_dim > 0:
                if self.encode_mode != "tilt":
                    raise ValueError(
                        "Must use tilt for heterogeneous reconstruction on ET capture!"
                    )

        if self.equivariance is not None:
            if self.z_dim == 0:
                raise ValueError(
                    "Cannot use equivariance with homogeneous reconstruction!."
                )
            if self.equivariance <= 0:
                raise ValueError("Regularization weight must be positive")
