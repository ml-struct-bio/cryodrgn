"""The training engine for cryoDRGN reconstruction models up to and including v3."""

import os
import pickle
from collections import OrderedDict
from datetime import datetime as dt
from typing import Any, Callable

import numpy as np
import torch
from torch import nn
import torch.nn.functional as F

import cryodrgn.config
from cryodrgn import ctf, dataset
from cryodrgn.models import lie_tools
from cryodrgn.models.losses import EquivarianceLoss
from cryodrgn.models.variational_autoencoder import unparallelize, HetOnlyVAE
from cryodrgn.models.neural_nets import get_decoder
from cryodrgn.models.pose_search import PoseSearch
from cryodrgn.trainers._base import ModelTrainer, ModelConfigurations


class HierarchicalPoseSearchConfigurations(ModelConfigurations):

    __slots__ = (
        "enc_only",
        "beta",
        "beta_control",
        "equivariance",
        "equivariance_start",
        "equivariance_stop",
        "l_ramp_epochs",
        "l_ramp_model",
        "reset_model_every",
        "reset_optim_every",
        "reset_optim_after_pretrain",
        "grid_niter",
        "ps_freq",
        "n_kept_poses",
        "pose_model_update_freq",
        "enc_layers",
        "enc_dim",
        "encode_mode",
        "enc_mask",
        "dec_layers",
        "dec_dim",
    )
    default_values = OrderedDict(
        {
            "enc_only": False,
            "beta": None,
            "beta_control": None,
            "equivariance": None,
            "equivariance_start": 100000,
            "equivariance_stop": 200000,
            "l_ramp_epochs": 0,
            "l_ramp_model": 0,
            "reset_model_every": None,
            "reset_optim_every": None,
            "grid_niter": 4,
            "ps_freq": 5,
            "n_kept_poses": 8,
            "pose_model_update_freq": None,
            "enc_layers": None,
            "enc_dim": None,
            "encode_mode": "resid",
            "enc_mask": None,
            "dec_layers": None,
            "dec_dim": None,
        }
    )

    def __init__(self, config_vals: dict[str, Any]) -> None:
        super().__init__(config_vals)

        assert (
            config_vals["model"] == "hps"
        ), f"Mismatched model {config_vals['model']} for HierarchicalSearchTrainer!"

        if self.dataset is None:
            if self.particles is None:
                raise ValueError(
                    "As dataset was not specified, please " "specify particles!"
                )
            if self.ctf is None:
                raise ValueError("As dataset was not specified, please " "specify ctf!")

        if self.beta is not None:
            if not self.z_dim:
                raise ValueError("Cannot use beta with homogeneous reconstruction!.")

            if not isinstance(self.beta, (int, float)) and not self.beta_control:
                raise ValueError(
                    f"Need to set beta control weight for schedule {self.beta}"
                )

        if self.tilt is None:
            if self.z_dim and self.use_real != (self.encode_mode == "conv"):
                raise ValueError(
                    "Using real space image is only available "
                    "for convolutional encoder in SPA heterogeneous reconstruction!"
                )
        else:
            if self.z_dim and self.encode_mode:
                raise ValueError(
                    "Must use tilt for heterogeneous reconstruction on ET capture!"
                )

        if self.enc_layers is None:
            self.enc_layers = self.hidden_layers
        if self.enc_dim is None:
            self.enc_dim = self.hidden_dim
        if self.dec_layers is None:
            self.dec_layers = self.hidden_layers
        if self.dec_dim is None:
            self.dec_dim = self.hidden_dim

        if self.equivariance is not None:
            if not self.z_dim:
                raise ValueError(
                    "Cannot use equivariance with homogeneous reconstruction!."
                )
            if self.equivariance <= 0:
                raise ValueError("Regularization weight must be positive")

        if self.volume_domain is None:
            self.volume_domain = "fourier" if self.use_gt_poses else "hartley"


class HierarchicalPoseSearchTrainer(ModelTrainer):

    config_cls = HierarchicalPoseSearchConfigurations

    @property
    def mask_dimensions(self) -> tuple[torch.Tensor, int]:
        if self.configs.z_dim:
            use_mask = self.configs.enc_mask or self.resolution // 2

            if use_mask > 0:
                assert use_mask <= self.resolution // 2
                enc_mask = self.lattice.get_circular_mask(use_mask)
                in_dim = enc_mask.sum().item()
            elif use_mask == -1:
                enc_mask = None
                in_dim = self.resolution**2
            else:
                raise RuntimeError(
                    f"Invalid argument for encoder mask radius {self.configs.enc_mask}"
                )

        else:
            enc_mask = None
            in_dim = None

        return enc_mask, in_dim

    def make_volume_model(self) -> nn.Module:
        self.configs: HierarchicalPoseSearchConfigurations

        if not self.configs.z_dim:
            model = get_decoder(
                in_dim=3,
                D=self.resolution,
                layers=self.configs.dec_layers,
                dim=self.configs.dec_dim,
                domain=self.configs.volume_domain,
                enc_type=self.configs.pe_type,
                enc_dim=self.configs.pe_dim,
                activation=self.activation,
                feat_sigma=self.configs.feat_sigma,
            )
        else:
            enc_mask, in_dim = self.mask_dimensions

            model = HetOnlyVAE(
                lattice=self.lattice,
                qlayers=self.configs.enc_layers,
                qdim=self.configs.enc_dim,
                players=self.configs.dec_layers,
                pdim=self.configs.dec_dim,
                in_dim=in_dim,
                z_dim=self.configs.z_dim,
                encode_mode=self.configs.encode_mode,
                enc_mask=enc_mask,
                enc_type=self.configs.pe_type,
                enc_dim=self.configs.pe_dim,
                domain=self.configs.volume_domain,
                activation=self.activation,
                feat_sigma=self.configs.feat_sigma,
            )

            enc_params = sum(
                p.numel() for p in model.encoder.parameters() if p.requires_grad
            )
            self.logger.info(f"{enc_params} parameters in encoder")
            dec_params = sum(
                p.numel() for p in model.decoder.parameters() if p.requires_grad
            )
            self.logger.info(f"{dec_params} parameters in decoder")

        all_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
        self.logger.info(f"{all_params} parameters in model")

        return model

    @classmethod
    def create_beta_schedule(
        cls, start_x: float, end_x: float, start_y: float, end_y: float
    ) -> Callable[[float], float]:
        min_y = min(start_y, end_y)
        max_y = max(start_y, end_y)
        coef = (end_y - start_y) / (end_x - start_x)

        return lambda x: np.clip((x - start_x) * coef + start_y, min_y, max_y).item(0)

    def __init__(self, configs: dict[str, Any]) -> None:
        super().__init__(configs)
        self.configs: HierarchicalPoseSearchConfigurations

        # set beta schedule
        if self.configs.z_dim:
            beta = self.configs.beta or self.configs.z_dim**-1
            # beta = 1.0 / args.ntilts

            if isinstance(beta, float):
                self.beta_schedule = lambda x: beta

            elif beta == "a":
                self.beta_schedule = self.create_beta_schedule(0, 1000000, 0.001, 15)
            elif beta == "b":
                self.beta_schedule = self.create_beta_schedule(200000, 800000, 5, 15)
            elif beta == "c":
                self.beta_schedule = self.create_beta_schedule(200000, 800000, 5, 18)
            elif beta == "d":
                self.beta_schedule = self.create_beta_schedule(1000000, 5000000, 5, 18)
            else:
                raise RuntimeError(f"Unrecognized beta schedule {beta=}!")
        else:
            self.beta_schedule = None

        if self.configs.encode_mode == "conv":
            if self.resolution - 1 != 64:
                raise ValueError("Image size must be 64x64 for convolutional encoder!")

        in_dim = self.mask_dimensions[1]
        if in_dim is not None and in_dim % 8 != 0:
            self.logger.warning(
                f"Warning: Masked input image dimension {in_dim=} "
                "is not a mutiple of 8 -- AMP training speedup is not optimized!"
            )

        self.equivariance_lambda = self.equivariance_loss = None
        if self.configs.z_dim:
            if self.configs.equivariance:
                self.equivariance_lambda = self.create_beta_schedule(
                    self.configs.equivariance_start,
                    self.configs.equivariance_stop,
                    0,
                    self.configs.equivariance,
                )
                self.equivariance_loss = EquivarianceLoss(
                    self.volume_model, self.resolution
                )

        if self.configs.refine_gt_poses:
            self.pose_optimizer = torch.optim.SparseAdam(
                list(self.pose_tracker.parameters()), lr=self.configs.pose_learning_rate
            )
        else:
            self.pose_optimizer = None

        if self.configs.pose:
            self.pose_search = None
            self.pose_model = None

        else:
            if self.configs.pose_model_update_freq:
                assert not self.configs.multigpu, "TODO"
                self.pose_model = self.make_volume_model()
                self.pose_model.to(self.device)
                self.pose_model.eval()
            else:
                self.pose_model = self.volume_model

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

    def end_epoch(self) -> None:
        eq_log = (
            "equivariance = {:.4f}, ".format(self.accum_losses["eq"] / self.image_count)
            if self.configs.equivariance
            else ""
        )
        avg_gen_loss = self.accum_losses["gen"] / self.image_count
        kld_loss = self.accum_losses["kld"] / self.image_count
        total_loss = self.accum_losses["total"] / self.image_count
        self.logger.info(
            f"# =====> Epoch: {self.current_epoch} "
            f"Average gen loss = {avg_gen_loss:.4}, KLD = {kld_loss:.4f},"
            f" {eq_log}total loss = {total_loss:.4f}; "
            f"Finished in {dt.now() - self.epoch_start_time}"
        )

    def train_batch(self, batch: dict[str, torch.Tensor]) -> None:
        y = batch["y"]
        ind = batch["indices"]

        y = y.to(self.device)
        ind_np = ind.cpu().numpy()
        B = y.size(0)
        if "tilt_indices" in batch:
            tilt_ind = batch["tilt_indices"]
            assert all(tilt_ind >= 0), tilt_ind
            tilt_ind.to(self.device)
        else:
            tilt_ind = None

        if self.equivariance_lambda is not None and self.equivariance_loss is not None:
            lamb = self.equivariance_lambda(self.total_images_seen)
            equivariance_tuple = (lamb, self.equivariance_loss)
        else:
            equivariance_tuple = None

        if self.configs.use_real:
            assert hasattr(self.data, "particles_real")
            yr = torch.from_numpy(self.data.particles_real[ind.numpy()]).to(
                self.device
            )  # type: ignore  # PYR02
        else:
            yr = None

        if self.pose_optimizer is not None:
            self.pose_optimizer.zero_grad()

        # getting the poses
        rot = trans = None
        if self.pose_search and (self.current_epoch - 1) % self.configs.ps_freq != 0:
            self.logger.info("Using previous iteration poses")
            rot = torch.tensor(
                self.predicted_rots[ind_np].astype(np.float32), device=self.device
            )

            if not self.configs.no_trans:
                trans = torch.tensor(
                    self.predicted_trans[ind_np].astype(np.float32), device=self.device
                )

        self.conf_search_particles += B
        if self.configs.pose_model_update_freq:
            if (
                self.current_epoch == 1
                or self.conf_search_particles > self.configs.pose_model_update_freq
            ):
                self.pose_model.load_state_dict(self.volume_model.state_dict())
                self.conf_search_particles = 0

        dose_filters = None
        if self.configs.tilt:
            if self.pose_tracker and not self.pose_search:
                rot, trans = self.pose_tracker.get_pose(tilt_ind.view(-1))

            ctf_param = (
                self.ctf_params[tilt_ind.view(-1)]
                if self.ctf_params is not None
                else None
            )
            y = y.view(-1, self.resolution, self.resolution, self.resolution)
            Apix = self.ctf_params[0, 0] if self.ctf_params is not None else None

            if self.configs.dose_per_tilt is not None:
                dose_filters = self.data.get_dose_filters(tilt_ind, self.lattice, Apix)

        else:
            if self.pose_tracker and not self.pose_search:
                rot, trans = self.pose_tracker.get_pose(ind)
            ctf_param = self.ctf_params[ind] if self.ctf_params is not None else None

        use_tilt = tilt_ind is not None
        use_ctf = ctf_param is not None

        if self.beta_schedule is not None:
            self.beta = self.beta_schedule(self.total_images_seen)

        ctf_i = None
        if use_ctf:
            freqs = self.lattice.freqs2d.unsqueeze(0).expand(
                B, *self.lattice.freqs2d.shape
            ) / ctf_param[:, 0].view(B, 1, 1)
            ctf_i = ctf.compute_ctf(freqs, *torch.split(ctf_param[:, 1:], 1, 1)).view(
                B, self.resolution, self.resolution
            )

        # VAE inference of z
        self.volume_optimizer.zero_grad()
        self.volume_model.train()
        input_ = (y, tilt_ind) if use_tilt else (y,)
        if ctf_i is not None:
            input_ = (x * ctf_i.sign() for x in input_)  # phase flip by the ctf

        lamb = None
        losses = dict()
        if self.configs.z_dim > 0:
            _model = unparallelize(self.volume_model)
            assert isinstance(_model, HetOnlyVAE)
            z_mu, z_logvar = _model.encode(*input_)
            z = _model.reparameterize(z_mu, z_logvar)

            if equivariance_tuple is not None:
                lamb, equivariance_loss = equivariance_tuple
                losses["eq_loss"] = equivariance_loss(y, z_mu)
        else:
            z_mu = z_logvar = z = None

        # execute a pose search if poses not given
        if rot is None:
            self.volume_model.eval()
            with torch.no_grad():
                rot, trans, base_pose = self.pose_search.opt_theta_trans(
                    y,
                    z=z,
                    images_tilt=None if self.configs.enc_only else tilt_ind,
                    ctf_i=ctf_i,
                )
            if self.configs.z_dim:
                self.volume_model.train()
            else:
                base_pose = base_pose.detach().cpu().numpy()
        else:
            base_pose = None

        # reconstruct circle of pixels instead of whole image
        L_model = self.lattice.D // 2
        if self.configs.l_ramp_epochs > 0:
            Lramp = self.configs.l_start + int(
                self.current_epoch
                / self.configs.l_ramp_epochs
                * (self.configs.l_end - self.configs.l_start)
            )
            self.pose_search.Lmin = min(Lramp, self.configs.l_start)
            self.pose_search.Lmax = min(Lramp, self.configs.l_end)

            if (
                self.current_epoch < self.configs.l_ramp_epochs
                and self.configs.l_ramp_model
            ):
                L_model = self.pose_search.Lmax
        mask = self.lattice.get_circular_mask(L_model)

        def gen_slice(R):
            if z is None:
                slice_ = self.volume_model(self.lattice.coords[mask] @ R).view(B, -1)
            else:
                slice_ = self.volume_model(self.lattice.coords[mask] @ R, z).view(B, -1)

            if ctf_i is not None:
                slice_ *= ctf_i.view(B, -1)[:, mask]
            return slice_

        def translate(img):
            img = self.lattice.translate_ht(img, trans.unsqueeze(1), mask)
            return img.view(B, -1)

        y = y.view(B, -1)[:, mask]
        if use_tilt:
            tilt_ind = tilt_ind.view(B, -1)[:, mask]
        if trans is not None:
            y = translate(y)
        if use_tilt:
            tilt_ind = translate(tilt_ind)

        if use_tilt:
            rot_loss = F.mse_loss(gen_slice(rot), y)
            tilt_loss = F.mse_loss(
                gen_slice(bnb.tilt @ rot), tilt_ind  # type: ignore  # noqa: F821
            )
            losses["gen"] = (rot_loss + tilt_loss) / 2
        else:
            losses["gen"] = F.mse_loss(gen_slice(rot), y)

        # latent loss
        if self.configs.z_dim:
            losses["kld"] = torch.mean(
                -0.5 * torch.sum(1 + z_logvar - z_mu.pow(2) - z_logvar.exp(), dim=1),
                dim=0,
            )

            if torch.isnan(losses["kld"]):
                self.logger.info(z_mu[0])
                self.logger.info(z_logvar[0])
                raise RuntimeError("KLD is nan")

            if self.configs.beta_control is None:
                losses["total"] = self.beta * losses["kld"] / mask.sum().float()
            else:
                losses["total"] = (self.beta - losses["kld"]) ** 2 / mask.sum().float()
                losses["total"] *= self.configs.beta_control

            losses["total"] += losses["gen"]
            if lamb is not None and "eq_loss" in losses:
                losses["total"] += lamb * losses["eq_loss"]
        else:
            losses["total"] = losses["gen"]

        losses["total"].backward()
        self.volume_optimizer.step()
        rot = rot.detach().cpu().numpy()
        if trans is not None:
            trans = trans.detach().cpu().numpy()

        if self.pose_optimizer is not None:
            if self.current_epoch >= self.configs.pretrain:
                self.pose_optimizer.step()

        ind_tilt = tilt_ind or ind
        self.predicted_rots[ind_tilt] = rot.reshape(-1, 3, 3)
        if trans is not None:
            self.predicted_trans[ind_tilt] = trans.reshape(-1, 2)
        if base_pose is not None:
            self.base_poses.append((ind_np, base_pose))

        # logging
        for loss_k, loss_val in losses.items():
            if loss_k in self.accum_losses:
                self.accum_losses[loss_k] += loss_val.item() * len(ind)
            else:
                self.accum_losses[loss_k] = loss_val.item() * len(ind)

    def make_batch_summary(self):
        eq_log = ""
        if self.equivariance_lambda is not None:
            eq_log = (
                f"equivariance={self.accum_losses['eq_loss']:.4f}, "
                f"lambda={self.equivariance_lambda:.4f}, "
            )

        self.logger.info(
            f"# [Train Epoch: {self.current_epoch}/{self.configs.num_epochs}] "
            f"[{self.epoch_images_seen}/{self.image_count} images] "
            f"gen loss={self.accum_losses['gen']:.4f}, "
            f"kld={self.accum_losses['kld']:.4f}, beta={self.beta:.4f}, "
            f"{eq_log}loss={self.accum_losses['total']:.4f}"
        )

    def preprocess_input(self, y, trans):
        """Center the image."""
        return self.lattice.translate_ht(
            y.view(y.size(0), -1), trans.unsqueeze(1)
        ).view(y.size(0), self.resolution, self.resolution)

    def pretrain_batch(self, batch: dict[str, torch.Tensor]) -> None:
        """Pretrain the decoder using random initial poses."""
        y = batch["y"]
        use_tilt = "tilt_indices" in batch
        if use_tilt:
            batch["tilt_indices"].to(self.device)
        y = y.to(self.device)
        B = y.size(0)

        self.volume_model.train()
        self.volume_optimizer.zero_grad()

        # reconstruct circle of pixels instead of whole image
        mask = self.lattice.get_circular_mask(self.lattice.D // 2)
        rot = lie_tools.random_SO3(B, device=self.device)

        if self.configs.z_dim > 0:
            z = torch.randn((B, self.configs.z_dim), device=self.device)

            def gen_slice(R):
                _model = unparallelize(self.volume_model)
                assert isinstance(_model, HetOnlyVAE)
                return _model.decode(self.lattice.coords[mask] @ R, z).view(B, -1)

        else:

            def gen_slice(R):
                slice_ = self.volume_model(self.lattice.coords[mask] @ R)
                return slice_.view(B, -1)

        y = y.view(B, -1)[:, mask]
        if use_tilt:
            yt = batch["tilt_indices"].view(B, -1)[:, mask]
            loss = 0.5 * F.mse_loss(gen_slice(rot), y) + 0.5 * F.mse_loss(
                gen_slice(self.configs.tilt @ rot), yt
            )
        else:
            loss = F.mse_loss(gen_slice(rot), y)

        loss.backward()
        self.volume_optimizer.step()
        self.accum_losses["total"] = loss.item()

    def make_epoch_summary(self):
        """Save model weights, latent encoding z, and decoder volumes"""
        out_mrc = os.path.join(self.outdir, f"reconstruct.{self.epoch_lbl}.mrc")
        out_weights = os.path.join(self.outdir, f"weights.{self.epoch_lbl}.pkl")
        out_poses = os.path.join(self.outdir, f"pose.{self.epoch_lbl}.pkl")
        out_conf = os.path.join(self.outdir, f"z.{self.epoch_lbl}.pkl")

        # save model weights
        torch.save(
            {
                "epoch": self.current_epoch,
                "model_state_dict": unparallelize(self.volume_model).state_dict(),
                "optimizer_state_dict": self.volume_optimizer.state_dict(),
                "search_pose": (self.predicted_rots, self.predicted_trans),
            },
            out_weights,
        )

        if self.configs.z_dim > 0:
            self.volume_model.eval()

            with torch.no_grad():
                assert not self.volume_model.training
                z_mu_all = []
                z_logvar_all = []
                data_generator = dataset.make_dataloader(
                    self.data,
                    batch_size=self.configs.batch_size,
                    shuffler_size=self.configs.shuffler_size,
                    shuffle=False,
                )

                for minibatch in data_generator:
                    ind = minibatch["indices"]
                    y = minibatch["y"].to(self.device)
                    yt = None
                    if self.configs.tilt:
                        yt = minibatch["tilt_indices"].to(self.device)
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

                    _model = unparallelize(self.volume_model)
                    assert isinstance(_model, HetOnlyVAE)
                    z_mu, z_logvar = _model.encode(*input_)
                    z_mu_all.append(z_mu.detach().cpu().numpy())
                    z_logvar_all.append(z_logvar.detach().cpu().numpy())

            # save z
            z_mu_all, z_logvar_all = np.vstack(z_mu_all), np.vstack(z_logvar_all)
            with open(out_conf, "wb") as f:
                pickle.dump(z_mu_all, f)
                pickle.dump(z_logvar_all, f)

        with open(out_poses, "wb") as f:
            pickle.dump(
                (self.predicted_rots, self.predicted_trans / self.model_resolution), f
            )

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

    @property
    def in_pose_search_step(self) -> bool:
        is_in_pose_search = False
        if self.pose_search:
            if (self.current_epoch - 1) % self.configs.ps_freq == 0:
                is_in_pose_search = True

        return is_in_pose_search


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

    cryodrgn.config.save(config, out_config)
