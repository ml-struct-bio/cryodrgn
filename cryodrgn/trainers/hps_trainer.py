import os
import pickle
from datetime import datetime as dt
from collections import OrderedDict
from typing import Any, Callable

import numpy as np
import torch
from torch import nn
from torch.nn.parallel import DataParallel
import torch.nn.functional as F

import cryodrgn.config
from cryodrgn import ctf, dataset, lie_tools
from cryodrgn.losses import EquivarianceLoss
from cryodrgn.models.variational_autoencoder import unparallelize, HetOnlyVAE
from cryodrgn.models.neural_nets import get_decoder
from cryodrgn.pose_search import PoseSearch
from cryodrgn.trainers._base import ModelTrainer, ModelConfigurations


class HierarchicalPoseSearchConfigurations(ModelConfigurations):

    __slots__ = (
        "enc_only",
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
        "grid_niter",
        "ps_freq",
        "n_kept_poses",
        "pose_model_update_freq",
        "enc_layers",
        "enc_dim",
        "encode_mode",
        "enc_mask",
        "use_real",
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
            "data_norm": None,
            "l_ramp_epochs": 0,
            "l_ramp_model": 0,
            "reset_model_every": None,
            "reset_optim_every": None,
            "reset_optim_after_pretrain": None,
            "grid_niter": 4,
            "ps_freq": 5,
            "n_kept_poses": 8,
            "pose_model_update_freq": None,
            "enc_layers": 3,
            "enc_dim": 1024,
            "encode_mode": "resid",
            "enc_mask": None,
            "use_real": False,
            "dec_layers": 3,
            "dec_dim": 1024,
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

        if self.equivariance is not None:
            if not self.z_dim:
                raise ValueError(
                    "Cannot use equivariance with homogeneous reconstruction!."
                )
            if self.equivariance <= 0:
                raise ValueError("Regularization weight must be positive")

        if self.volume_domain is None:
            self.volume_domain = "fourier"


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
                3,
                self.resolution,
                self.configs.dec_layers,
                self.configs.dec_dim,
                self.configs.volume_domain,
                self.configs.pe_type,
                self.configs.pe_dim,
                activation=self.activation,
                feat_sigma=self.configs.feat_sigma,
            )
        else:
            enc_mask, in_dim = self.mask_dimensions

            model = HetOnlyVAE(
                self.lattice,
                self.configs.enc_layers,
                self.configs.enc_dim,
                self.configs.dec_layers,
                self.configs.dec_dim,
                in_dim,
                self.configs.z_dim,
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

    @property
    def global_it(self) -> int:
        return self.particle_count * self.current_epoch + self.batch_it

    def train_epoch(self):
        self.configs: HierarchicalPoseSearchConfigurations
        te = dt.now()

        accum_losses = {"total": 0, "gen": 0, "kld": 0, "eq": 0}
        all_poses = list()
        base_poses = list()
        self.batch_it = 0

        # reset the model if asked for
        if self.configs.reset_model_every:
            if self.current_epoch % self.configs.reset_model_every == 0:
                self.logger.info(">> Resetting model")
                self.volume_model = self.make_volume_model()

        # reset the optimizer if asked for
        if self.configs.reset_optim_every:
            if self.current_epoch % self.configs.reset_optim_every == 0:
                self.logger.info(">> Resetting optim")
                self.volume_optimizer = self.volume_optim_class(
                    self.volume_model.parameters(),
                    lr=self.configs.learning_rate,
                    weight_decay=self.configs.weight_decay,
                )

        for batch, tilt_ind, ind in self.data_iterator:
            ind_np = ind.cpu().numpy()
            y = batch.to(self.device)
            self.batch_it += len(batch[0])

            if tilt_ind is not None:
                tilt_ind = tilt_ind.to(self.device)
                assert all(tilt_ind >= 0), tilt_ind

            lamb = None
            if (
                self.equivariance_lambda is not None
                and self.equivariance_loss is not None
            ):
                lamb = self.equivariance_lambda(self.global_it)
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
            batch_poses = None
            if self.pose_search and self.current_epoch % self.configs.ps_freq != 0:
                self.logger.info("Using previous iteration poses")
                batch_poses = [
                    torch.tensor(x[ind_np], device=self.device)
                    for x in self.sorted_poses
                ]  # type: ignore

            self.conf_search_particles += len(batch[0])
            if self.configs.pose_model_update_freq:
                if (
                    self.current_epoch == 0
                    or self.conf_search_particles > self.configs.pose_model_update_freq
                ):
                    self.pose_model.load_state_dict(self.volume_model.state_dict())
                    self.conf_search_particles = 0

            dose_filters = None
            if self.configs.tilt:
                if self.pose_tracker and not self.pose_search:
                    batch_poses = self.pose_tracker.get_pose(tilt_ind.view(-1))

                ctf_param = (
                    self.ctf_params[tilt_ind.view(-1)]
                    if self.ctf_params is not None
                    else None
                )
                y = y.view(-1, self.resolution, self.resolution, self.resolution)
                Apix = self.ctf_params[0, 0] if self.ctf_params is not None else None

                if self.configs.dose_per_tilt is not None:
                    dose_filters = self.data.get_dose_filters(
                        tilt_ind, self.lattice, Apix
                    )

            else:
                if self.pose_tracker and not self.pose_search:
                    batch_poses = self.pose_tracker.get_pose(ind)
                ctf_param = (
                    self.ctf_params[ind] if self.ctf_params is not None else None
                )

            # train the model
            losses, new_poses, new_base_poses = self.train_step(
                batch,
                tilt_ind,
                equivariance_tuple,
                rot=batch_poses[0],
                trans=batch_poses[1],
                ctf_params=ctf_param,
            )

            if self.pose_optimizer is not None:
                if self.current_epoch >= self.configs.pretrain:
                    self.pose_optimizer.step()

            all_poses.append((ind.cpu().numpy(), new_poses))
            if new_base_poses:
                base_poses.append((ind_np, new_base_poses))

            # logging
            for k in losses:
                accum_losses[k] += losses[k] * len(ind)

        eq_log = (
            "equivariance = {:.4f}, ".format(accum_losses["eq"] / self.image_count)
            if self.configs.equivariance
            else ""
        )
        avg_gen_loss = accum_losses["gen"] / self.image_count
        kld_loss = accum_losses["kld"] / self.image_count
        total_loss = accum_losses["total"] / self.image_count
        self.logger.info(
            f"# =====> Epoch: {self.current_epoch} "
            f"Average gen loss = {avg_gen_loss:.4}, KLD = {kld_loss:.4f},"
            f" {eq_log}total loss = {total_loss:.4f}; "
            f"Finished in {dt.now() - te,}"
        )

        # sort poses if necessary
        if self.pose_search:
            if all_poses:
                ind = [x[0] for x in all_poses]
                ind = np.concatenate(ind)
                rot = [x[1][0] for x in all_poses]
                rot = np.concatenate(rot)
                rot = rot[np.argsort(ind)]

                if len(all_poses[0][1]) == 2:
                    trans = [x[1][1] for x in all_poses]
                    trans = np.concatenate(trans)
                    trans = trans[np.argsort(ind)]
                    self.sorted_poses = (rot, trans)
                else:
                    self.sorted_poses = (rot,)
            else:
                self.sorted_poses = tuple()

    def train_step(
        self,
        y,
        tilt_ind,
        equivariance=None,
        rot=None,
        trans=None,
        ctf_params=None,
    ) -> tuple[dict[str, float], tuple[torch.Tensor, torch.Tensor], tuple]:
        use_tilt = tilt_ind is not None
        use_ctf = ctf_params is not None
        B = y.size(0)

        if self.beta_schedule is not None:
            beta = self.beta_schedule(self.global_it)
        else:
            beta = None

        ctf_i = None
        if use_ctf:
            freqs = self.lattice.freqs2d.unsqueeze(0).expand(
                B, *self.lattice.freqs2d.shape
            ) / ctf_params[:, 0].view(B, 1, 1)
            ctf_i = ctf.compute_ctf(freqs, *torch.split(ctf_params[:, 1:], 1, 1)).view(
                B, self.resolution, self.resolution
            )

        if trans is not None:
            y = self.preprocess_input(y, trans)

        # VAE inference of z
        self.volume_optimizer.zero_grad()
        self.volume_model.train()
        input_ = (y, tilt_ind) if use_tilt else (y,)
        if ctf_i is not None:
            input_ = (x * ctf_i.sign() for x in input_)  # phase flip by the ctf

        _model = unparallelize(self.volume_model)
        assert isinstance(_model, HetOnlyVAE)
        z_mu, z_logvar = _model.encode(*input_)
        z = _model.reparameterize(z_mu, z_logvar)

        lamb = None
        losses = dict()
        if equivariance is not None:
            lamb, equivariance_loss = equivariance
            losses["eq_loss"] = equivariance_loss(y, z_mu)

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
            losses = {"gen": (rot_loss + tilt_loss) / 2}
        else:
            losses = {"gen": F.mse_loss(gen_slice(rot), y)}

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
                losses["total"] = beta * losses["kld"] / mask.sum().float()
            else:
                losses["total"] = (beta - losses["kld"]) ** 2 / mask.sum().float()
                losses["total"] *= self.configs.beta_control

            losses["total"] += losses["gen"]
            if lamb is not None and "eq_loss" in losses:
                losses["total"] += lamb * losses["eq_loss"]

        losses["total"].backward()
        self.volume_optimizer.step()
        rot = rot.detach().cpu().numpy()
        if trans is not None:
            trans = trans.detach().cpu().numpy()

        if self.batch_it % self.configs.log_interval == 0:
            eq_log = (
                f"equivariance={losses['eq_loss']:.4f}, lambda={lamb:.4f}, "
                if "eq_loss" in losses and lamb is not None
                else ""
            )
            self.logger.info(
                f"# [Train Epoch: {self.current_epoch}/{self.configs.num_epochs}] "
                f"[{self.batch_it}/{self.image_count} images] "
                f"gen loss={losses['gen']:.4f}, "
                f"kld={losses['kld']:.4f}, beta={beta:.4f}, "
                f"{eq_log}loss={losses['total']:.4f}"
            )

        return {k: loss.item() for k, loss in losses.items()}, (rot, trans), base_pose

    def preprocess_input(self, y, trans):
        """Center the image."""
        return self.lattice.translate_ht(
            y.view(y.size(0), -1), trans.unsqueeze(1)
        ).view(y.size(0), self.resolution, self.resolution)

    def pretrain_step(self, batch):
        if self.z_dim > 0:
            y, yt = batch
            use_tilt = yt is not None
            B = y.size(0)

            self.model.train()
            self.optim.zero_grad()

            rot = lie_tools.random_SO3(B, device=y.device)
            z = torch.randn((B, self.z_dim), device=y.device)

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

    def make_summary(
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

        if self.z_dim > 0:
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
