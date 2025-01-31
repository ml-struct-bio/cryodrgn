"""The training engine for cryoDRGN reconstruction models up to and including v3.

This `HierarchicalPoseSearchTrainer` class (and its associated configuration class
`HierarchicalPoseSearchConfigurations`) are used to replicate the behaviour of the
cryoDRGN reconstruction algorithm as featured in commands such as `cryodrgn train_vae`,
`cryodrgn abinit_homo`, etc. in versions of the package before v4.0.0 — that is, before
the merge of the cryoDRGN codebase with the offshoot DRGN-AI amortized inference
reconstruction algorithm.

"""
import os
import pickle
from dataclasses import dataclass
from typing import Any, Callable
import numpy as np
import torch
from torch import nn
import torch.nn.functional as F

from cryodrgn import ctf, dataset
from cryodrgn.models import lie_tools
from cryodrgn.models.losses import EquivarianceLoss
from cryodrgn.models.variational_autoencoder import unparallelize, HetOnlyVAE
from cryodrgn.models.neural_nets import get_decoder
from cryodrgn.models.pose_search import PoseSearch
from cryodrgn.trainers import (
    ReconstructionModelTrainer,
    ReconstructionModelConfigurations,
)
from cryodrgn.mrcfile import write_mrc
from cryodrgn.pose import PoseTracker

try:
    import apex.amp as amp  # type: ignore  # PYR01
except ImportError:
    pass


@dataclass
class HierarchicalPoseSearchConfigurations(ReconstructionModelConfigurations):

    # a parameter belongs to this configuration set if and only if it has a default
    # value defined here, note that children classes inherit these from parents
    model = "hps"

    # specifying size and type of model encoder and decoder
    enc_layers: int = None
    enc_dim: int = None
    encode_mode: str = "resid"
    enc_mask: bool = None
    tilt_enc_only: bool = False
    dec_layers: int = None
    dec_dim: int = None
    # how often pose search is done and other pose search parameters
    ps_freq: int = 5
    n_kept_poses: int = 8
    pose_model_update_freq: int = None
    grid_niter: int = 4
    # other learning model parameters
    beta: float = None
    beta_control: float = None
    equivariance: bool = None
    equivariance_start: int = 100000
    equivariance_stop: int = 200000
    l_ramp_epochs: int = None
    l_ramp_model: int = 0
    # resetting every certain number of epochs
    reset_model_every: int = None
    reset_optim_every: int = None

    def __post_init__(self) -> None:
        super().__post_init__()
        assert self.model == "hps"

        if self.capture_setup is not None:
            if self.capture_setup == "et":
                self.subtomo_averaging = True
                self.shuffler_size = 0
                self.num_workers = 0
                self.t_extent = 0.0

        if self.dataset is None:
            if self.particles is None:
                raise ValueError(
                    "As dataset was not specified, please specify particles!"
                )

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
            if self.pose_estimation == "fixed":
                self.volume_domain = "fourier"
            else:
                self.volume_domain = "hartley"

        if self.l_ramp_epochs is None:
            self.l_ramp_epochs = 25 if self.z_dim == 0 else 0


class HierarchicalPoseSearchTrainer(ReconstructionModelTrainer):

    configs: HierarchicalPoseSearchConfigurations
    config_cls = HierarchicalPoseSearchConfigurations
    label = "cDRGN v3 training"

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

    def make_reconstruction_model(self) -> nn.Module:
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

        # set beta schedule
        if self.configs.z_dim:
            beta = self.configs.beta or self.configs.z_dim**-1

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
        # also check e.g. enc_mask dim?
        if self.configs.amp:
            if self.configs.enc_dim % 8 != 0:
                self.logger.warning(
                    f"Encoder hidden layer dimension {self.configs.enc_dim} not "
                    f"divisible by 8 and thus not optimal for AMP training!"
                )
            if self.configs.dec_dim % 8 != 0:
                self.logger.warning(
                    f"Decoder hidden layer dimension {self.configs.dec_dim} not "
                    f"divisible by 8 and thus not optimal for AMP training!"
                )
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
                    self.reconstruction_model, self.resolution
                )

        self.pose_tracker = None
        self.pose_search = None
        self.pose_model = None
        self.do_pretrain = self.configs.pose_estimation == "abinit"

        if self.configs.poses:
            self.pose_tracker = PoseTracker.load(
                infile=self.configs.poses,
                Nimg=self.image_count,
                D=self.resolution,
                emb_type="s2s2" if self.configs.pose_estimation == "refine" else None,
                ind=self.ind,
                device=self.device,
            )

        if self.configs.pose_estimation == "refine":
            self.pose_optimizer = torch.optim.SparseAdam(
                list(self.pose_tracker.parameters()), lr=self.configs.pose_learning_rate
            )
        elif self.configs.pose_estimation == "abinit":
            if self.configs.pose_model_update_freq:
                assert not self.configs.multigpu, "TODO"
                self.pose_model = self.make_reconstruction_model()
                self.pose_model.to(self.device)
                self.pose_model.eval()
            else:
                self.pose_model = self.reconstruction_model

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

    def begin_epoch(self) -> None:
        if self.configs.l_ramp_epochs > 0 and self.pose_search is not None:
            Lramp = self.configs.l_start + int(
                self.current_epoch
                / self.configs.l_ramp_epochs
                * (self.configs.l_end - self.configs.l_start)
            )
            self.pose_search.Lmin = min(Lramp, self.configs.l_start)
            self.pose_search.Lmax = min(Lramp, self.configs.l_end)

    def train_batch(self, batch: dict[str, torch.Tensor]) -> tuple:
        y = batch["y"]
        ind = batch["indices"]

        y = y.to(self.device)
        ind_np = ind.cpu().numpy()
        B = y.size(0)

        if self.configs.tilt:
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

            if self.configs.dose_per_tilt is not None:
                dose_filters = self.data.get_dose_filters(
                    tilt_ind, self.lattice, self.apix
                )

        else:
            ctf_param = self.ctf_params[ind] if self.ctf_params is not None else None

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

        if not self.pose_search:
            rot, trans = self.pose_tracker.get_pose(ind)
        else:
            rot, trans = None, None

        lamb = None
        losses = dict()
        if self.configs.z_dim > 0:
            if trans is not None:
                y_trans = self.lattice.translate_ht(
                    y.view(B, -1), trans.unsqueeze(1)
                ).view(B, self.resolution, self.resolution)
            else:
                y_trans = y.clone()

            if yr is not None:
                input_ = (yr, tilt_ind) if tilt_ind is not None else (yr,)
            else:
                input_ = (y_trans, tilt_ind) if tilt_ind is not None else (y_trans,)

            if ctf_i is not None:
                input_ = (x * ctf_i.sign() for x in input_)  # phase flip by the ctf

            _model = unparallelize(self.reconstruction_model)
            assert isinstance(_model, HetOnlyVAE)
            z_mu, z_logvar = _model.encode(*input_)
            z = _model.reparameterize(z_mu, z_logvar)

            if equivariance_tuple is not None:
                lamb, equivariance_loss = equivariance_tuple
                losses["eq"] = equivariance_loss(y, z_mu)
        else:
            z_mu = z_logvar = z = None

        # getting the poses; execute a pose search if poses not given
        if self.pose_search:
            if self.in_pose_search_step:
                self.base_pose = None
                self.reconstruction_model.eval()

                with torch.no_grad():
                    rot, trans, self.base_pose = self.pose_search.opt_theta_trans(
                        y,
                        z=z,
                        images_tilt=None if self.configs.tilt_enc_only else tilt_ind,
                        init_poses=self.base_pose,
                        ctf_i=ctf_i,
                    )
                if self.configs.z_dim > 0:
                    self.reconstruction_model.train()

            else:
                if self.epoch_batch_count == 1:
                    self.logger.info("Using previous iteration's learned poses...")

                if isinstance(self.base_pose, torch.Tensor):
                    self.base_pose = self.base_pose.detach().cpu().numpy()
                rot = torch.tensor(
                    self.predicted_rots[ind_np].astype(np.float32), device=self.device
                )
                if not self.configs.no_trans:
                    trans = torch.tensor(
                        self.predicted_trans[ind_np].astype(np.float32),
                        device=self.device,
                    )
                else:
                    trans = None

        self.conf_search_particles += B
        if self.configs.pose_model_update_freq:
            if (
                self.current_epoch == 1
                or self.conf_search_particles > self.configs.pose_model_update_freq
            ):
                self.pose_model.load_state_dict(self.reconstruction_model.state_dict())
                self.conf_search_particles = 0

        # reconstruct circle of pixels instead of whole image
        L_model = self.lattice.D // 2
        if (
            self.current_epoch < self.configs.l_ramp_epochs
            and self.configs.l_ramp_model
        ):
            L_model = self.pose_search.Lmax
        mask = self.lattice.get_circular_mask(L_model)

        def gen_slice(R):
            lat_coords = self.lattice.coords[mask] / self.lattice.extent / 2

            if z is None:
                slice_ = self.reconstruction_model(lat_coords @ R).view(B, -1)
            else:
                slice_ = self.reconstruction_model(lat_coords @ R, z).view(B, -1)

            if ctf_i is not None:
                slice_ *= ctf_i.view(B, -1)[:, mask]

            return slice_

        def translate(img):
            img_trans = img.view(B, -1)[:, mask]

            if trans is not None:
                img_trans = self.lattice.translate_ht(
                    img_trans, trans.unsqueeze(1), mask
                )

            return img_trans

        if self.configs.tilt:
            tilt_ind = tilt_ind.view(B, -1)[:, mask]
            tilt_ind = translate(tilt_ind)

            if dose_filters:
                rot_loss = F.mse_loss(gen_slice(rot), y)
            else:
                y_recon = torch.mul(gen_slice(rot), dose_filters[:, mask])
                rot_loss = F.mse_loss(y_recon, y)

            tilt_loss = F.mse_loss(
                gen_slice(bnb.tilt @ rot), tilt_ind  # type: ignore  # noqa: F821
            )
            losses["gen"] = (rot_loss + tilt_loss) / 2

        else:
            losses["gen"] = F.mse_loss(gen_slice(rot), translate(y).view(B, -1))

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
            if lamb is not None and "eq" in losses:
                losses["total"] += lamb * losses["eq"]
        else:
            losses["total"] = losses["gen"]

        if self.configs.amp:
            if self.scaler is not None:  # torch mixed precision
                self.scaler.scale(losses["total"]).backward()
                self.scaler.step(self.reconstruction_optimizer)
                self.scaler.update()
            else:  # apex.amp mixed precision
                with amp.scale_loss(
                    losses["total"], self.reconstruction_optimizer
                ) as scaled_loss:
                    scaled_loss.backward()
                self.reconstruction_optimizer.step()
        else:
            losses["total"].backward()
            self.reconstruction_optimizer.step()

        return losses, tilt_ind, ind, rot, trans, z_mu, z_logvar

    def get_configs(self) -> dict[str, Any]:
        """Retrieves all given and inferred configurations for downstream use."""
        configs = super().get_configs()
        enc_mask = self.configs.enc_mask or self.resolution // 2

        configs["model_args"].update(
            dict(
                qlayers=self.configs.enc_layers,
                qdim=self.configs.enc_dim,
                players=self.configs.dec_layers,
                pdim=self.configs.dec_dim,
                encode_mode=self.configs.encode_mode,
                enc_mask=enc_mask,
            )
        )

        return configs

    def print_batch_summary(self, losses: dict[str, float]) -> None:
        """Create a summary at the end of a training batch and print it to the log."""
        eq_log = ""
        if self.equivariance_lambda is not None and "eq" in losses:
            eq_log = (
                f"equivariance={losses['eq']:.4f}, "
                f"lambda={self.equivariance_lambda:.4f}, "
            )

        batch_str = (
            f"### Epoch [{self.current_epoch}/{self.configs.num_epochs}], "
            f"Batch [{self.epoch_batch_count}] "
            f"({self.epoch_images_seen}/{self.image_count} images); "
            f"gen loss={losses['gen']:.6g}, "
        )
        if "kld" in losses:
            batch_str += f"kld={losses['kld']:.4f}, "
        if self.beta is not None:
            batch_str += f"beta={self.beta:.4f}, "

        self.logger.info(batch_str + f"{eq_log}")

    def pretrain_batch(self, batch: dict[str, torch.Tensor]) -> None:
        """Pretrain the decoder using random initial poses."""
        y = batch["y"]
        if self.configs.tilt:
            batch["tilt_indices"].to(self.device)
        y = y.to(self.device)
        B = y.size(0)

        self.reconstruction_model.train()
        self.reconstruction_optimizer.zero_grad()

        # reconstruct circle of pixels instead of whole image
        mask = self.lattice.get_circular_mask(self.lattice.D // 2)
        rot = lie_tools.random_SO3(B, device=self.device)

        if self.configs.z_dim > 0:
            z = torch.randn((B, self.configs.z_dim), device=self.device)

            def gen_slice(R):
                _model = unparallelize(self.reconstruction_model)
                assert isinstance(_model, HetOnlyVAE)
                return _model.decode(self.lattice.coords[mask] @ R, z).view(B, -1)

        else:

            def gen_slice(R):
                slice_ = self.reconstruction_model(self.lattice.coords[mask] @ R)
                return slice_.view(B, -1)

        y = y.view(B, -1)[:, mask]
        if self.configs.tilt:
            yt = batch["tilt_indices"].view(B, -1)[:, mask]
            loss = 0.5 * F.mse_loss(gen_slice(rot), y) + 0.5 * F.mse_loss(
                gen_slice(self.configs.tilt @ rot), yt
            )
        else:
            loss = F.mse_loss(gen_slice(rot), y)

        loss.backward()
        self.reconstruction_optimizer.step()
        self.accum_losses["total"] += loss.item() * B

    def save_epoch_data(self):
        """Save model weights, latent encoding z, and decoder volumes"""
        out_weights = os.path.join(self.outdir, f"weights.{self.epoch_lbl}.pkl")
        out_poses = os.path.join(self.outdir, f"pose.{self.epoch_lbl}.pkl")
        out_conf = os.path.join(self.outdir, f"conf.{self.epoch_lbl}.pkl")

        # save a reconstructed volume by evaluating model on a 3d lattice
        if self.configs.z_dim == 0:
            out_mrc = os.path.join(self.outdir, f"reconstruct.{self.epoch_lbl}.mrc")
            self.reconstruction_model.eval()
            vol = self.reconstruction_model.eval_volume(
                coords=self.lattice.coords,
                D=self.lattice.D,
                extent=self.lattice.extent,
                norm=self.data.norm,
            )
            write_mrc(out_mrc, vol, Apix=self.apix)

        # save model weights
        torch.save(
            {
                "epoch": self.current_epoch,
                "model_state_dict": unparallelize(
                    self.reconstruction_model
                ).state_dict(),
                "optimizer_state_dict": self.reconstruction_optimizer.state_dict(),
                "search_pose": (self.predicted_rots, self.predicted_trans),
            },
            out_weights,
        )

        # If we are doing heterogeneous reconstruction, also save latent conformations
        if self.configs.z_dim > 0:
            self.reconstruction_model.eval()

            with torch.no_grad():
                assert not self.reconstruction_model.training
                z_mu_all = []
                z_logvar_all = []
                data_generator = dataset.make_dataloader(
                    self.data,
                    batch_size=self.configs.batch_size,
                    shuffler_size=self.configs.shuffler_size,
                    shuffle=False,
                    seed=self.configs.seed,
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

                    if self.pose_search is None and self.pose_tracker is not None:
                        if self.pose_tracker.trans is not None:
                            y = self.lattice.translate_ht(
                                y.view(B, -1), self.pose_tracker.trans[ind].unsqueeze(1)
                            ).view(B, D, D)

                    input_ = (y, yt) if yt is not None else (y,)
                    if c is not None:
                        input_ = (x * c.sign() for x in input_)  # phase flip by the ctf

                    _model = unparallelize(self.reconstruction_model)
                    assert isinstance(_model, HetOnlyVAE)
                    z_mu, z_logvar = _model.encode(*input_)
                    z_mu_all.append(z_mu.detach().cpu().numpy())
                    z_logvar_all.append(z_logvar.detach().cpu().numpy())

            # save z
            z_mu_all, z_logvar_all = np.vstack(z_mu_all), np.vstack(z_logvar_all)
            with open(out_conf, "wb") as f:
                pickle.dump(z_mu_all, f)
                pickle.dump(z_logvar_all, f)

        if self.configs.pose_estimation == "refine":
            self.pose_tracker.save(out_poses)
        elif self.configs.pose_estimation == "abinit":
            with open(out_poses, "wb") as f:
                pickle.dump(
                    (self.predicted_rots, self.predicted_trans / self.model_resolution),
                    f,
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
        """Whether the current epoch of training will include searching over poses."""
        is_in_pose_search = False
        if self.pose_search:
            if (self.current_epoch - 1) % self.configs.ps_freq == 0:
                is_in_pose_search = True

        return is_in_pose_search
