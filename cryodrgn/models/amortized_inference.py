"""Pytorch models implementing amortized inference."""

import numpy as np
import torch
import torch.nn as nn
import time
from cryodrgn.models.neural_nets import half_linear, single_linear
from cryodrgn import fft, masking
from cryodrgn.models import lie_tools, pose_search_amortized


class DRGNai(nn.Module):
    def __init__(
        self,
        lattice,
        output_mask,
        n_particles_dataset,
        n_tilts_dataset,
        cnn_params,
        conf_regressor_params,
        hyper_volume_params,
        resolution_encoder=64,
        no_trans=False,
        use_gt_poses=False,
        use_gt_trans=False,
        will_use_point_estimates=False,
        ps_params=None,
        verbose_time=False,
        pretrain_with_gt_poses=False,
        n_tilts_pose_search=1,
    ):
        """
        lattice: Lattice
        output_mask: Mask
        n_particles_dataset: int
        n_tilts_dataset: int
        cnn_params: dict
            conf: bool
            depth_cnn: int
            channels_cnn: int
            kernel_size_cnn: int
        conf_regressor_params: dict
            z_dim: int
            std_z_init: float
            variational: bool
        hyper_volume_params: dict
            n_layers: int
            hidden_dim: it
            pe_type: str
            pe_dim: int
            feat_sigma: float
            domain: str
            extent: float
            pe_type_conf: str or None
        resolution_encoder: int
        no_trans: bool
        use_gt_poses: bool
        use_gt_trans: bool
        will_use_point_estimates: bool
        ps_params: dict
            'l_min': int
            'l_max': int
            't_extent': float
            't_n_grid': int
            'niter': int
            'nkeptposes': int
            'base_healpy': int
            't_xshift': float
            't_yshift': float
        verbose_time: bool
        pretrain_with_gt_poses: bool
        n_tilts_pose_search: int
        """
        super(DRGNai, self).__init__()
        self.lattice = lattice
        self.D = lattice.D
        self.output_mask = output_mask
        self.verbose_time = verbose_time
        self.n_tilts_pose_search = n_tilts_pose_search

        # will be moved to the local gpu of each replica
        self.coords = nn.Parameter(
            self.lattice.coords, requires_grad=False
        )  # [D * D, 3]
        self.freqs2d = nn.Parameter(self.lattice.freqs2d, requires_grad=False)
        if ps_params is not None:
            self.base_shifts = nn.Parameter(
                pose_search_amortized.get_base_shifts(ps_params), requires_grad=False
            )
            self.base_rot = nn.Parameter(
                pose_search_amortized.get_base_rot(ps_params), requires_grad=False
            )
            self.so3_base_quat = nn.Parameter(
                pose_search_amortized.get_so3_base_quat(ps_params), requires_grad=False
            )
            self.base_inplane = nn.Parameter(
                pose_search_amortized.get_base_inplane(ps_params), requires_grad=False
            )

        self.no_trans = no_trans
        self.z_dim = conf_regressor_params["z_dim"]
        self.variational_conf = conf_regressor_params["variational"]
        self.std_z_init = conf_regressor_params["std_z_init"]

        self.pose_only = False
        self.use_point_estimates = False
        self.pretrain = False
        self.is_in_pose_search_step = False
        self.use_point_estimates_conf = False

        # pose
        if not use_gt_poses and will_use_point_estimates:
            self.pose_table = PoseTable(
                n_tilts_dataset, self.no_trans, self.D, use_gt_trans
            )

        # conformation
        if self.z_dim > 0:
            if cnn_params["conf"]:
                self.conf_cnn = SharedCNN(
                    resolution_encoder
                    if resolution_encoder is not None
                    else self.D - 1,
                    cnn_params["depth_cnn"],
                    cnn_params["channels_cnn"],
                    cnn_params["kernel_size_cnn"],
                    1,
                )

                final_channels = self.conf_cnn.final_channels
                final_size = self.conf_cnn.final_size

                self.conf_regressor = ConfRegressor(
                    final_channels,
                    final_size,
                    conf_regressor_params["z_dim"],
                    conf_regressor_params["std_z_init"],
                    conf_regressor_params["variational"],
                )

            else:
                self.conf_table = ConfTable(
                    n_particles_dataset,
                    self.z_dim,
                    conf_regressor_params["variational"],
                    conf_regressor_params["std_z_init"],
                )

        self.use_gt_poses = use_gt_poses
        self.use_gt_trans = use_gt_trans
        self.pretrain_with_gt_poses = pretrain_with_gt_poses

        # pose search parameters
        self.ps_params = ps_params
        self.trans_search_factor = None
        if ps_params is not None and ps_params["no_trans_search_at_pose_search"]:
            self.trans_search_factor = 0.0

        # hyper-volume
        if not hyper_volume_params["explicit_volume"]:
            self.hypervolume = HyperVolume(
                self.D,
                self.z_dim,
                hyper_volume_params["n_layers"],
                hyper_volume_params["hidden_dim"],
                hyper_volume_params["pe_type"],
                hyper_volume_params["pe_dim"],
                hyper_volume_params["feat_sigma"],
                hyper_volume_params["domain"],
                pe_type_conf=hyper_volume_params["pe_type_conf"],
            )
        else:
            self.hypervolume = VolumeExplicit(
                self.D, hyper_volume_params["domain"], hyper_volume_params["extent"]
            )

    def update_trans_search_factor(self, ratio):
        if self.trans_search_factor is not None:
            self.trans_search_factor = ratio

    def forward(self, in_dict):
        """
        Merges encoding and decoding steps for data parallelization.

        in_dict: dict (not in pose supervision step)
            index: [batch_size]
            y: [batch_size(, n_tilts), D, D]
            y_real: [batch_size(, n_tilts), D - 1, D - 1]
            R: [batch_size(, n_tilts), 3, 3]
            t: [batch_size(, n_tilts), 2]
            tilt_index: [batch_size(, n_tilts)]
            ctf: [batch_size(, n_tilts), D, D]
        or (in pose supervision step)
            y_real: [batch_size, D - 1, D - 1]
            indices: [batch_size]
            R: [batch_size, 3, 3]
            t: [batch_size, 2]

        out_dict:
            y_pred: [batch_size, n_pts]
            y_gt_processed: [batch_size, n_pts]
            R: [batch_size, 3, 3]
            t: [batch_size, 2]
            z: [batch_size, z_dim]
            z_logvar: [batch_size, z_dim]
            time_encoder: [1]
            time_decoder: [1]
        """
        device = self.coords.device
        if self.verbose_time:
            torch.cuda.synchronize(device)
        start_time_encoder = time.time()
        latent_variables_dict = self.encode(in_dict, ctf=in_dict["ctf"])

        if in_dict["tilt_indices"] is not None:
            in_dict["tilt_indices"] = in_dict["tilt_indices"].reshape(-1)

        if self.verbose_time:
            torch.cuda.synchronize(device)
        start_time_decoder = time.time()
        y_pred, y_gt_processed, times, latent_variables_dict = self.decode(
            latent_variables_dict, in_dict["ctf"], in_dict["y"]
        )
        if self.verbose_time:
            torch.cuda.synchronize()
        end_time = time.time()
        out_dict = {"y_pred": y_pred, "y_gt_processed": y_gt_processed}
        if self.verbose_time:
            out_dict["time_encoder"] = (
                torch.tensor([start_time_decoder - start_time_encoder])
                .float()
                .to(device)
            )
            out_dict["time_decoder"] = (
                torch.tensor([end_time - start_time_decoder]).float().to(device)
            )
            out_dict["time_decoder_coords"] = (
                torch.tensor([times["coords"]]).float().to(device)
            )
            out_dict["time_decoder_query"] = (
                torch.tensor([times["query"]]).float().to(device)
            )
        for key in latent_variables_dict.keys():
            out_dict[key] = latent_variables_dict[key]
        return out_dict

    @staticmethod
    def process_y_real(in_dict):
        y_real = in_dict["y_real"]
        return y_real[..., None, :, :]

    def encode(self, in_dict, ctf=None):
        """
        in_dict: dict
            index: [batch_size]
            y: [batch_size(, n_tilts), D, D]
            y_real: [batch_size(, n_tilts), D - 1, D - 1]
            R: [batch_size(, n_tilts), 3, 3]
            t: [batch_size(, n_tilts), 2]
            tilt_index: [batch_size( * n_tilts)]
        ctf: [batch_size(, n_tilts), D, D]

        output: dict
            R: [batch_size(, n_tilts), 3, 3]
            t: [batch_size(, n_tilts), 2]
            z: [batch_size, z_dim]
            z_logvar: [batch_size, z_dim]
        """
        latent_variables_dict = {}
        z = None
        device = self.coords.device
        batch_size = in_dict["y"].shape[0]

        # conformation
        if self.z_dim > 0:
            # pretrain and pose only
            if self.pose_only:
                z = self.std_z_init * torch.randn(
                    (batch_size, self.z_dim), dtype=torch.float32, device=device
                )
                conf_dict = {"z": z}
                if self.variational_conf:
                    logvar = torch.ones(
                        (batch_size, self.z_dim), dtype=torch.float32, device=device
                    )
                    conf_dict["z_logvar"] = logvar
            # amortized inference
            elif not self.use_point_estimates_conf:
                y_real = self.process_y_real(in_dict)
                particles_real = y_real.mean(1) if y_real.ndim == 5 else y_real
                conf_features = self.conf_cnn(particles_real)
                conf_dict = self.conf_regressor(conf_features)
            # latent optimization
            else:
                conf_dict = self.conf_table(in_dict)
            z = conf_dict["z"]
            for key in conf_dict:
                latent_variables_dict[key] = conf_dict[key]

        # use gt poses
        if self.use_gt_poses or (self.pretrain and self.pretrain_with_gt_poses):
            rots = in_dict["R"]
            pose_dict = {"R": rots}
            if not self.no_trans:
                trans = in_dict["t"]
                pose_dict["t"] = trans

        # random poses
        # TODO: refactor to make pretraining a method of the trainer like for HPS?
        elif self.pretrain:
            in_dim = in_dict["y"].shape[:-2]
            device = in_dict["y"].device
            pose_dict = {"R": lie_tools.random_rotmat(np.prod(in_dim), device=device)}
            pose_dict["R"] = pose_dict["R"].reshape(*in_dim, 3, 3)
            if not self.no_trans:
                pose_dict["t"] = torch.zeros((*in_dim, 2)).float().to(device)

        # use pose search
        elif self.is_in_pose_search_step:
            self.hypervolume.eval()
            rot, trans, _ = pose_search_amortized.opt_theta_trans(
                self,
                in_dict["y"],
                self.lattice,
                self.ps_params,
                z=z,
                ctf_i=ctf,
                gt_trans=in_dict["t"]
                if not self.no_trans and self.use_gt_trans
                else None,
                trans_search_factor=self.trans_search_factor,
            )
            pose_dict = {"R": rot, "indices": in_dict["indices"]}
            if not self.no_trans:
                pose_dict["t"] = trans
            self.hypervolume.train()

        # use point estimates
        else:
            assert self.use_point_estimates
            pose_dict = self.pose_table(in_dict)

        for key in pose_dict:
            latent_variables_dict[key] = pose_dict[key]

        return latent_variables_dict

    def decode(self, latent_variables_dict, ctf_local, y_gt):
        """
        latent_variables_dict: dict
            R: [batch_size(, n_tilts), 3, 3]
            t: [batch_size(, n_tilts), 2]
            z: [batch_size, z_dim]
            z_logvar: [batch_size, z_dim]
        ctf_local: [batch_size(, n_tilts), D, D]
        y: [batch_size(, n_tilts), D, D]

        output: [batch_size(, n_tilts), n_pts], [batch_size, n_pts], dict ('coords': float, 'query': float)
        """
        rots = latent_variables_dict["R"].to(self.coords.device)
        in_shape = latent_variables_dict["R"].shape[:-2]
        z = None

        # sample conformations
        if self.z_dim > 0:
            if self.variational_conf:
                z = sample_conf(
                    latent_variables_dict["z"], latent_variables_dict["z_logvar"]
                )
            else:
                z = latent_variables_dict["z"]

        # generate slices
        device = self.coords.device
        if self.verbose_time:
            torch.cuda.synchronize(device)
        start_time_coords = time.time()
        x = self.coords[self.output_mask.binary_mask] @ rots.to(
            self.coords.device
        )  # batch_size(, n_tilts), n_pts, 3
        if self.verbose_time:
            torch.cuda.synchronize(device)
        start_time_query = time.time()
        y_pred = self.hypervolume(x, z)  # batch_size(, n_tilts), n_pts
        if self.verbose_time:
            torch.cuda.synchronize(device)
        end_time_query = time.time()
        times = {
            "coords": start_time_query - start_time_coords,
            "query": end_time_query - start_time_query,
        }

        # apply ctf
        y_pred = self.apply_ctf(y_pred, ctf_local)  # batch_size(, n_tilts), n_pts

        # apply translations (to gt)
        if not self.no_trans:
            trans = latent_variables_dict["t"][..., None, :].reshape(-1, 1, 2)
            y_gt_processed = self.lattice.translate_ht(
                y_gt.reshape(-1, self.lattice.D**2), trans, freqs2d=self.freqs2d
            ).reshape(*in_shape, -1)
            y_gt_processed = y_gt_processed[..., self.output_mask.binary_mask]
        else:
            y_gt_processed = y_gt.reshape(*in_shape, -1)
            y_gt_processed = y_gt_processed[..., self.output_mask.binary_mask]

        return y_pred, y_gt_processed, times, latent_variables_dict

    def eval_on_slice(self, x, z=None):
        """
        x: [batch_size, (nq, ) n_pts, 3]
        z: [batch_size, z_dim]

        output: [..., n_pts]
        """
        if x.dim() == 4:
            batch_size, nq, n_pts, _3 = x.shape
            x = x.reshape(batch_size, nq * n_pts, 3)
            y_pred = self.hypervolume(x, z)
            y_pred = y_pred.reshape(batch_size, nq, n_pts)
        else:
            y_pred = self.hypervolume(x, z)

        return y_pred

    def eval_volume(
        self,
        coords=None,
        resolution=None,
        extent=None,
        norm=None,
        zval=None,
        radius=None,
    ):
        use_coords = coords or self.lattice.coords
        use_resolution = resolution or self.lattice.D
        use_extent = extent or self.lattice.extent

        return self.hypervolume.eval_volume(
            coords=use_coords,
            resolution=use_resolution,
            extent=use_extent,
            norm=norm,
            zval=zval,
            radius=radius or self.output_mask.current_radius,
            z_dim=self.z_dim,
        )

    def apply_ctf(self, y_pred, ctf_local):
        """
        y_pred: [batch_size(, n_tilts), n_pts]
        ctf_local: [batch_size(, n_tilts), D, D]

        output: [batch_size(, n_tilts), n_pts]
        """
        ctf_local = ctf_local.reshape(*ctf_local.shape[:-2], -1)[
            ..., self.output_mask.binary_mask
        ]
        y_pred = ctf_local * y_pred
        return y_pred


def sample_conf(z_mu, z_logvar):
    """
    z_mu: [batch_size, z_dim]
    z_logvar: [batch_size, z_dim]

    output: [batch_size, z_dim]
    """
    # std = nn.Softplus(beta=2)(.5 * z_logvar)
    # std = nn.Softplus(beta=1)(z_logvar)
    std = torch.exp(0.5 * z_logvar)
    eps = torch.randn_like(std)
    z = eps * std + z_mu

    return z


class SharedCNN(nn.Module):
    def __init__(
        self,
        resolution,
        depth,
        channels,
        kernel_size,
        in_channels,
        nl=nn.ReLU,
        coord_conv=False,
        dropout=False,
        radial_average=False,
    ):
        """
        resolution: int
        depth: int
        channels: int
        kernel_size: int
        in_channels: int
        coord_conv: bool
        dropout: bool
        radial_average: bool
        """
        super(SharedCNN, self).__init__()

        cnn = []

        if radial_average:
            cnn.append(RadialAverager())
            final_size = resolution // 2
        else:
            final_size = resolution
        if coord_conv:
            cnn.append(AddCoords(final_size))
            in_channels = in_channels + 3
        else:
            in_channels = in_channels
        out_channels = channels
        for i in range(depth):
            ks = min(kernel_size, final_size)
            if dropout and i > 0:
                cnn.append(nn.Dropout2d())
            cnn.append(
                nn.Conv2d(
                    in_channels,
                    out_channels,
                    ks,
                    padding="same",
                    padding_mode="reflect",
                )
            )
            in_channels = out_channels
            cnn.append(nl())
            if 2 * in_channels <= 2048:
                out_channels = 2 * in_channels
            else:
                out_channels = in_channels
            if dropout:
                cnn.append(nn.Dropout2d())
            cnn.append(
                nn.Conv2d(
                    in_channels,
                    out_channels,
                    ks,
                    padding="same",
                    padding_mode="reflect",
                )
            )
            in_channels = out_channels
            cnn.append(nn.GroupNorm(channels, in_channels))
            if i < depth - 1:
                cnn.append(nl())
            else:
                cnn.append(nn.Tanh())
            if final_size // 2 > 0:
                cnn.append(nn.AvgPool2d(2))
                final_size = final_size // 2

        self.cnn = nn.Sequential(*cnn)

        self.final_size = final_size
        self.final_channels = in_channels

    def forward(self, y_real):
        """
        y_real: [..., d, D - 1, D - 1]

        output: [..., final_channels, final_size, final_size]
        """
        in_dims = y_real.shape[:-3]
        d = y_real.shape[-3]
        res = y_real.shape[-2]
        return self.cnn(y_real.reshape(np.prod(in_dims), d, res, res)).reshape(
            *in_dims, self.final_channels, self.final_size, self.final_size
        )


class RadialAverager(nn.Module):
    def __init__(self):
        super(RadialAverager, self).__init__()

    @staticmethod
    def forward(y_real):
        """
        y_real: [batch_size, d, D - 1, D - 1]

        output: [batch_size, d, (D - 1) // 2, (D - 1) // 2]
        """
        res = y_real.shape[-1]
        y_real_avg = torch.mean(
            torch.cat(
                [
                    y_real[..., None],
                    torch.flip(y_real, [-1, -2])[..., None],
                    torch.flip(torch.transpose(y_real, -2, -1), [-2])[..., None],
                    torch.flip(torch.transpose(y_real, -2, -1), [-1])[..., None],
                ],
                -1,
            ),
            -1,
        )
        return y_real_avg[..., : res // 2, : res // 2]


class AddCoords(nn.Module):
    def __init__(self, resolution, radius_channel=True):
        """
        resolution: int
        radius_channel: bool
        """
        super(AddCoords, self).__init__()
        self.radius_channel = radius_channel

        xx_ones = torch.ones([1, resolution], dtype=torch.int32)
        xx_ones = xx_ones.unsqueeze(-1)

        xx_range = torch.arange(resolution, dtype=torch.int32).unsqueeze(0)
        xx_range = xx_range.unsqueeze(1)

        xx_channel = torch.matmul(xx_ones, xx_range)
        xx_channel = xx_channel.unsqueeze(-1)

        yy_ones = torch.ones([1, resolution], dtype=torch.int32)
        yy_ones = yy_ones.unsqueeze(1)

        yy_range = torch.arange(resolution, dtype=torch.int32).unsqueeze(0)
        yy_range = yy_range.unsqueeze(-1)

        yy_channel = torch.matmul(yy_range, yy_ones)
        yy_channel = yy_channel.unsqueeze(-1)

        xx_channel = xx_channel.permute(0, 3, 1, 2)
        yy_channel = yy_channel.permute(0, 3, 1, 2)

        xx_channel = xx_channel.float() / (resolution - 1)
        yy_channel = yy_channel.float() / (resolution - 1)

        xx_channel = xx_channel - 0.5
        yy_channel = yy_channel - 0.5

        self.xx_channel = nn.Parameter(xx_channel, requires_grad=False)
        self.yy_channel = nn.Parameter(yy_channel, requires_grad=False)

        self.radius_calc = None
        if radius_channel:
            self.radius_calc = nn.Parameter(
                torch.sqrt(torch.pow(xx_channel, 2) + torch.pow(yy_channel, 2)),
                requires_grad=False,
            )

    def forward(self, x):
        """
        x: [batch_size, d, D - 1, D - 1]

        output: [batch_size, d + 2/3, D - 1, D - 1]
        """
        batch_size = x.shape[0]

        xx_channel = self.xx_channel.repeat(batch_size, 1, 1, 1)
        yy_channel = self.yy_channel.repeat(batch_size, 1, 1, 1)

        out = torch.cat([x, xx_channel, yy_channel], dim=1)

        if self.radius_channel:
            out = torch.cat([out, self.radius_calc.repeat(batch_size, 1, 1, 1)], dim=1)

        return out


class ConfTable(nn.Module):
    def __init__(self, n_imgs, z_dim, variational, std_z_init):
        """
        n_imgs: int
        z_dim: int
        variational: bool
        """
        super(ConfTable, self).__init__()
        self.variational = variational
        self.conf_init = torch.tensor(
            std_z_init * np.random.randn(n_imgs, z_dim)
        ).float()
        self.table_conf = nn.Parameter(self.conf_init, requires_grad=True)
        if variational:
            logvar_init = torch.tensor(np.ones((n_imgs, z_dim))).float()
            self.table_logvar = nn.Parameter(logvar_init, requires_grad=True)

    def initialize(self, conf):
        """
        conf: [n_imgs, z_dim] (numpy)
        """
        state_dict = self.state_dict()
        state_dict["table_conf"] = torch.tensor(conf).float()
        self.load_state_dict(state_dict)

    def forward(self, in_dict):
        """
        in_dict: dict
            index: [batch_size]
            y: [batch_size(, n_tilts), D, D]
            y_real: [batch_size(, n_tilts), D - 1, D - 1]
            R: [batch_size(, n_tilts), 3, 3]
            t: [batch_size(, n_tilts), 2]
            tilt_index: [batch_size( * n_tilts)]

        output: dict
            z: [batch_size, z_dim]
            z_logvar: [batch_size, z_dim] if variational and not pose_only
        """
        conf = self.table_conf[in_dict["indices"]]
        conf_dict = {"z": conf}
        if self.variational:
            logvar = self.table_logvar[in_dict["indices"]]
            conf_dict["z_logvar"] = logvar
        return conf_dict

    def reset(self):
        state_dict = self.state_dict()
        state_dict["table_conf"] = self.conf_init / 10.0
        self.load_state_dict(state_dict)


class PoseTable(nn.Module):
    def __init__(self, n_imgs, no_trans, resolution, use_gt_trans):
        """
        n_imgs: int
        no_trans: bool
        resolution: int
        use_gt_trans: bool
        """
        super(PoseTable, self).__init__()
        s2s2_init = torch.tensor(
            np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0])
            .reshape(1, 6)
            .repeat(n_imgs, axis=0)
        ).float()
        self.table_s2s2 = nn.Parameter(s2s2_init, requires_grad=True)
        self.no_trans = no_trans
        self.resolution = resolution
        self.use_gt_trans = use_gt_trans
        if not self.no_trans and not self.use_gt_trans:
            trans_init = torch.tensor(np.zeros((n_imgs, 2))).float()
            self.table_trans = nn.Parameter(trans_init, requires_grad=True)

    def initialize(self, rots, trans):
        """
        rots: [n_imgs, 3, 3] (numpy)
        trans: [n_imgs, 2] (numpy)
        """
        state_dict = self.state_dict()
        # rots must contain "corrected" rotations
        state_dict["table_s2s2"] = lie_tools.rotmat_to_s2s2(torch.tensor(rots).float())
        if "table_trans" in state_dict:
            # trans must be order 1
            state_dict["table_trans"] = torch.tensor(trans).float()
        self.load_state_dict(state_dict)

    def forward(self, in_dict):
        """
        in_dict: dict
            index: [batch_size]
            y: [batch_size(, n_tilts), D, D]
            y_real: [batch_size(, n_tilts), D - 1, D - 1]
            R: [batch_size(, n_tilts), 3, 3]
            t: [batch_size(, n_tilts), 2]
            tilt_index: [batch_size( * n_tilts)]

        output: dict
            R: [batch_size(, n_tilts), 3, 3]
            t: [batch_size(, n_tilts), 2]
        """
        rots_s2s2 = self.table_s2s2[in_dict["tilt_indices"]]
        rots_matrix = lie_tools.s2s2_to_rotmat(rots_s2s2)
        pose_dict = {"R": rots_matrix}
        if not self.no_trans:
            if not self.use_gt_trans:
                pose_dict["t"] = self.table_trans[in_dict["tilt_indices"]]
            else:
                pose_dict["t"] = in_dict["t"]
        if in_dict["y"].ndim == 4:
            pose_dict["R"] = pose_dict["R"].reshape(*in_dict["y"].shape[:-2], 3, 3)
            if not self.no_trans:
                pose_dict["t"] = pose_dict["t"].reshape(*in_dict["y"].shape[:-2], 2)
        return pose_dict


class ConfRegressor(nn.Module):
    def __init__(self, channels, kernel_size, z_dim, std_z_init, variational):
        """
        channels: int
        kernel_size: int
        z_dim: int
        std_z_init: float
        variational: bool
        """
        super(ConfRegressor, self).__init__()
        self.z_dim = z_dim
        self.variational = variational
        self.std_z_init = std_z_init
        if variational:
            out_features = 2 * z_dim
        else:
            out_features = z_dim
        self.out_features = out_features
        self.regressor = nn.Conv2d(channels, out_features, kernel_size, padding="valid")

    def forward(self, shared_features):
        """
        shared_features: [..., channels, kernel_size, kernel_size]

        output: dict
            z: [..., z_dim]
            z_logvar: [..., z_dim] if variational and not pose_only
        """
        in_dim = shared_features.shape[:-3]
        c = shared_features.shape[-3]
        ks = shared_features.shape[-2]
        z_full = self.regressor(shared_features.reshape(-1, c, ks, ks)).reshape(
            np.prod(in_dim), self.out_features
        )
        if self.variational:
            conf_dict = {
                "z": z_full[:, : self.z_dim],
                "z_logvar": nn.Tanh()(z_full[:, self.z_dim :] / 10.0) * 10.0,
            }
        else:
            conf_dict = {"z": z_full}
        return conf_dict


class HyperVolume(nn.Module):
    def __init__(
        self,
        resolution,
        z_dim,
        n_layers,
        hidden_dim,
        pe_type,
        pe_dim,
        feat_sigma,
        domain,
        pe_type_conf=None,
    ):
        """
        resolution: int
        z_dim: int
        n_layers: int
        hidden_dim: int
        pe_type: str
        pe_dim: int
        feat_sigma: float
        domain: str
        """
        super(HyperVolume, self).__init__()
        self.pe_type = pe_type
        self.pe_dim = pe_dim
        if pe_type == "gaussian":
            rand_freqs = torch.randn((3 * pe_dim, 3), dtype=torch.float) * feat_sigma
            self.rand_freqs = nn.Parameter(rand_freqs, requires_grad=False)
            x_pe_dim = 3 * 2 * pe_dim
        else:
            raise NotImplementedError
        self.pe_type_conf = pe_type_conf
        if pe_type_conf is None:
            z_pe_dim = z_dim
        elif pe_type_conf == "geom":
            min_freq = -4
            n_freqs = 4
            geom_freqs_conf = (
                2.0 ** torch.arange(min_freq, min_freq + n_freqs, dtype=torch.float)
                * np.pi
            )
            self.geom_freqs_conf = nn.Parameter(geom_freqs_conf, requires_grad=False)
            z_pe_dim = z_dim * 2 * n_freqs
        else:
            raise NotImplementedError

        self.D = resolution
        self.z_dim = z_dim
        self.n_layers = n_layers
        self.hidden_dim = hidden_dim
        self.feat_sigma = feat_sigma
        self.domain = domain

        in_features = x_pe_dim + z_pe_dim
        if domain == "hartley":
            self.mlp = ResidualLinearMLP(in_features, n_layers, hidden_dim, 1)
        else:
            raise NotImplementedError

    def forward(self, x, z):
        """
        x: [batch_size(, n_tilts), n_pts, 3]
        z: [batch_size, z_dim] or None

        output: [batch_size(, n_tilts), n_pts]
        """
        batch_size_in = x.shape[0]
        n_pts = x.shape[-2]
        subtomogram_averaging = x.dim() == 4
        if self.pe_type == "gaussian":
            x = self.random_fourier_encoding(x)
        if z is not None:
            if self.pe_type_conf == "geom":
                z = self.geom_fourier_encoding_conf(z)
            if subtomogram_averaging:
                n_tilts = x.shape[1]
                z_expand = z[:, None, None].expand(-1, n_tilts, n_pts, -1)
            else:
                z_expand = z[:, None].expand(-1, n_pts, -1)
            x = torch.cat([x, z_expand], -1)
        if subtomogram_averaging:
            n_tilts = x.shape[1]
            out_shape = (batch_size_in, n_tilts, n_pts)
        else:
            out_shape = (batch_size_in, n_pts)
        return self.mlp(x).reshape(*out_shape)

    def random_fourier_encoding(self, x):
        """
        x: [batch_size(, n_tilts), n_pts, 3]

        output: [batch_size(, n_tilts), n_pts, 3 * 2 * pe_dim]
        """
        freqs = self.rand_freqs.reshape(1, 1, -1, 3) * (self.D // 2)
        kx_ky_kz = x[..., None, :] * freqs
        k = kx_ky_kz.sum(-1)
        s = torch.sin(k)
        c = torch.cos(k)
        x_encoded = torch.cat([s, c], -1)
        return x_encoded

    def geom_fourier_encoding_conf(self, z):
        """
        z: [batch_size, z_dim]

        output: [batch_size, z_dim * 2 * pe_dim]
        """
        in_dims = z.shape[:-1]
        s = torch.sin(z[..., None] * self.geom_freqs_conf)  # [..., z_dim, pe_dim]
        c = torch.cos(z[..., None] * self.geom_freqs_conf)  # [..., z_dim, pe_dim]
        z_encoded = torch.cat([s, c], -1).reshape(*in_dims, -1)
        return z_encoded

    def get_building_params(self):
        building_params = {
            "resolution": self.D,
            "z_dim": self.z_dim,
            "n_layers": self.n_layers,
            "hidden_dim": self.hidden_dim,
            "pe_type": self.pe_type,
            "pe_dim": self.pe_dim,
            "feat_sigma": self.feat_sigma,
            "domain": self.domain,
            "pe_type_conf": self.pe_type_conf,
        }
        return building_params

    def eval_volume(
        self,
        coords,
        resolution,
        extent,
        norm,
        zval=None,
        radius=None,
        z_dim=None,
    ):
        """
        lattice: Lattice
        z_dim: int
        norm: (mean, std)
        zval: [z_dim]
        radius: int
        """
        z_dim = z_dim or self.z_dim
        radius_normalized = extent * 2 * radius / resolution

        z = None
        if zval is not None:
            z = torch.tensor(zval, dtype=torch.float32, device=coords.device).reshape(
                1, z_dim
            )

        volume = np.zeros((resolution, resolution, resolution), dtype=np.float32)
        assert not self.training
        dzs = np.linspace(-extent, extent, resolution, endpoint=True, dtype=np.float32)

        with torch.no_grad():
            for i, dz in enumerate(dzs):
                x = coords + torch.tensor([0, 0, dz], device=coords.device)
                x = x.reshape(1, -1, 3)
                y = self(x, z)

                slice_radius = int(
                    np.sqrt(max(radius_normalized**2 - dz**2, 0.0)) * resolution
                )
                slice_mask = masking.CircularMask(
                    slice_radius, coords, resolution, extent
                ).binary_mask

                y[0, ~slice_mask] = 0.0
                y = y.view(resolution, resolution).detach().cpu().numpy()
                volume[i] = y

            # remove last +k freq for inverse FFT
            volume = volume * norm[1] + norm[0]
            volume_real = fft.ihtn_center(torch.Tensor(volume[0:-1, 0:-1, 0:-1]))

        return volume_real


class VolumeExplicit(nn.Module):
    def __init__(self, resolution, domain, extent):
        """
        resolution: int
        domain: str
        extent: float
        """
        super(VolumeExplicit, self).__init__()
        assert domain == "hartley"
        self.D = resolution
        self.domain = domain
        self.extent = extent

        self.volume = nn.Parameter(
            1e-5
            * torch.tensor(np.random.randn(resolution, resolution, resolution)).float(),
            requires_grad=True,
        )

    def forward(self, x, z):
        """
        x: [batch_size, n_pts, 3] in [-extent, extent]
        z: None

        output: [batch_size, n_pts]
        """
        assert (
            z is None
        ), "Explicit volume(s) do not support heterogeneous reconstruction."
        batch_size_in = x.shape[0]
        out = torch.nn.functional.grid_sample(
            1e2 * self.volume[None, None].repeat(batch_size_in, 1, 1, 1, 1),
            x[:, None, None, :, :] / (2.0 * self.extent) * 2,
            mode="bilinear",
            padding_mode="zeros",
            align_corners=False,
        )

        return out.reshape(batch_size_in, -1)

    def get_building_params(self):
        building_params = {
            "resolution": self.D,
            "domain": self.domain,
            "extent": self.extent,
        }

        return building_params


class ResidualLinear(nn.Module):
    def __init__(self, n_in, n_out):
        super(ResidualLinear, self).__init__()
        self.linear = nn.Linear(n_in, n_out)

    def forward(self, x):
        z = self.linear(x) + x
        return z


class MyLinear(nn.Linear):
    def forward(self, x):
        if x.dtype == torch.half:
            return half_linear(x, self.weight, self.bias)
        else:
            return single_linear(x, self.weight, self.bias)


class ResidualLinearMLP(nn.Module):
    def __init__(self, in_dim, n_layers, hidden_dim, out_dim, nl=nn.ReLU):
        super(ResidualLinearMLP, self).__init__()
        layers = [
            ResidualLinear(in_dim, hidden_dim)
            if in_dim == hidden_dim
            else nn.Linear(in_dim, hidden_dim),
            nl(),
        ]
        for n in range(n_layers):
            layers.append(ResidualLinear(hidden_dim, hidden_dim))
            layers.append(nl())
        layers.append(
            ResidualLinear(hidden_dim, out_dim)
            if out_dim == hidden_dim
            else MyLinear(hidden_dim, out_dim)
        )
        self.main = nn.Sequential(*layers)

    def forward(self, x):
        """
        x: [..., in_dim]

        output: [..., out_dim]
        """
        flat = x.view(-1, x.shape[-1])
        ret_flat = self.main(flat)
        ret = ret_flat.view(*x.shape[:-1], ret_flat.shape[-1])
        return ret


class MyDataParallel(nn.DataParallel):
    def __getattr__(self, name):
        try:
            return super().__getattr__(name)
        except AttributeError:
            return getattr(self.module, name)
