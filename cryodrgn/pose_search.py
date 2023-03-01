import logging
import numpy as np
import torch
import torch.nn.functional as F
from typing import Optional, Union, Tuple
from cryodrgn import lie_tools, shift_grid, so3_grid
from cryodrgn.models import unparallelize, HetOnlyVAE
from cryodrgn.lattice import Lattice
import torch.nn as nn

logger = logging.getLogger(__name__)


def rot_2d(angle: float, outD: int, device: torch.device) -> torch.Tensor:
    rot = torch.zeros((outD, outD), device=device)
    rot[0, 0] = np.cos(angle)
    rot[0, 1] = -np.sin(angle)
    rot[1, 0] = np.sin(angle)
    rot[1, 1] = np.cos(angle)
    return rot


def to_tensor(x: Union[np.ndarray, torch.Tensor, None]):
    if isinstance(x, np.ndarray):
        x = torch.from_numpy(x)
    return x


def interpolate(img: torch.Tensor, coords: torch.Tensor) -> torch.Tensor:
    # print(f"Interpolating {img.shape} {coords.shape}")
    assert len(coords.shape) == 2
    assert coords.shape[-1] == 2
    grid = coords * 2  # careful here! grid_sample expects [-1,1] instead of [-0.5,0.5]
    grid = grid[None, None, ...].expand(img.shape[0], -1, -1, -1)

    res = (
        F.grid_sample(
            img.unsqueeze(1),
            grid,
        )
        .squeeze(2)
        .squeeze(1)
    )

    return res


FAST_INPLANE = True


class PoseSearch:
    """Pose search"""

    def __init__(
        self,
        model: nn.Module,
        lattice: Lattice,
        Lmin: int,
        Lmax: int,
        tilt=None,
        base_healpy: int = 1,
        t_extent: int = 5,
        t_ngrid: int = 7,
        niter: int = 5,
        nkeptposes: int = 24,
        loss_fn: str = "msf",
        t_xshift: int = 0,
        t_yshift: int = 0,
        device: Optional[torch.device] = None,
    ):
        self.model = model
        self.lattice = lattice
        self.base_healpy = base_healpy
        self.so3_base_quat = so3_grid.grid_SO3(base_healpy)
        self.base_quat = (
            so3_grid.s2_grid_SO3(base_healpy) if FAST_INPLANE else self.so3_base_quat
        )
        self.so3_base_rot = lie_tools.quaternions_to_SO3(
            to_tensor(self.so3_base_quat)
        ).to(device)
        self.base_rot = lie_tools.quaternions_to_SO3(to_tensor(self.base_quat)).to(
            device
        )

        self.nbase = len(self.base_quat)
        self.base_inplane = so3_grid.grid_s1(base_healpy)
        self.base_shifts = torch.tensor(
            shift_grid.base_shift_grid(
                base_healpy - 1, t_extent, t_ngrid, xshift=t_xshift, yshift=t_yshift
            ),
            device=device,
        ).float()
        self.t_extent = t_extent
        self.t_ngrid = t_ngrid

        self.Lmin = Lmin
        self.Lmax = Lmax
        self.niter = niter
        self.tilt = tilt
        self.nkeptposes = nkeptposes
        self.loss_fn = loss_fn
        self._so3_neighbor_cache = {}  # for memoization
        self._shift_neighbor_cache = {}  # for memoization

        self.device = device

    def eval_grid(
        self,
        *,
        images: torch.Tensor,
        rot: torch.Tensor,
        z: Optional[torch.Tensor],
        NQ: int,
        L: int,
        images_tilt: Optional[torch.Tensor] = None,
        angles_inplane: Optional[np.ndarray] = None,
        ctf_i: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        images: B x T x Npix
        rot: (NxQ) x 3 x 3 rotation matrics (N=1 for base grid, N=B for incremental grid)
        NQ: number of slices evaluated for each image
        L: radius of fourier components to evaluate
        """
        B = images.size(0)
        mask = self.lattice.get_circular_mask(L)
        coords = self.lattice.coords[mask]  # .to(rot.device)
        YX = coords.size(-2)
        device = next(self.model.parameters()).device
        if ctf_i is not None:
            ctf_i = ctf_i.view(B, 1, 1, -1)[..., mask]  # Bx1x1xYX

        def compute_err(images, rot):
            adj_angles_inplane = None
            if angles_inplane is not None:
                # apply a random in-plane rotation from the set
                # to avoid artifacts due to grid alignment
                rand_a = angles_inplane[np.random.randint(len(angles_inplane))]
                rand_inplane_rot = rot_2d(rand_a, 3, rot.device)
                rot = rand_inplane_rot @ rot
                adj_angles_inplane = angles_inplane - rand_a

            x = coords @ rot
            if z is not None:
                _model = unparallelize(self.model)
                assert isinstance(_model, HetOnlyVAE)
                x = _model.cat_z(x, z)
            x = x.to(device)
            # logger.info(f"Evaluating model on {x.shape} = {x.nelement() // 3} points")
            with torch.no_grad():
                y_hat = self.model(x)
                y_hat = y_hat.float()
            y_hat = y_hat.view(
                -1, 1, NQ, YX
            )  # 1x1xNQxYX for base grid, Bx1x8xYX for incremental grid
            if ctf_i is not None:
                y_hat = y_hat * ctf_i
            if adj_angles_inplane is not None:
                y_hat = self.rotate_images(y_hat, adj_angles_inplane, L)
            images = images.unsqueeze(2)  # BxTx1xYX
            if self.loss_fn == "mse":
                err = (images - y_hat).pow(2).sum(-1)  # BxTxQ
            elif self.loss_fn == "msf":
                B, T, _, Npix = images.shape
                Npix = images.shape[-1]
                dots = (
                    images.view(B, -1, Npix)
                    @ y_hat.view(y_hat.shape[0], -1, Npix).transpose(-1, -2)
                ).view(B, T, -1)
                norm = (y_hat * y_hat).sum(-1) / 2

                err = -dots + norm  # BxTxQ

                # err1 = -(images * y_hat).sum(-1) + (y_hat * y_hat).sum(-1) / 2   # BxTxQ
                # delta = (err1 - err).abs().max() / (err1 + err).mean() < 1e-2
            elif self.loss_fn == "cor":
                err = -(images * y_hat).sum(-1) / y_hat.std(-1)
            else:
                raise NotImplementedError(f"Unknown loss_fn: {self.loss_fn}")
            return err

        err = compute_err(images, rot)
        if images_tilt is not None:
            err_tilt = compute_err(images_tilt, self.tilt @ rot)
            err += err_tilt
        return err  # BxTxQ

    def mask_images(self, images, L):
        """
        images: B x NY x NX x 2
        Returns: B x Npix at resolution L
        """
        B = images.size(0)
        mask = self.lattice.get_circular_mask(L)
        return images.view(B, -1)[:, mask]

    def translate_images(
        self, images: torch.Tensor, shifts: torch.Tensor, L: int
    ) -> torch.Tensor:
        """
        images: B x NY x NX
        shifts: B x T x 2 or B
        Returns: B x T x Npix at resolution L
        """

        B = images.size(0)
        mask = self.lattice.get_circular_mask(L)
        res = self.lattice.translate_ht(images.view(B, -1)[:, mask], shifts, mask)

        return res

    def rotate_images(
        self, images: torch.Tensor, angles: np.ndarray, L: int
    ) -> torch.Tensor:
        B, d1, NQ, YX = images.shape
        BNQ = B * NQ
        squeezed_images = images.view(BNQ, YX)
        D = self.lattice.D
        res = torch.zeros((B * NQ, len(angles), YX), device=images.device)
        # B x NQ x YX
        mask = self.lattice.get_circular_mask(L)

        rot_matrices = torch.stack([rot_2d(a, 2, images.device) for a in angles], dim=0)
        lattice_coords = self.lattice.coords[mask][:, :2]
        rot_coords = lattice_coords @ rot_matrices

        full_images = torch.zeros((BNQ, D, D), device=images.device)
        full_images.view(BNQ, D * D)[:, mask] = squeezed_images

        for angle_idx, interp_coords in enumerate(rot_coords):
            interpolated = interpolate(full_images, interp_coords)
            assert squeezed_images.shape == interpolated.shape
            # IMPORTANT TRICK HERE!
            interpolated *= squeezed_images.std(-1, keepdim=True) / interpolated.std(
                -1, keepdim=True
            )  # FIXME

            res[:, angle_idx] = interpolated

        return res.view(B, 1, NQ * len(angles), YX)

    def get_neighbor_so3(self, quat: np.ndarray, s2i: int, s1i: int, res: int):
        """Memoization of so3_grid.get_neighbor."""
        key = (int(s2i), int(s1i), int(res))
        if key not in self._so3_neighbor_cache:
            self._so3_neighbor_cache[key] = so3_grid.get_neighbor(quat, s2i, s1i, res)
        # FIXME: will this cache get too big? maybe don't do it when res is too
        return self._so3_neighbor_cache[key]

    def get_neighbor_shift(self, x, y, res):
        """Memoization of shift_grid.get_neighbor."""
        key = (int(x), int(y), int(res))
        if key not in self._shift_neighbor_cache:
            self._shift_neighbor_cache[key] = shift_grid.get_neighbor(
                x, y, res - 1, self.t_extent, self.t_ngrid
            )
        # FIXME: will this cache get too big? maybe don't do it when res is too
        return self._shift_neighbor_cache[key]

    def subdivide(
        self, quat: np.ndarray, q_ind: np.ndarray, cur_res: int
    ) -> Tuple[np.ndarray, np.ndarray, torch.Tensor]:
        """
        Subdivides poses for next resolution level

        Inputs:
            quat (N x 4 tensor): quaternions
            q_ind (N x 2 np.array): index of current S2xS1 grid
            cur_res (int): Current resolution level

        Returns:
            quat  (N x 8 x 4) np.array
            q_ind (N x 8 x 2) np.array
            rot   (N*8 x 3 x 3) tensor
        """
        N = quat.shape[0]

        assert len(quat.shape) == 2 and quat.shape == (N, 4), quat.shape
        assert len(q_ind.shape) == 2 and q_ind.shape == (N, 2), q_ind.shape

        # get neighboring SO3 elements at next resolution level -- todo: make this an array operation
        neighbors = [
            self.get_neighbor_so3(quat[i], q_ind[i][0], q_ind[i][1], cur_res)
            for i in range(len(quat))
        ]
        quat = np.array([x[0] for x in neighbors])  # Bx8x4
        q_ind = np.array([x[1] for x in neighbors])  # Bx8x2
        rot = lie_tools.quaternions_to_SO3(torch.from_numpy(quat).view(-1, 4)).to(
            self.device
        )

        assert len(quat.shape) == 3 and quat.shape == (N, 8, 4), quat.shape
        assert len(q_ind.shape) == 3 and q_ind.shape == (N, 8, 2), q_ind.shape
        assert len(rot.shape) == 3 and rot.shape == (N * 8, 3, 3), rot.shape

        return quat, q_ind, rot

    def keep_matrix(self, loss: torch.Tensor, B: int, max_poses: int) -> torch.Tensor:
        """
        Inputs:
            loss (B, T, Q): tensor of losses for each translation and rotation.

        Returns:
            keep (3, B * max_poses): bool tensor of rotations to keep, along with the best translation for each
        """
        shape = loss.shape
        assert len(shape) == 3
        best_loss, best_trans_idx = loss.min(1)
        flat_loss = best_loss.view(B, -1)
        flat_idx = flat_loss.topk(max_poses, dim=-1, largest=False, sorted=True)[1]
        # add the batch index in, to make it completely flat
        flat_idx += (
            torch.arange(B, device=loss.device).unsqueeze(1) * flat_loss.shape[1]
        )
        flat_idx = flat_idx.view(-1)

        keep_idx = torch.empty(
            len(shape), B * max_poses, dtype=torch.long, device=loss.device
        )
        keep_idx[0] = flat_idx // shape[2]
        keep_idx[2] = flat_idx % shape[2]
        keep_idx[1] = best_trans_idx[keep_idx[0], keep_idx[2]]
        return keep_idx

    def getL(self, iter_: int) -> int:
        L = self.Lmin + int(iter_ / self.niter * (self.Lmax - self.Lmin))
        return min(L, self.lattice.D // 2)
        # return min(self.Lmin * 2 ** iter_, self.Lmax)

    def opt_theta_trans(
        self,
        images: torch.Tensor,
        z: Optional[torch.Tensor] = None,
        images_tilt: Optional[torch.Tensor] = None,
        init_poses: Optional[torch.Tensor] = None,
        ctf_i=None,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        images = to_tensor(images)  # type: ignore
        images_tilt = to_tensor(images_tilt)
        init_poses = to_tensor(init_poses)
        z = to_tensor(z)
        device = images.device
        do_tilt = images_tilt is not None

        B = images.size(0)
        assert not self.model.training

        loss = rot = None
        if init_poses is None:
            # Expand the base grid B times if each image has a different z
            if z is not None:
                base_rot = self.base_rot.expand(
                    B, *self.base_rot.shape
                )  # B x 576 x 3 x 3
            else:
                base_rot = self.base_rot  # 576 x 3 x 3
            base_rot = base_rot.to(device)
            # Compute the loss for all poses
            L = self.getL(0)
            loss = self.eval_grid(
                images=self.translate_images(images, self.base_shifts, L),
                rot=base_rot,
                z=z,
                NQ=self.nbase,
                L=L,
                images_tilt=self.translate_images(images_tilt, self.base_shifts, L)
                if do_tilt
                else None,
                angles_inplane=self.base_inplane if FAST_INPLANE else None,
                ctf_i=ctf_i,
            )
            keepB, keepT, keepQ = self.keep_matrix(
                loss, B, self.nkeptposes
            ).cpu()  # B x -1
        else:
            # careful, overwrite the old batch index which is now invalid
            keepB = (
                torch.arange(B, device=init_poses.device)
                .unsqueeze(1)
                .repeat(1, self.nkeptposes)
                .view(-1)
            )
            keepT, keepQ = init_poses.reshape(-1, 2).t()

        new_init_poses = (
            torch.cat((keepT, keepQ), dim=-1)
            .view(2, B, self.nkeptposes)
            .permute(1, 2, 0)
        )

        quat = self.so3_base_quat[keepQ]
        q_ind = so3_grid.get_base_ind(keepQ, self.base_healpy)  # Np x 2
        trans = self.base_shifts[keepT]
        shifts = self.base_shifts.clone()
        for iter_ in range(1, self.niter + 1):
            keepB8 = (
                keepB.unsqueeze(1).repeat(1, 8).view(-1)
            )  # repeat each element 8 times
            zb = z[keepB8] if z is not None else None

            L = self.getL(iter_)
            quat, q_ind, rot = self.subdivide(quat, q_ind, iter_ + self.base_healpy - 1)
            shifts /= 2
            trans = trans.unsqueeze(1) + shifts.unsqueeze(0)  # FIXME: scale
            rot = rot.to(device)
            loss = self.eval_grid(
                images=self.translate_images(
                    images[keepB], trans, L
                ),  # (B*24, 4, Npoints)
                rot=rot,
                z=zb,
                NQ=8,
                L=L,
                images_tilt=self.translate_images(images_tilt[keepB], trans, L)
                if do_tilt
                else None,  # (B*24, 4, Npoints)
                ctf_i=ctf_i[keepB] if ctf_i is not None else ctf_i,
            )  # sum(NP), 8

            # nkeptposes = 1
            # nkeptposes = max(1, math.ceil(self.nkeptposes / 2 ** (iter_-1)))
            nkeptposes = self.nkeptposes if iter_ < self.niter else 1

            keepBN, keepT, keepQ = self.keep_matrix(
                loss, B, nkeptposes
            ).cpu()  # B x (self.Nkeptposes*32)
            keepB = keepBN * B // loss.shape[0]  # FIXME: expain
            assert (
                len(keepB) == B * nkeptposes
            ), f"{len(keepB)} != {B} x {nkeptposes} at iter {iter_}"
            quat = quat[keepBN, keepQ]
            q_ind = q_ind[keepBN, keepQ]
            trans = trans[keepBN, keepT]

        assert loss is not None
        bestBN, bestT, bestQ = self.keep_matrix(loss, B, 1).cpu()
        assert len(bestBN) == B
        if self.niter == 0:
            best_rot = self.so3_base_rot[bestQ].to(device)
            best_trans = self.base_shifts[bestT].to(device)
        else:
            assert rot is not None
            best_rot = rot.view(-1, 8, 3, 3)[bestBN, bestQ]
            best_trans = trans.to(device)

        return best_rot, best_trans, new_init_poses
