"""Pytorch models"""
from typing import Optional, Tuple, Type, Sequence, Any
import numpy as np
import torch
from torch import Tensor
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.parameter import Parameter
from torch.nn.parallel import DataParallel
from cryodrgn import fft, lie_tools, utils
import cryodrgn.config
from cryodrgn.lattice import Lattice

Norm = Sequence[Any]  # mean, std


def unparallelize(model: nn.Module) -> nn.Module:
    if isinstance(model, DataParallelDecoder):
        return model.dp.module
    if isinstance(model, DataParallel):
        return model.module
    return model


class HetOnlyVAE(nn.Module):
    # No pose inference
    def __init__(
        self,
        lattice: Lattice,
        qlayers: int,
        qdim: int,
        players: int,
        pdim: int,
        in_dim: int,
        zdim: int = 1,
        encode_mode: str = "resid",
        enc_mask=None,
        enc_type="linear_lowf",
        enc_dim=None,
        domain="fourier",
        activation=nn.ReLU,
        feat_sigma: Optional[float] = None,
        tilt_params={},
    ):
        super(HetOnlyVAE, self).__init__()
        self.lattice = lattice
        self.zdim = zdim
        self.in_dim = in_dim
        self.enc_mask = enc_mask
        if encode_mode == "conv":
            self.encoder = ConvEncoder(qdim, zdim * 2)
        elif encode_mode == "resid":
            self.encoder = ResidLinearMLP(
                in_dim,
                qlayers,
                qdim,
                zdim * 2,
                activation,  # nlayers  # hidden_dim  # out_dim
            )
        elif encode_mode == "mlp":
            self.encoder = MLP(
                in_dim, qlayers, qdim, zdim * 2, activation  # hidden_dim  # out_dim
            )  # in_dim -> hidden_dim
        elif encode_mode == "tilt":
            self.encoder = TiltEncoder(
                in_dim,
                qlayers,
                qdim,
                tilt_params["t_emb_dim"],  # embedding dim
                tilt_params["ntilts"],  # number of encoded tilts
                tilt_params["tlayers"],
                tilt_params["tdim"],
                zdim * 2,  # outdim
                activation,
            )
        else:
            raise RuntimeError("Encoder mode {} not recognized".format(encode_mode))
        self.encode_mode = encode_mode
        self.decoder = get_decoder(
            3 + zdim,
            lattice.D,
            players,
            pdim,
            domain,
            enc_type,
            enc_dim,
            activation,
            feat_sigma,
        )

    @classmethod
    def load(cls, config, weights=None, device=None):
        """Instantiate a model from a config.yaml

        Inputs:
            config (str, dict): Path to config.yaml or loaded config.yaml
            weights (str): Path to weights.pkl
            device: torch.device object

        Returns:
            HetOnlyVAE instance, Lattice instance
        """
        cfg = cryodrgn.config.load(config)

        c = cfg["lattice_args"]
        lat = Lattice(c["D"], extent=c["extent"], device=device)
        c = cfg["model_args"]
        if c["enc_mask"] > 0:
            enc_mask = lat.get_circular_mask(c["enc_mask"])
            in_dim = int(enc_mask.sum())
        else:
            assert c["enc_mask"] == -1
            enc_mask = None
            in_dim = lat.D**2
        activation = {"relu": nn.ReLU, "leaky_relu": nn.LeakyReLU}[c["activation"]]
        model = HetOnlyVAE(
            lat,
            c["qlayers"],
            c["qdim"],
            c["players"],
            c["pdim"],
            in_dim,
            c["zdim"],
            encode_mode=c["encode_mode"],
            enc_mask=enc_mask,
            enc_type=c["pe_type"],
            enc_dim=c["pe_dim"],
            domain=c["domain"],
            activation=activation,
            feat_sigma=c["feat_sigma"],
            tilt_params=c.get("tilt_params", {}),
        )
        if weights is not None:
            ckpt = torch.load(weights, map_location=device)
            model.load_state_dict(ckpt["model_state_dict"])
        if device is not None:
            model.to(device)
        return model, lat

    def reparameterize(self, mu, logvar):
        if not self.training:
            return mu
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return eps * std + mu

    def encode(self, *img) -> Tuple[Tensor, Tensor]:
        img = (x.view(x.shape[0], -1) for x in img)
        if self.enc_mask is not None:
            img = (x[:, self.enc_mask] for x in img)
        z = self.encoder(*img)
        return z[:, : self.zdim], z[:, self.zdim :]

    def cat_z(self, coords, z) -> Tensor:
        """
        coords: Bx...x3
        z: Bxzdim
        """
        assert coords.size(0) == z.size(0), (coords.shape, z.shape)
        z = z.view(z.size(0), *([1] * (coords.ndimension() - 2)), self.zdim)
        z = torch.cat((coords, z.expand(*coords.shape[:-1], self.zdim)), dim=-1)
        return z

    def decode(self, coords, z=None) -> torch.Tensor:
        """
        coords: BxNx3 image coordinates
        z: Bxzdim latent coordinate
        """
        decoder = self.decoder
        assert isinstance(decoder, nn.Module)
        retval = decoder(self.cat_z(coords, z) if z is not None else coords)
        return retval

    # Need forward func for DataParallel -- TODO: refactor
    def forward(self, *args, **kwargs):
        return self.decode(*args, **kwargs)


def load_decoder(config, weights=None, device=None):
    """
    Instantiate a decoder model from a config.yaml

    Inputs:
        config (str, dict): Path to config.yaml or loaded config.yaml
        weights (str): Path to weights.pkl
        device: torch.device object

    Returns a decoder model
    """
    cfg = cryodrgn.config.load(config)
    c = cfg["model_args"]
    D = cfg["lattice_args"]["D"]
    activation = {"relu": nn.ReLU, "leaky_relu": nn.LeakyReLU}[c["activation"]]
    model = get_decoder(
        3,
        D,
        c["layers"],
        c["dim"],
        c["domain"],
        c["pe_type"],
        c["pe_dim"],
        activation,
        c["feat_sigma"],
    )
    if weights is not None:
        ckpt = torch.load(weights)
        model.load_state_dict(ckpt["model_state_dict"])
    if device is not None:
        model.to(device)
    return model


class Decoder(nn.Module):
    def eval_volume(
        self,
        coords: Tensor,
        D: int,
        extent: float,
        norm: Norm,
        zval: Optional[np.ndarray] = None,
    ) -> Tensor:
        """
        Evaluate the model on a DxDxD volume
        Inputs:
            coords: lattice coords on the x-y plane (D^2 x 3)
            D: size of lattice
            extent: extent of lattice [-extent, extent]
            norm: data normalization
            zval: value of latent (zdim x 1)
        """
        raise NotImplementedError

    def get_voxel_decoder(self) -> Optional["Decoder"]:
        return None


class DataParallelDecoder(Decoder):
    def __init__(self, decoder: Decoder):
        super(DataParallelDecoder, self).__init__()
        self.dp = torch.nn.parallel.DataParallel(decoder)

    def eval_volume(self, *args, **kwargs):
        module = self.dp.module
        assert isinstance(module, Decoder)
        return module.eval_volume(*args, **kwargs)

    def forward(self, *args, **kwargs):
        return self.dp.module.forward(*args, **kwargs)

    def state_dict(self, *args, **kwargs):
        return self.dp.module.state_dict(*args, **kwargs)


class PositionalDecoder(Decoder):
    def __init__(
        self,
        in_dim,
        D,
        nlayers,
        hidden_dim,
        activation,
        enc_type="linear_lowf",
        enc_dim=None,
        feat_sigma: Optional[float] = None,
    ):
        super(PositionalDecoder, self).__init__()
        assert in_dim >= 3
        self.zdim = in_dim - 3
        self.D = D
        self.D2 = D // 2
        self.DD = 2 * (D // 2)
        self.enc_dim = self.D2 if enc_dim is None else enc_dim
        self.enc_type = enc_type
        self.in_dim = 3 * (self.enc_dim) * 2 + self.zdim
        self.decoder = ResidLinearMLP(self.in_dim, nlayers, hidden_dim, 1, activation)

        if enc_type == "gaussian":
            # We construct 3 * self.enc_dim random vector frequences, to match the original positional encoding:
            # In the positional encoding we produce self.enc_dim features for each of the x,y,z dimensions,
            # whereas in gaussian encoding we produce self.enc_dim features each with random x,y,z components
            #
            # Each of the random feats is the sine/cosine of the dot product of the coordinates with a frequency
            # vector sampled from a gaussian with std of feat_sigma
            rand_freqs = (
                torch.randn((3 * self.enc_dim, 3), dtype=torch.float) * feat_sigma
            )
            # make rand_feats a parameter so it is saved in the checkpoint, but do not perform SGD on it
            self.rand_freqs = Parameter(rand_freqs, requires_grad=False)
        else:
            self.rand_feats = None

    def positional_encoding_geom(self, coords):
        """Expand coordinates in the Fourier basis with geometrically spaced wavelengths from 2/D to 2pi"""
        if self.enc_type == "gaussian":
            return self.random_fourier_encoding(coords)
        freqs = torch.arange(self.enc_dim, dtype=torch.float, device=coords.device)
        if self.enc_type == "geom_ft":
            freqs = (
                self.DD * np.pi * (2.0 / self.DD) ** (freqs / (self.enc_dim - 1))
            )  # option 1: 2/D to 1
        elif self.enc_type == "geom_full":
            freqs = (
                self.DD
                * np.pi
                * (1.0 / self.DD / np.pi) ** (freqs / (self.enc_dim - 1))
            )  # option 2: 2/D to 2pi
        elif self.enc_type == "geom_lowf":
            freqs = self.D2 * (1.0 / self.D2) ** (
                freqs / (self.enc_dim - 1)
            )  # option 3: 2/D*2pi to 2pi
        elif self.enc_type == "geom_nohighf":
            freqs = self.D2 * (2.0 * np.pi / self.D2) ** (
                freqs / (self.enc_dim - 1)
            )  # option 4: 2/D*2pi to 1
        elif self.enc_type == "linear_lowf":
            return self.positional_encoding_linear(coords)
        else:
            raise RuntimeError("Encoding type {} not recognized".format(self.enc_type))
        freqs = freqs.view(*[1] * len(coords.shape), -1)  # 1 x 1 x D2
        coords = coords.unsqueeze(-1)  # B x 3 x 1
        k = coords[..., 0:3, :] * freqs  # B x 3 x D2
        s = torch.sin(k)  # B x 3 x D2
        c = torch.cos(k)  # B x 3 x D2
        x = torch.cat([s, c], -1)  # B x 3 x D
        x = x.view(*coords.shape[:-2], self.in_dim - self.zdim)  # B x in_dim-zdim
        if self.zdim > 0:
            x = torch.cat([x, coords[..., 3:, :].squeeze(-1)], -1)
            assert x.shape[-1] == self.in_dim
        return x

    def random_fourier_encoding(self, coords):
        assert self.rand_freqs is not None
        # k = coords . rand_freqs
        # expand rand_freqs with singleton dimension along the batch dimensions
        # e.g. dim (1, ..., 1, n_rand_feats, 3)
        freqs = self.rand_freqs.view(*[1] * (len(coords.shape) - 1), -1, 3) * self.D2

        kxkykz = coords[..., None, 0:3] * freqs  # compute the x,y,z components of k
        k = kxkykz.sum(-1)  # compute k
        s = torch.sin(k)
        c = torch.cos(k)
        x = torch.cat([s, c], -1)
        x = x.view(*coords.shape[:-1], self.in_dim - self.zdim)
        if self.zdim > 0:
            x = torch.cat([x, coords[..., 3:]], -1)
            assert x.shape[-1] == self.in_dim
        return x

    def positional_encoding_linear(self, coords):
        """Expand coordinates in the Fourier basis, i.e. cos(k*n/N), sin(k*n/N), n=0,...,N//2"""
        freqs = torch.arange(1, self.D2 + 1, dtype=torch.float, device=coords.device)
        freqs = freqs.view(*[1] * len(coords.shape), -1)  # 1 x 1 x D2
        coords = coords.unsqueeze(-1)  # B x 3 x 1
        k = coords[..., 0:3, :] * freqs  # B x 3 x D2
        s = torch.sin(k)  # B x 3 x D2
        c = torch.cos(k)  # B x 3 x D2
        x = torch.cat([s, c], -1)  # B x 3 x D
        x = x.view(*coords.shape[:-2], self.in_dim - self.zdim)  # B x in_dim-zdim
        if self.zdim > 0:
            x = torch.cat([x, coords[..., 3:, :].squeeze(-1)], -1)
            assert x.shape[-1] == self.in_dim
        return x

    def forward(self, coords: Tensor) -> Tensor:
        """Input should be coordinates from [-.5,.5]"""
        assert (coords[..., 0:3].abs() - 0.5 < 1e-4).all()
        return self.decoder(self.positional_encoding_geom(coords))

    def eval_volume(
        self,
        coords: Tensor,
        D: int,
        extent: float,
        norm: Norm,
        zval: Optional[np.ndarray] = None,
    ) -> Tensor:
        """
        Evaluate the model on a DxDxD volume

        Inputs:
            coords: lattice coords on the x-y plane (D^2 x 3)
            D: size of lattice
            extent: extent of lattice [-extent, extent]
            norm: data normalization
            zval: value of latent (zdim x 1)
        """
        # Note: extent should be 0.5 by default, except when a downsampled
        # volume is generated
        assert extent <= 0.5
        zdim = 0
        z = torch.tensor([])
        if zval is not None:
            zdim = len(zval)
            z = torch.tensor(zval, dtype=torch.float32, device=coords.device)

        vol_f = torch.zeros((D, D, D), dtype=torch.float32)
        assert not self.training
        # evaluate the volume by zslice to avoid memory overflows
        for i, dz in enumerate(
            np.linspace(-extent, extent, D, endpoint=True, dtype=np.float32)
        ):
            x = coords + torch.tensor([0, 0, dz], device=coords.device)
            if zval is not None:
                x = torch.cat((x, z.expand(x.shape[0], zdim)), dim=-1)
            with torch.no_grad():
                y = self.forward(x)
                y = y.view(D, D)
            vol_f[i] = y
        vol_f = vol_f * norm[1] + norm[0]
        vol = fft.ihtn_center(
            vol_f[0:-1, 0:-1, 0:-1]
        )  # remove last +k freq for inverse FFT
        return vol


class FTPositionalDecoder(Decoder):
    def __init__(
        self,
        in_dim: int,
        D: int,
        nlayers: int,
        hidden_dim: int,
        activation: Type,
        enc_type: str = "linear_lowf",
        enc_dim: Optional[int] = None,
        feat_sigma: Optional[float] = None,
    ):
        super(FTPositionalDecoder, self).__init__()
        assert in_dim >= 3
        self.zdim = in_dim - 3
        self.D = D
        self.D2 = D // 2
        self.DD = 2 * (D // 2)
        self.enc_type = enc_type
        self.enc_dim = self.D2 if enc_dim is None else enc_dim
        self.in_dim = 3 * (self.enc_dim) * 2 + self.zdim
        self.decoder = ResidLinearMLP(self.in_dim, nlayers, hidden_dim, 2, activation)

        if enc_type == "gaussian":
            # We construct 3 * self.enc_dim random vector frequences, to match the original positional encoding:
            # In the positional encoding we produce self.enc_dim features for each of the x,y,z dimensions,
            # whereas in gaussian encoding we produce self.enc_dim features each with random x,y,z components
            #
            # Each of the random feats is the sine/cosine of the dot product of the coordinates with a frequency
            # vector sampled from a gaussian with std of feat_sigma
            rand_freqs = (
                torch.randn((3 * self.enc_dim, 3), dtype=torch.float) * feat_sigma
            )
            # make rand_feats a parameter so it is saved in the checkpoint, but do not perform SGD on it
            self.rand_freqs = Parameter(rand_freqs, requires_grad=False)
        else:
            self.rand_feats = None

    def positional_encoding_geom(self, coords: Tensor) -> Tensor:
        """Expand coordinates in the Fourier basis with geometrically spaced wavelengths from 2/D to 2pi"""
        if self.enc_type == "gaussian":
            return self.random_fourier_encoding(coords)
        freqs = torch.arange(self.enc_dim, dtype=torch.float, device=coords.device)
        if self.enc_type == "geom_ft":
            freqs = (
                self.DD * np.pi * (2.0 / self.DD) ** (freqs / (self.enc_dim - 1))
            )  # option 1: 2/D to 1
        elif self.enc_type == "geom_full":
            freqs = (
                self.DD
                * np.pi
                * (1.0 / self.DD / np.pi) ** (freqs / (self.enc_dim - 1))
            )  # option 2: 2/D to 2pi
        elif self.enc_type == "geom_lowf":
            freqs = self.D2 * (1.0 / self.D2) ** (
                freqs / (self.enc_dim - 1)
            )  # option 3: 2/D*2pi to 2pi
        elif self.enc_type == "geom_nohighf":
            freqs = self.D2 * (2.0 * np.pi / self.D2) ** (
                freqs / (self.enc_dim - 1)
            )  # option 4: 2/D*2pi to 1
        elif self.enc_type == "linear_lowf":
            return self.positional_encoding_linear(coords)
        else:
            raise RuntimeError("Encoding type {} not recognized".format(self.enc_type))
        freqs = freqs.view(*[1] * len(coords.shape), -1)  # 1 x 1 x D2
        coords = coords.unsqueeze(-1)  # B x 3 x 1
        k = coords[..., 0:3, :] * freqs  # B x 3 x D2
        s = torch.sin(k)  # B x 3 x D2
        c = torch.cos(k)  # B x 3 x D2
        x = torch.cat([s, c], -1)  # B x 3 x D
        x = x.view(*coords.shape[:-2], self.in_dim - self.zdim)  # B x in_dim-zdim
        if self.zdim > 0:
            x = torch.cat([x, coords[..., 3:, :].squeeze(-1)], -1)
            assert x.shape[-1] == self.in_dim
        return x

    def random_fourier_encoding(self, coords):
        assert self.rand_freqs is not None
        # k = coords . rand_freqs
        # expand rand_freqs with singleton dimension along the batch dimensions
        # e.g. dim (1, ..., 1, n_rand_feats, 3)
        freqs = self.rand_freqs.view(*[1] * (len(coords.shape) - 1), -1, 3) * self.D2

        kxkykz = coords[..., None, 0:3] * freqs  # compute the x,y,z components of k
        k = kxkykz.sum(-1)  # compute k
        s = torch.sin(k)
        c = torch.cos(k)
        x = torch.cat([s, c], -1)
        x = x.view(*coords.shape[:-1], self.in_dim - self.zdim)
        if self.zdim > 0:
            x = torch.cat([x, coords[..., 3:]], -1)
            assert x.shape[-1] == self.in_dim
        return x

    def positional_encoding_linear(self, coords: Tensor) -> Tensor:
        """Expand coordinates in the Fourier basis, i.e. cos(k*n/N), sin(k*n/N), n=0,...,N//2"""
        freqs = torch.arange(1, self.D2 + 1, dtype=torch.float, device=coords.device)
        freqs = freqs.view(*[1] * len(coords.shape), -1)  # 1 x 1 x D2
        coords = coords.unsqueeze(-1)  # B x 3 x 1
        k = coords[..., 0:3, :] * freqs  # B x 3 x D2
        s = torch.sin(k)  # B x 3 x D2
        c = torch.cos(k)  # B x 3 x D2
        x = torch.cat([s, c], -1)  # B x 3 x D
        x = x.view(*coords.shape[:-2], self.in_dim - self.zdim)  # B x in_dim-zdim
        if self.zdim > 0:
            x = torch.cat([x, coords[..., 3:, :].squeeze(-1)], -1)
            assert x.shape[-1] == self.in_dim
        return x

    def forward(self, lattice: Tensor) -> Tensor:
        """
        Call forward on central slices only
            i.e. the middle pixel should be (0,0,0)

        lattice: B x N x 3+zdim
        """
        # if ignore_DC = False, then the size of the lattice will be odd (since it
        # includes the origin), so we need to evaluate one additional pixel
        c = lattice.shape[-2] // 2  # top half
        cc = c + 1 if lattice.shape[-2] % 2 == 1 else c  # include the origin
        assert abs(lattice[..., 0:3].mean()) < 1e-4, "{} != 0.0".format(
            lattice[..., 0:3].mean()
        )
        image = torch.empty(lattice.shape[:-1], device=lattice.device)
        top_half = self.decode(lattice[..., 0:cc, :])
        image[..., 0:cc] = top_half[..., 0] - top_half[..., 1]
        # the bottom half of the image is the complex conjugate of the top half
        image[..., cc:] = (top_half[..., 0] + top_half[..., 1])[
            ..., np.arange(c - 1, -1, -1)
        ]
        return image

    def decode(self, lattice: Tensor):
        """Return FT transform"""
        assert (lattice[..., 0:3].abs() - 0.5 < 1e-4).all()
        # convention: only evalute the -z points
        w = lattice[..., 2] > 0.0
        new_lattice = lattice.clone()
        # negate lattice coordinates where z > 0
        new_lattice[..., 0:3][w] *= -1
        result = self.decoder(self.positional_encoding_geom(new_lattice))
        # replace with complex conjugate to get correct values for original lattice positions
        result[..., 1][w] *= -1
        return result

    def eval_volume(
        self,
        coords: Tensor,
        D: int,
        extent: float,
        norm: Norm,
        zval: Optional[np.ndarray] = None,
    ) -> Tensor:
        """
        Evaluate the model on a DxDxD volume

        Inputs:
            coords: lattice coords on the x-y plane (D^2 x 3)
            D: size of lattice
            extent: extent of lattice [-extent, extent]
            norm: data normalization
            zval: value of latent (zdim x 1)
        """
        assert extent <= 0.5
        zdim = 0
        z = torch.tensor([])
        if zval is not None:
            zdim = len(zval)
            z = torch.tensor(zval, dtype=torch.float32, device=coords.device)

        vol_f = torch.zeros((D, D, D), dtype=torch.float32)
        assert not self.training
        # evaluate the volume by zslice to avoid memory overflows
        for i, dz in enumerate(
            np.linspace(-extent, extent, D, endpoint=True, dtype=np.float32)
        ):
            x = coords + torch.tensor([0, 0, dz], device=coords.device)
            keep = x.pow(2).sum(dim=1) <= extent**2
            x = x[keep]
            if zval is not None:
                x = torch.cat((x, z.expand(x.shape[0], zdim)), dim=-1)
            with torch.no_grad():
                if dz == 0.0:
                    y = self.forward(x)
                else:
                    y = self.decode(x)
                    y = y[..., 0] - y[..., 1]
                slice_ = torch.zeros(D**2, device="cpu")
                slice_[keep] = y.cpu()
                slice_ = slice_.view(D, D)
            vol_f[i] = slice_
        vol_f = vol_f * norm[1] + norm[0]
        vol = fft.ihtn_center(
            vol_f[:-1, :-1, :-1]
        )  # remove last +k freq for inverse FFT
        return vol


class FTSliceDecoder(Decoder):
    """
    Evaluate a central slice out of a 3D FT of a model, returns representation in
    Hartley reciprocal space

    Exploits the symmetry of the FT where F*(x,y) = F(-x,-y) and only
    evaluates half of the lattice. The decoder is f(x,y,z) => real, imag
    """

    def __init__(self, in_dim: int, D: int, nlayers: int, hidden_dim: int, activation):
        """D: image width or height"""
        super(FTSliceDecoder, self).__init__()
        self.decoder = ResidLinearMLP(in_dim, nlayers, hidden_dim, 2, activation)
        D2 = int(D / 2)

        # various pixel indices to keep track of for forward_even
        self.center = D2 * D + D2
        self.extra = np.arange(
            (D2 + 1) * D, D**2, D
        )  # bottom-left column without conjugate pair
        # evalute the top half of the image up through the center pixel
        # and extra bottom-left column (todo: just evaluate a D-1 x D-1 image so
        # we don't have to worry about this)
        self.all_eval = np.concatenate((np.arange(self.center + 1), self.extra))

        # pixel indices for the top half of the image up to (but not incl)
        # the center pixel and excluding the top row and left-most column
        i, j = np.meshgrid(np.arange(1, D), np.arange(1, D2 + 1))
        self.top = (j * D + i).ravel()[:-D2]

        # pixel indices for bottom half of the image after the center pixel
        # excluding left-most column and given in reverse order
        i, j = np.meshgrid(np.arange(1, D), np.arange(D2, D))
        self.bottom_rev = (j * D + i).ravel()[D2:][::-1].copy()

        self.D = D
        self.D2 = D2

    def forward(self, lattice):
        """
        Call forward on central slices only
            i.e. the middle pixel should be (0,0,0)

        lattice: B x N x 3+zdim
        """
        assert lattice.shape[-2] % 2 == 1
        c = lattice.shape[-2] // 2  # center pixel
        assert lattice[..., c, 0:3].sum() == 0.0, "{} != 0.0".format(
            lattice[..., c, 0:3].sum()
        )
        assert abs(lattice[..., 0:3].mean()) < 1e-4, "{} != 0.0".format(
            lattice[..., 0:3].mean()
        )
        image = torch.empty(lattice.shape[:-1], device=lattice.device)
        top_half = self.decode(lattice[..., 0 : c + 1, :])
        image[..., 0 : c + 1] = top_half[..., 0] - top_half[..., 1]
        # the bottom half of the image is the complex conjugate of the top half
        image[..., c + 1 :] = (top_half[..., 0] + top_half[..., 1])[
            ..., np.arange(c - 1, -1, -1)
        ]
        return image

    def forward_even(self, lattice):
        """Extra bookkeeping with extra row/column for an even sized DFT"""
        image = torch.empty(lattice.shape[:-1], device=lattice.device)
        top_half = self.decode(lattice[..., self.all_eval, :])
        image[..., self.all_eval] = top_half[..., 0] - top_half[..., 1]
        # the bottom half of the image is the complex conjugate of the top half
        image[..., self.bottom_rev] = (
            top_half[..., self.top, 0] + top_half[..., self.top, 1]
        )
        return image

    def decode(self, lattice):
        """Return FT transform"""
        # convention: only evalute the -z points
        w = lattice[..., 2] > 0.0
        new_lattice = lattice.clone()
        # negate lattice coordinates where z > 0

        new_lattice[..., 0:3][w] *= -1
        result = self.decoder(new_lattice)
        # replace with complex conjugate to get correct values for original lattice positions
        result[..., 1][w] *= -1
        return result

    def eval_volume(
        self,
        coords: Tensor,
        D: int,
        extent: float,
        norm: Norm,
        zval: Optional[np.ndarray] = None,
    ) -> Tensor:
        """
        Evaluate the model on a DxDxD volume

        Inputs:
            coords: lattice coords on the x-y plane (D^2 x 3)
            D: size of lattice
            extent: extent of lattice [-extent, extent]
            norm: data normalization
            zval: value of latent (zdim x 1)
        """
        if zval is not None:
            zdim = len(zval)
            z = torch.zeros(D**2, zdim, dtype=torch.float32)
            z += torch.tensor(zval, dtype=torch.float32, device=coords.device)
        else:
            z = None

        vol_f = torch.zeros((D, D, D), dtype=torch.float32)
        assert not self.training
        # evaluate the volume by zslice to avoid memory overflows
        for i, dz in enumerate(
            np.linspace(-extent, extent, D, endpoint=True, dtype=np.float32)
        ):
            x = coords + torch.tensor([0, 0, dz], device=coords.device)
            if zval is not None:
                assert z is not None
                x = torch.cat((x, z), dim=-1)
            with torch.no_grad():
                y = self.decode(x)
                y = y[..., 0] - y[..., 1]
                y = y.view(D, D).cpu()
            vol_f[i] = y
        vol_f = vol_f * norm[1] + norm[0]
        vol_f = utils.zero_sphere(vol_f)
        vol = fft.ihtn_center(
            vol_f[:-1, :-1, :-1]
        )  # remove last +k freq for inverse FFT
        return vol


def get_decoder(
    in_dim: int,
    D: int,
    layers: int,
    dim: int,
    domain: str,
    enc_type: str,
    enc_dim: Optional[int] = None,
    activation: Type = nn.ReLU,
    feat_sigma: Optional[float] = None,
) -> Decoder:
    if enc_type == "none":
        if domain == "hartley":
            model = ResidLinearMLP(in_dim, layers, dim, 1, activation)
        else:
            model = FTSliceDecoder(in_dim, D, layers, dim, activation)
    else:
        model_t = PositionalDecoder if domain == "hartley" else FTPositionalDecoder
        model = model_t(
            in_dim,
            D,
            layers,
            dim,
            activation,
            enc_type=enc_type,
            enc_dim=enc_dim,
            feat_sigma=feat_sigma,
        )
    return model


class VAE(nn.Module):
    def __init__(
        self,
        lattice,
        qlayers: int,
        qdim: int,
        players: int,
        pdim: int,
        encode_mode: str = "mlp",
        no_trans: bool = False,
        enc_mask: Optional[Tensor] = None,
    ):
        super(VAE, self).__init__()
        self.lattice = lattice
        self.D = lattice.D
        if enc_mask is not None:
            self.in_dim = (
                lattice.D * lattice.D if enc_mask is None else int(enc_mask.sum())
            )
        self.enc_mask = enc_mask
        assert qlayers > 2
        if encode_mode == "conv":
            self.encoder = ConvEncoder(qdim, qdim)
        elif encode_mode == "resid":
            self.encoder = ResidLinearMLP(
                self.in_dim,
                qlayers - 2,  # -2 bc we add 2 more layers in the homeomorphic encoer
                qdim,  # hidden_dim
                qdim,  # out_dim
                nn.ReLU,
            )  # in_dim -> hidden_dim
        elif encode_mode == "mlp":
            self.encoder = MLP(
                self.in_dim, qlayers - 2, qdim, qdim, nn.ReLU  # hidden_dim  # out_dim
            )  # in_dim -> hidden_dim
        else:
            raise RuntimeError("Encoder mode {} not recognized".format(encode_mode))
        # predict rotation and translation in two completely separate NNs
        # self.so3_encoder = SO3reparameterize(qdim) # hidden_dim -> SO(3) latent variable
        # self.trans_encoder = ResidLinearMLP(nx*ny, 5, qdim, 4, nn.ReLU)

        # or predict rotation/translations from intermediate encoding
        self.so3_encoder = SO3reparameterize(
            qdim, 1, qdim
        )  # hidden_dim -> SO(3) latent variable
        self.trans_encoder = ResidLinearMLP(qdim, 1, qdim, 4, nn.ReLU)

        self.decoder = FTSliceDecoder(3, self.D, players, pdim, nn.ReLU)
        self.no_trans = no_trans

    def reparameterize(self, mu: Tensor, logvar: Tensor) -> Tensor:
        if not self.training:
            return mu
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return eps * std + mu

    def encode(self, img) -> Tuple[Tensor, Tensor, Optional[Tensor], Optional[Tensor]]:
        """img: BxDxD"""
        img = img.view(img.size(0), -1)
        if self.enc_mask is not None:
            img = img[:, self.enc_mask]
        enc = nn.ReLU()(self.encoder(img))
        z_mu, z_std = self.so3_encoder(enc)
        if self.no_trans:
            tmu, tlogvar = (None, None)
        else:
            z = self.trans_encoder(enc)
            tmu, tlogvar = z[:, :2], z[:, 2:]
        return z_mu, z_std, tmu, tlogvar

    def eval_volume(self, norm) -> Tensor:
        return self.decoder.eval_volume(
            self.lattice.coords, self.D, self.lattice.extent, norm
        )

    def decode(self, rot):
        # transform lattice by rot.T
        x = self.lattice.coords @ rot  # R.T*x
        y_hat = self.decoder(x)
        y_hat = y_hat.view(-1, self.D, self.D)
        return y_hat

    def forward(self, img: Tensor):
        z_mu, z_std, tmu, tlogvar = self.encode(img)
        rot, w_eps = self.so3_encoder.sampleSO3(z_mu, z_std)
        # transform lattice by rot and predict image
        y_hat = self.decode(rot)
        if not self.no_trans:
            # translate image by t
            assert tmu is not None and tlogvar is not None
            B = img.size(0)
            t = self.reparameterize(tmu, tlogvar)
            t = t.unsqueeze(1)  # B x 1 x 2
            img = self.lattice.translate_ht(img.view(B, -1), t)
            img = img.view(B, self.D, self.D)
        return y_hat, img, z_mu, z_std, w_eps, tmu, tlogvar


class TiltVAE(nn.Module):
    def __init__(
        self, lattice, tilt, qlayers, qdim, players, pdim, no_trans=False, enc_mask=None
    ):
        super(TiltVAE, self).__init__()
        self.lattice = lattice
        self.D = lattice.D
        self.in_dim = lattice.D * lattice.D if enc_mask is None else enc_mask.sum()
        self.enc_mask = enc_mask
        assert qlayers > 3
        self.encoder = ResidLinearMLP(self.in_dim, qlayers - 3, qdim, qdim, nn.ReLU)
        self.so3_encoder = SO3reparameterize(
            2 * qdim, 3, qdim
        )  # hidden_dim -> SO(3) latent variable
        self.trans_encoder = ResidLinearMLP(2 * qdim, 2, qdim, 4, nn.ReLU)
        self.decoder = FTSliceDecoder(3, self.D, players, pdim, nn.ReLU)
        assert tilt.shape == (3, 3), "Rotation matrix input required"
        self.tilt = torch.tensor(tilt)
        self.no_trans = no_trans

    def reparameterize(self, mu, logvar):
        if not self.training:
            return mu
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return eps * std + mu

    def eval_volume(self, norm) -> Tensor:
        return self.decoder.eval_volume(
            self.lattice.coords, self.D, self.lattice.extent, norm
        )

    def encode(self, img, img_tilt):
        img = img.view(img.size(0), -1)
        img_tilt = img_tilt.view(img_tilt.size(0), -1)
        if self.enc_mask is not None:
            img = img[:, self.enc_mask]
            img_tilt = img_tilt[:, self.enc_mask]
        enc1 = self.encoder(img)
        enc2 = self.encoder(img_tilt)
        enc = torch.cat((enc1, enc2), -1)  # then nn.ReLU?
        z_mu, z_std = self.so3_encoder(enc)
        rot, w_eps = self.so3_encoder.sampleSO3(z_mu, z_std)
        if self.no_trans:
            tmu, tlogvar, t = (None, None, None)
        else:
            z = self.trans_encoder(enc)
            tmu, tlogvar = z[:, :2], z[:, 2:]
            t = self.reparameterize(tmu, tlogvar)
        return z_mu, z_std, w_eps, rot, tmu, tlogvar, t

    def forward(self, img, img_tilt):
        B = img.size(0)
        z_mu, z_std, w_eps, rot, tmu, tlogvar, t = self.encode(img, img_tilt)
        if not self.no_trans:
            assert t is not None
            t = t.unsqueeze(1)  # B x 1 x 2
            img = self.lattice.translate_ht(img.view(B, -1), -t)
            img_tilt = self.lattice.translate_ht(img_tilt.view(B, -1), -t)
            img = img.view(B, self.D, self.D)
            img_tilt = img_tilt.view(B, self.D, self.D)

        # rotate lattice by rot.T
        x = self.lattice.coords @ rot  # R.T*x
        y_hat = self.decoder(x)
        y_hat = y_hat.view(-1, self.D, self.D)

        # tilt series pair
        x = self.lattice.coords @ self.tilt @ rot
        y_hat2 = self.decoder(x)
        y_hat2 = y_hat2.view(-1, self.D, self.D)
        return y_hat, y_hat2, img, img_tilt, z_mu, z_std, w_eps, tmu, tlogvar


# fixme: this is half-deprecated (not used in TiltVAE, but still used in tilt BNB)
class TiltEncoder(nn.Module):
    def __init__(
        self,
        in_dim,
        nlayers,
        hidden_dim,
        out_dim,
        ntilts,
        nlayers2,
        hidden_dim2,
        out_dim2,
        activation,
    ):
        super(TiltEncoder, self).__init__()
        self.encoder1 = ResidLinearMLP(in_dim, nlayers, hidden_dim, out_dim, activation)
        self.encoder2 = ResidLinearMLP(
            out_dim * ntilts, nlayers2, hidden_dim2, out_dim2, activation
        )
        self.in_dim = in_dim
        self.in_dim2 = out_dim * ntilts

    def forward(self, x):
        x = self.encoder1(x)
        z = self.encoder2(x.view(-1, self.in_dim2))
        return z


class ResidLinearMLP(Decoder):
    def __init__(
        self,
        in_dim: int,
        nlayers: int,
        hidden_dim: int,
        out_dim: int,
        activation: Type,
    ):
        super(ResidLinearMLP, self).__init__()
        layers = [
            ResidLinear(in_dim, hidden_dim)
            if in_dim == hidden_dim
            else MyLinear(in_dim, hidden_dim),
            activation(),
        ]
        for n in range(nlayers):
            layers.append(ResidLinear(hidden_dim, hidden_dim))
            layers.append(activation())
        layers.append(
            ResidLinear(hidden_dim, out_dim)
            if out_dim == hidden_dim
            else MyLinear(hidden_dim, out_dim)
        )
        self.main = nn.Sequential(*layers)

    def forward(self, x):
        flat = x.view(-1, x.shape[-1])
        ret_flat = self.main(flat)
        ret = ret_flat.view(*x.shape[:-1], ret_flat.shape[-1])
        return ret

    def eval_volume(
        self, coords: Tensor, D: int, extent: float, norm: Norm, zval=None
    ) -> Tensor:
        """
        Evaluate the model on a DxDxD volume

        Inputs:
            coords: lattice coords on the x-y plane (D^2 x 3)
            D: size of lattice
            extent: extent of lattice [-extent, extent]
            norm: data normalization
            zval: value of latent (zdim x 1)
        """
        # Note: extent should be 0.5 by default, except when a downsampled
        # volume is generated
        if zval is not None:
            zdim = len(zval)
            z = torch.zeros(D**2, zdim, dtype=torch.float32, device=coords.device)
            z += torch.tensor(zval, dtype=torch.float32, device=coords.device)

        vol_f = torch.zeros((D, D, D), dtype=torch.float32)
        assert not self.training
        # evaluate the volume by zslice to avoid memory overflows
        for i, dz in enumerate(
            np.linspace(-extent, extent, D, endpoint=True, dtype=np.float32)
        ):
            x = coords + torch.tensor([0, 0, dz], device=coords.device)
            if zval is not None:
                x = torch.cat((x, zval), dim=-1)
            with torch.no_grad():
                y = self.forward(x)
                y = y.view(D, D).cpu()
            vol_f[i] = y
        vol_f = vol_f * norm[1] + norm[0]
        vol = fft.ihtn_center(
            vol_f[0:-1, 0:-1, 0:-1]
        )  # remove last +k freq for inverse FFT
        return vol


def half_linear(input, weight, bias):
    # print('half', input.shape, weight.shape)
    return F.linear(input, weight.half(), bias.half())


def single_linear(input, weight, bias):
    # print('single', input.shape, weight.shape)
    # assert input.shape[0] < 10000

    return F.linear(input, weight, bias)


class MyLinear(nn.Linear):
    def forward(self, input):
        if input.dtype == torch.half:
            return half_linear(
                input, self.weight, self.bias
            )  # F.linear(input, self.weight.half(), self.bias.half())
        else:
            return single_linear(
                input, self.weight, self.bias
            )  # F.linear(input, self.weight, self.bias)


class ResidLinear(nn.Module):
    def __init__(self, nin, nout):
        super(ResidLinear, self).__init__()
        self.linear = MyLinear(nin, nout)
        # self.linear = nn.utils.weight_norm(MyLinear(nin, nout))

    def forward(self, x):
        z = self.linear(x) + x
        return z


class MLP(nn.Module):
    def __init__(
        self,
        in_dim: int,
        nlayers: int,
        hidden_dim: int,
        out_dim: int,
        activation: Type,
    ):
        super(MLP, self).__init__()
        layers = [MyLinear(in_dim, hidden_dim), activation()]
        for n in range(nlayers):
            layers.append(MyLinear(hidden_dim, hidden_dim))
            layers.append(activation())
        layers.append(MyLinear(hidden_dim, out_dim))
        self.main = nn.Sequential(*layers)

    def forward(self, x):
        return self.main(x)


# Adapted from soumith DCGAN
class ConvEncoder(nn.Module):
    def __init__(self, hidden_dim, out_dim):
        super(ConvEncoder, self).__init__()
        ndf = hidden_dim
        self.main = nn.Sequential(
            # input is 1 x 64 x 64
            nn.Conv2d(1, ndf, 4, 2, 1, bias=False),
            nn.LeakyReLU(0.2, inplace=True),
            # state size. (ndf) x 32 x 32
            nn.Conv2d(ndf, ndf * 2, 4, 2, 1, bias=False),
            # nn.BatchNorm2d(ndf * 2),
            nn.LeakyReLU(0.2, inplace=True),
            # state size. (ndf*2) x 16 x 16
            nn.Conv2d(ndf * 2, ndf * 4, 4, 2, 1, bias=False),
            # nn.BatchNorm2d(ndf * 4),
            nn.LeakyReLU(0.2, inplace=True),
            # state size. (ndf*4) x 8 x 8
            nn.Conv2d(ndf * 4, ndf * 8, 4, 2, 1, bias=False),
            # nn.BatchNorm2d(ndf * 8),
            nn.LeakyReLU(0.2, inplace=True),
            # state size. (ndf*8) x 4 x 4
            nn.Conv2d(ndf * 8, out_dim, 4, 1, 0, bias=False),
            # state size. out_dims x 1 x 1
        )

    def forward(self, x):
        x = x.view(-1, 1, 64, 64)
        x = self.main(x)
        return x.view(x.size(0), -1)  # flatten


class SO3reparameterize(nn.Module):
    """Reparameterize R^N encoder output to SO(3) latent variable"""

    def __init__(self, input_dims, nlayers: int, hidden_dim: int):
        super().__init__()
        if nlayers is not None:
            self.main = ResidLinearMLP(input_dims, nlayers, hidden_dim, 9, nn.ReLU)
        else:
            self.main = MyLinear(input_dims, 9)

        # start with big outputs
        # self.s2s2map.weight.data.uniform_(-5,5)
        # self.s2s2map.bias.data.uniform_(-5,5)

    def sampleSO3(
        self, z_mu: torch.Tensor, z_std: torch.Tensor
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Reparameterize SO(3) latent variable
        # z represents mean on S2xS2 and variance on so3, which enocdes a Gaussian distribution on SO3
        # See section 2.5 of http://ethaneade.com/lie.pdf
        """
        # resampling trick
        if not self.training:
            return z_mu, z_std
        eps = torch.randn_like(z_std)
        w_eps = eps * z_std
        rot_eps = lie_tools.expmap(w_eps)
        # z_mu = lie_tools.quaternions_to_SO3(z_mu)
        rot_sampled = z_mu @ rot_eps
        return rot_sampled, w_eps

    def forward(self, x) -> Tuple[torch.Tensor, torch.Tensor]:
        z = self.main(x)
        z1 = z[:, :3].double()
        z2 = z[:, 3:6].double()
        z_mu = lie_tools.s2s2_to_SO3(z1, z2).float()
        logvar = z[:, 6:]
        z_std = torch.exp(0.5 * logvar)  # or could do softplus
        return z_mu, z_std
