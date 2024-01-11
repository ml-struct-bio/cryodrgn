"""Pytorch models implenting neural networks."""

import numpy as np
import torch
from torch import Tensor
from torch import nn
import torch.nn.functional as F
from torch.nn.parameter import Parameter
from typing import Optional, Tuple, Type, Sequence, Any

from cryodrgn import fft
from cryodrgn import lie_tools
from cryodrgn import utils
import cryodrgn.config

Norm = Sequence[Any]  # mean, std


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


def load_decoder(config, weights=None, device=None) -> Decoder:
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
            ResidualLinear(in_dim, hidden_dim)
            if in_dim == hidden_dim
            else MyLinear(in_dim, hidden_dim),
            activation(),
        ]
        for n in range(nlayers):
            layers.append(ResidualLinear(hidden_dim, hidden_dim))
            layers.append(activation())
        layers.append(
            ResidualLinear(hidden_dim, out_dim)
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


class ResidualLinear(nn.Module):
    def __init__(self, nin, nout):
        super(ResidualLinear, self).__init__()
        self.linear = MyLinear(nin, nout)
        # self.linear = nn.utils.weight_norm(MyLinear(nin, nout))

    def forward(self, x):
        z = self.linear(x) + x
        return z


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


class MyDataParallel(nn.DataParallel):
    def __getattr__(self, name):
        try:
            return super().__getattr__(name)
        except AttributeError:
            return getattr(self.module, name)
