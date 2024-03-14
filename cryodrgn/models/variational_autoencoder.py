"""Pytorch models implementing variational autoencoders."""

import torch
import torch.nn as nn
from torch.nn.parallel import DataParallel
from typing import Optional, Tuple, Type
from cryodrgn.models.neural_nets import (
    DataParallelDecoder,
    MyLinear,
    ResidLinearMLP,
    TiltEncoder,
    get_decoder,
)
from cryodrgn.lattice import Lattice

import logging

logger = logging.getLogger(__name__)


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
        z_dim: int = 1,
        encode_mode: str = "resid",
        enc_mask=None,
        enc_type="linear_lowf",
        enc_dim=None,
        domain="fourier",
        activation=nn.ReLU,
        feat_sigma: Optional[float] = None,
        tilt_params=None,
    ):
        super(HetOnlyVAE, self).__init__()
        self.lattice = lattice
        self.z_dim = z_dim
        self.in_dim = in_dim
        self.enc_mask = enc_mask
        if encode_mode == "conv":
            self.encoder = ConvEncoder(qdim, z_dim * 2)
        elif encode_mode == "resid":
            self.encoder = ResidLinearMLP(
                in_dim,
                qlayers,
                qdim,
                z_dim * 2,
                activation,  # nlayers  # hidden_dim  # out_dim
            )
        elif encode_mode == "mlp":
            self.encoder = MLP(
                in_dim, qlayers, qdim, z_dim * 2, activation  # hidden_dim  # out_dim
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
                z_dim * 2,  # outdim
                activation,
            )
        else:
            raise RuntimeError("Encoder mode {} not recognized".format(encode_mode))
        self.encode_mode = encode_mode
        self.decoder = get_decoder(
            3 + z_dim,
            lattice.D,
            players,
            pdim,
            domain,
            enc_type,
            enc_dim,
            activation,
            feat_sigma,
        )

    def reparameterize(self, mu, logvar):
        if not self.training:
            return mu
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return eps * std + mu

    def encode(self, *img) -> Tuple[torch.Tensor, torch.Tensor]:
        img = (x.view(x.shape[0], -1) for x in img)
        if self.enc_mask is not None:
            img = (x[:, self.enc_mask] for x in img)
        z = self.encoder(*img)
        return z[:, : self.z_dim], z[:, self.z_dim :]

    def cat_z(self, coords, z) -> torch.Tensor:
        """
        coords: B x ... x 3
        z: B x z_dim
        """
        assert coords.size(0) == z.size(0), (coords.shape, z.shape)
        z = z.view(z.size(0), *([1] * (coords.ndimension() - 2)), self.z_dim)
        z = torch.cat((coords, z.expand(*coords.shape[:-1], self.z_dim)), dim=-1)
        return z

    def decode(self, coords, z=None) -> torch.Tensor:
        """
        coords: BxNx3 image coordinates
        z: Bxz_dim latent coordinate
        """
        decoder = self.decoder
        assert isinstance(decoder, nn.Module)
        retval = decoder(self.cat_z(coords, z) if z is not None else coords)
        return retval

    # Need forward func for DataParallel -- TODO: refactor
    def forward(self, *args, **kwargs):
        return self.decode(*args, **kwargs)

    def eval_volume(
        self,
        coords=None,
        resolution=None,
        extent=None,
        norm=(0.0, 1.0),
        zval=None,
        **vol_args,
    ) -> torch.Tensor:
        """
        Evaluate the model on a DxDxD volume

        Inputs:
            coords: lattice coords on the x-y plane (D^2 x 3)
            D: size of lattice
            extent: extent of lattice [-extent, extent]
            norm: data normalization
            zval: value of latent (z_dim x 1)
        """

        if not hasattr(self.decoder, "eval_volume"):
            raise NotImplementedError

        if coords is None:
            coords = self.lattice.coords
        if resolution is None:
            resolution = self.lattice.D
        if extent is None:
            extent = self.lattice.extent

        # TODO: kludge because VAE and drgnai models have different loading APIs
        for k, v in vol_args.items():
            logger.info(f"ignoring argument {k}={v} in VAE volume generation")

        return self.decoder.eval_volume(coords, resolution, extent, norm, zval)


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
