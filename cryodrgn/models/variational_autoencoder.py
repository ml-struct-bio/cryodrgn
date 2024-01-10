"""Pytorch models implementing variational autoencoders."""

import torch
from torch import Tensor
import torch.nn as nn
from torch.nn.parallel import DataParallel
from typing import Optional, Tuple, Type
from cryodrgn.models.neural_nets import (
    DataParallelDecoder,
    half_linear,
    single_linear,
    ResidLinearMLP,
    TiltEncoder,
    get_decoder,
)
from cryodrgn import lie_tools
from cryodrgn.lattice import Lattice


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


class ResidualLinear(nn.Module):
    def __init__(self, nin, nout):
        super(ResidualLinear, self).__init__()
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
