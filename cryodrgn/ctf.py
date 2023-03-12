from typing import Optional
import numpy as np
import seaborn as sns
import torch
import logging
from cryodrgn import utils

logger = logging.getLogger(__name__)


def compute_ctf(
    freqs: torch.Tensor,
    dfu: torch.Tensor,
    dfv: torch.Tensor,
    dfang: torch.Tensor,
    volt: torch.Tensor,
    cs: torch.Tensor,
    w: torch.Tensor,
    phase_shift: torch.Tensor = torch.Tensor([0]),
    bfactor: Optional[torch.Tensor] = None,
) -> torch.Tensor:
    """
    Compute the 2D CTF

    Input:
        freqs (np.ndarray) Nx2 or BxNx2 tensor of 2D spatial frequencies
        dfu (float or Bx1 tensor): DefocusU (Angstrom)
        dfv (float or Bx1 tensor): DefocusV (Angstrom)
        dfang (float or Bx1 tensor): DefocusAngle (degrees)
        volt (float or Bx1 tensor): accelerating voltage (kV)
        cs (float or Bx1 tensor): spherical aberration (mm)
        w (float or Bx1 tensor): amplitude contrast ratio
        phase_shift (float or Bx1 tensor): degrees
        bfactor (float or Bx1 tensor): envelope fcn B-factor (Angstrom^2)
    """
    assert freqs.shape[-1] == 2
    # convert units
    volt = volt * 1000
    cs = cs * 10**7
    dfang = dfang * np.pi / 180
    phase_shift = phase_shift * np.pi / 180

    # lam = sqrt(h^2/(2*m*e*Vr)); Vr = V + (e/(2*m*c^2))*V^2
    lam = 12.2639 / (volt + 0.97845e-6 * volt**2) ** 0.5
    x = freqs[..., 0]
    y = freqs[..., 1]
    ang = torch.atan2(y, x)
    s2 = x**2 + y**2
    df = 0.5 * (dfu + dfv + (dfu - dfv) * torch.cos(2 * (ang - dfang)))
    gamma = (
        2 * np.pi * (-0.5 * df * lam * s2 + 0.25 * cs * lam**3 * s2**2)
        - phase_shift
    )
    ctf = (1 - w**2) ** 0.5 * torch.sin(gamma) - w * torch.cos(gamma)
    if bfactor is not None:
        ctf *= torch.exp(-bfactor / 4 * s2)
    return ctf


def compute_ctf_np(
    freqs: np.ndarray,
    dfu: float,
    dfv: float,
    dfang: float,
    volt: float,
    cs: float,
    w: float,
    phase_shift: float = 0,
    bfactor: Optional[float] = None,
) -> np.ndarray:
    """
    Compute the 2D CTF

    Input:
        freqs (np.ndarray) Nx2 array of 2D spatial frequencies
        dfu (float): DefocusU (Angstrom)
        dfv (float): DefocusV (Angstrom)
        dfang (float): DefocusAngle (degrees)
        volt (float): accelerating voltage (kV)
        cs (float): spherical aberration (mm)
        w (float): amplitude contrast ratio
        phase_shift (float): degrees
        bfactor (float): envelope fcn B-factor (Angstrom^2)
    """
    # convert units
    volt = volt * 1000
    cs = cs * 10**7
    dfang = dfang * np.pi / 180
    phase_shift = phase_shift * np.pi / 180

    # lam = sqrt(h^2/(2*m*e*Vr)); Vr = V + (e/(2*m*c^2))*V^2
    lam = 12.2639 / np.sqrt(volt + 0.97845e-6 * volt**2)
    x = freqs[:, 0]
    y = freqs[:, 1]
    ang = np.arctan2(y, x)
    s2 = x**2 + y**2
    df = 0.5 * (dfu + dfv + (dfu - dfv) * np.cos(2 * (ang - dfang)))
    gamma = (
        2 * np.pi * (-0.5 * df * lam * s2 + 0.25 * cs * lam**3 * s2**2)
        - phase_shift
    )
    ctf = np.sqrt(1 - w**2) * np.sin(gamma) - w * np.cos(gamma)
    if bfactor is not None:
        ctf *= np.exp(-bfactor / 4 * s2)
    return np.require(ctf, dtype=freqs.dtype)


def print_ctf_params(params: np.ndarray) -> None:
    assert len(params) == 9
    logger.info("Image size (pix)  : {}".format(int(params[0])))
    logger.info("A/pix             : {}".format(params[1]))
    logger.info("DefocusU (A)      : {}".format(params[2]))
    logger.info("DefocusV (A)      : {}".format(params[3]))
    logger.info("Dfang (deg)       : {}".format(params[4]))
    logger.info("voltage (kV)      : {}".format(params[5]))
    logger.info("cs (mm)           : {}".format(params[6]))
    logger.info("w                 : {}".format(params[7]))
    logger.info("Phase shift (deg) : {}".format(params[8]))


def plot_ctf(D: int, Apix: float, ctf_params: np.ndarray) -> None:
    assert len(ctf_params) == 7

    freqs = (
        np.stack(
            np.meshgrid(
                np.linspace(-0.5, 0.5, D, endpoint=False),
                np.linspace(-0.5, 0.5, D, endpoint=False),
            ),
            -1,
        )
        / Apix
    )
    freqs = freqs.reshape(-1, 2)
    c = compute_ctf_np(freqs, *ctf_params)
    sns.heatmap(c.reshape(D, D))


def load_ctf_for_training(D: int, ctf_params_pkl: str) -> np.ndarray:
    assert D % 2 == 0
    ctf_params = utils.load_pkl(ctf_params_pkl)
    assert ctf_params.shape[1] == 9
    # Replace original image size with current dimensions
    Apix = ctf_params[0, 0] * ctf_params[0, 1] / D
    ctf_params[:, 0] = D
    ctf_params[:, 1] = Apix
    print_ctf_params(ctf_params[0])
    # Slice out the first column (D)
    return ctf_params[:, 1:]
