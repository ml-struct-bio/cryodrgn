"""Utility functions shared between various cryoDRGN operations and commands."""

from collections.abc import Hashable
import functools
import os
import subprocess
import pickle
import yaml
import logging
from typing import Tuple
import numpy as np
import torch

logger = logging.getLogger(__name__)


def meshgrid_2d(lo, hi, n, endpoint=False):
    """
    Torch-compatible implementation of:
    np.meshgrid(
            np.linspace(-0.5, 0.5, D, endpoint=endpoint),
            np.linspace(-0.5, 0.5, D, endpoint=endpoint),
        )
    Torch doesn't support the 'endpoint' argument (always assumed True)
    and the behavior of torch.meshgrid is different unless the 'indexing' argument is supplied.
    """
    if endpoint:
        values = torch.linspace(lo, hi, n)
    else:
        values = torch.linspace(lo, hi, n + 1)[:-1]

    return torch.meshgrid(values, values, indexing="xy")


class memoized(object):
    """Decorator. Caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned
    (not reevaluated).
    """

    def __init__(self, func):
        self.func = func
        self.cache = {}

    def __call__(self, *args):
        if not isinstance(args, Hashable):
            # uncacheable. a list, for instance.
            # better to not cache than blow up.
            return self.func(*args)
        if args in self.cache:
            return self.cache[args]
        else:
            value = self.func(*args)
            self.cache[args] = value
            return value

    def __repr__(self):
        """Return the function's docstring."""
        return self.func.__doc__

    def __get__(self, obj, objtype):
        """Support instance methods."""
        return functools.partial(self.__call__, obj)


def find_latest_output(outdir: str, outlbl: str = "weights") -> str:
    """Get the output file corresponding to the latest saved training epoch."""

    return os.path.join(
        outdir,
        sorted(
            [
                f
                for f in os.listdir(outdir)
                if os.path.isfile(os.path.join(outdir, f)) and f.startswith(outlbl)
            ],
            key=lambda x: int(x.split(".")[-2]) if x.split(".")[-2].isnumeric() else -1,
        )[-1],
    )


def load_pkl(pkl: str):
    with open(pkl, "rb") as f:
        x = pickle.load(f)
    return x


def save_pkl(data, out_pkl: str, mode: str = "wb") -> None:
    if mode == "wb" and os.path.exists(out_pkl):
        logger.warning(f"Warning: {out_pkl} already exists. Overwriting.")
    with open(out_pkl, mode) as f:
        pickle.dump(data, f)  # type: ignore


def load_yaml(yamlfile: str) -> dict:
    with open(yamlfile, "r") as f:
        return yaml.safe_load(f)


def save_yaml(data, out_yamlfile: str, mode: str = "w"):
    if mode == "w" and os.path.exists(out_yamlfile):
        logger.warning(f"Warning: {out_yamlfile} already exists. Overwriting.")
    with open(out_yamlfile, mode) as f:
        yaml.dump(data, f)


def create_basedir(out: str) -> None:
    """Create the parent directory of a path if necessary."""
    os.makedirs(os.path.dirname(out), exist_ok=True)


def warn_file_exists(out: str) -> None:
    """Notify user if an output file or directory already exists."""
    if os.path.exists(out):
        logger.warning(f"Warning: {out} already exists. Overwriting.")


def run_command(cmd: str) -> tuple[str, str]:
    try:
        cmd_out = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, check=True
        )
    except subprocess.CalledProcessError as e:
        raise ValueError(f"Command {cmd} failed:\n{e.stderr}")

    return cmd_out.stdout, cmd_out.stderr


def R_from_eman(a: np.ndarray, b: np.ndarray, y: np.ndarray) -> np.ndarray:
    a *= np.pi / 180.0
    b *= np.pi / 180.0
    y *= np.pi / 180.0
    ca, sa = np.cos(a), np.sin(a)
    cb, sb = np.cos(b), np.sin(b)
    cy, sy = np.cos(y), np.sin(y)
    Ra = np.array([[ca, -sa, 0], [sa, ca, 0], [0, 0, 1]])
    Rb = np.array([[1, 0, 0], [0, cb, -sb], [0, sb, cb]])
    Ry = np.array(([cy, -sy, 0], [sy, cy, 0], [0, 0, 1]))
    R = np.dot(np.dot(Ry, Rb), Ra)
    # handling EMAN convention mismatch for where the origin of an image is (bottom right vs top right)
    R[0, 1] *= -1
    R[1, 0] *= -1
    R[1, 2] *= -1
    R[2, 1] *= -1
    return R


def R_from_relion(euler: np.ndarray) -> np.ndarray:
    """Produce a rotation matrix from Euler angles returned by RELION."""

    a = euler[:, 0] * np.pi / 180.0
    b = euler[:, 1] * np.pi / 180.0
    y = euler[:, 2] * np.pi / 180.0
    nsamp = euler.shape[0]
    ca, sa = np.cos(a), np.sin(a)
    cb, sb = np.cos(b), np.sin(b)
    cy, sy = np.cos(y), np.sin(y)

    r_amat = np.array(
        [
            [ca, -sa, np.repeat(0, nsamp)],
            [sa, ca, np.repeat(0, nsamp)],
            [np.repeat(0, nsamp), np.repeat(0, nsamp), np.repeat(1, nsamp)],
        ],
    )
    r_bmat = np.array(
        [
            [cb, np.repeat(0, nsamp), -sb],
            [np.repeat(0, nsamp), np.repeat(1, nsamp), np.repeat(0, nsamp)],
            [sb, np.repeat(0, nsamp), cb],
        ]
    )
    r_ymat = np.array(
        [
            [cy, -sy, np.repeat(0, nsamp)],
            [sy, cy, np.repeat(0, nsamp)],
            [np.repeat(0, nsamp), np.repeat(0, nsamp), np.repeat(1, nsamp)],
        ]
    )

    rmat = np.matmul(np.matmul(r_ymat.T, r_bmat.T), r_amat.T)
    rmat[:, 0, 2] *= -1
    rmat[:, 2, 0] *= -1
    rmat[:, 1, 2] *= -1
    rmat[:, 2, 1] *= -1

    return rmat


def R_from_relion_scipy(euler_: np.ndarray, degrees: bool = True) -> np.ndarray:
    """Nx3 array of RELION euler angles to rotation matrix"""
    from scipy.spatial.transform import Rotation as RR

    euler = euler_.copy()
    if euler.shape == (3,):
        euler = euler.reshape(1, 3)
    euler[:, 0] += 90
    euler[:, 2] -= 90
    f = np.ones((3, 3))
    f[0, 1] = -1
    f[1, 0] = -1
    f[1, 2] = -1
    f[2, 1] = -1
    rot = RR.from_euler("zxz", euler, degrees=degrees).as_matrix() * f
    return rot


def R_to_relion_scipy(rot: np.ndarray, degrees: bool = True) -> np.ndarray:
    """Nx3x3 rotation matrices to RELION euler angles"""
    from scipy.spatial.transform import Rotation as RR

    if rot.shape == (3, 3):
        rot = rot.reshape(1, 3, 3)
    assert len(rot.shape) == 3, "Input must have dim Nx3x3"
    f = np.ones((3, 3))
    f[0, 1] = -1
    f[1, 0] = -1
    f[1, 2] = -1
    f[2, 1] = -1
    euler = RR.from_matrix(rot * f).as_euler("zxz", degrees=True)
    euler[:, 0] -= 90
    euler[:, 2] += 90
    euler += 180
    euler %= 360
    euler -= 180
    if not degrees:
        euler *= np.pi / 180
    return euler


def xrot(tilt_deg):
    """Return rotation matrix associated with rotation over the x-axis"""
    theta = tilt_deg * np.pi / 180
    tilt = np.array(
        [
            [1.0, 0.0, 0.0],
            [0, np.cos(theta), -np.sin(theta)],
            [0, np.sin(theta), np.cos(theta)],
        ]
    )
    return tilt


@memoized
def _zero_sphere_helper(D: int) -> Tuple[np.ndarray, np.ndarray]:
    xx = np.linspace(-1, 1, D, endpoint=True if D % 2 == 1 else False)
    z, y, x = np.meshgrid(xx, xx, xx, indexing="ij")
    coords = np.stack((x, y, z), -1)
    r = np.sum(coords**2, axis=-1) ** 0.5
    retval = np.where(r > 1)
    return retval


def zero_sphere(vol: np.ndarray) -> np.ndarray:
    """Zero values of @vol outside the sphere"""
    assert len(set(vol.shape)) == 1, "volume must be a cube"
    D = vol.shape[0]
    tmp = _zero_sphere_helper(D)
    logger.debug("Zeroing {} pixels".format(len(tmp[0])))
    vol[tmp] = 0
    return vol


def assert_pkl_close(pkl_a: str, pkl_b: str, atol: float = 1e-4) -> None:
    a = pickle.load(open(pkl_a, "rb"))
    b = pickle.load(open(pkl_b, "rb"))
    if isinstance(a, tuple):
        for _a, _b in zip(a, b):
            assert np.linalg.norm(_a - _b) < atol
    else:
        assert np.linalg.norm(a - b) < atol
