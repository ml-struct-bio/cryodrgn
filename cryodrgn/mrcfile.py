"""Utilities for reading and writing .mrc/.mrcs files.

Example usage
-------------
> from cryodrgn.mrcfile import parse_mrc, write_mrc
> img_array, header = parse_mrc("particles.mrcs")
> img_array = img_array[::5, :, :]  # take every fifth image
> write_mrc("new-particles.mrcs", img_array, header)

"""
import sys
import struct
from collections import OrderedDict
from typing import Tuple, Union, Optional, Callable
from typing_extensions import Self
import numpy as np
import torch

import logging

logger = logging.getLogger(__name__)


class MRCHeader:
    """A class for representing the headers of .mrc files which store metadata.

    See ref:
        MRC2014: Extensions to the MRC format header for electron cryo-microscopy and
                 tomography
    and:
        https://www.ccpem.ac.uk/mrc_format/mrc2014.php

    """

    ENDIANNESS = "="
    FIELDS = [
        "nx",
        "ny",
        "nz",  # int
        "mode",  # int
        "nxstart",
        "nystart",
        "nzstart",  # int
        "mx",
        "my",
        "mz",  # int
        "xlen",
        "ylen",
        "zlen",  # float
        "alpha",
        "beta",
        "gamma",  # float
        "mapc",
        "mapr",
        "maps",  # int
        "amin",
        "amax",
        "amean",  # float
        "ispg",
        "next",
        "creatid",  # int, int, short, [pad 10]
        "nversion",  # int, [pad 20]
        "nint",
        "nreal",  # short, [pad 20]
        "imodStamp",
        "imodFlags",  # int
        "idtype",
        "lens",
        "nd1",
        "nd2",
        "vd1",
        "vd2",  # short
        "tilt_ox",
        "tilt_oy",
        "tilt_oz",  # float
        "tilt_cx",
        "tilt_cy",
        "tilt_cz",  # float
        "xorg",
        "yorg",
        "zorg",  # float
        "cmap",
        "stamp",
        "rms",  # char[4], float
        "nlabl",
        "labels",
    ]  # int, char[10][80]
    FSTR = "3ii3i3i3f3f3i3f2ih10xi16x2h20x2i6h6f3f4s4sfi800s"

    # Mappings for number formats used by .mrc files to number formats used by numpy
    # Note that (u)int32 is treated equivalent to float32 here
    DTYPE_FOR_MODE = {
        0: np.uint8,
        1: np.int16,
        2: np.float32,
        3: "2h",  # complex number from 2 shorts
        4: np.complex64,
        6: np.uint16,
        12: np.float16,
        16: "3B",
        17: np.int8,
    }  # RBG values
    MODE_FOR_DTYPE = {vv: kk for kk, vv in DTYPE_FOR_MODE.items()}

    MACHST_OFFSET = 213
    MACHST_FOR_ENDIANNESS = {"<": b"\x44\x44\x00\x00", ">": b"\x11\x11\x00\x00"}
    ENDIANNESS_FOR_MACHST = {v: k for k, v in MACHST_FOR_ENDIANNESS.items()}

    def __init__(self, header_values, extended_header=b""):
        self.fields = OrderedDict(zip(self.FIELDS, header_values))
        self.extended_header = extended_header
        self.D = self.fields["nx"]
        self.N = self.fields["nz"]

        if self.fields["mode"] not in self.DTYPE_FOR_MODE:
            raise ValueError(
                f"This file contains a Data Type mode label `{self.fields['mode']}` "
                f"not found in the dictionary of recognized mode to dtype mappings:\n"
                f"{self.DTYPE_FOR_MODE}"
            )
        self.dtype = self.DTYPE_FOR_MODE[self.fields["mode"]]

    def __str__(self):
        return f"Header: {self.fields}\nExtended header: {self.extended_header}"

    @classmethod
    def parse(cls, fname: str) -> Self:
        """Create a `MRCHeader` object by reading in the header from a .mrc(s) file."""

        with open(fname, "rb") as f:
            f.seek(cls.MACHST_OFFSET)
            cls.ENDIANNESS = cls.ENDIANNESS_FOR_MACHST.get(f.read(2), "=")

            f.seek(0)
            # prepend endianness specifier to python struct specification
            STRUCT = struct.Struct(cls.ENDIANNESS + cls.FSTR)
            header = cls(STRUCT.unpack(f.read(1024)))

            extbytes = header.fields["next"]
            extended_header = f.read(extbytes)
            header.extended_header = extended_header

        return header

    @classmethod
    def make_default_header(
        cls,
        nz: Optional[int] = None,
        ny: Optional[int] = None,
        nx: Optional[int] = None,
        data: Optional[Union[np.ndarray, torch.Tensor]] = None,
        dtype: Optional[Union[str, np.dtype]] = None,
        is_vol: bool = True,
        Apix: float = 1.0,
        xorg: float = 0.0,
        yorg: float = 0.0,
        zorg: float = 0.0,
    ) -> Self:
        if dtype is not None:
            data_dtype = np.dtype(dtype)
        else:
            data_dtype = np.dtype("float32")  # default to np.float 32 mode

        if data is not None:
            nz, ny, nx = data.shape

            if dtype is None:
                if isinstance(data, torch.Tensor):
                    try:
                        data_dtype = np.dtype(str(data.dtype).split(".")[1])
                    except TypeError:
                        data_dtype = np.dtype("float32")
                else:
                    data_dtype = data.dtype

        if data_dtype in cls.MODE_FOR_DTYPE:
            use_mode = cls.MODE_FOR_DTYPE[data_dtype]
        elif data_dtype.type in cls.MODE_FOR_DTYPE:
            use_mode = cls.MODE_FOR_DTYPE[data_dtype.type]
        else:
            use_mode = 2

        assert nz is not None
        assert ny is not None
        assert nx is not None
        ispg = 1 if is_vol else 0

        if is_vol:
            if data is None:
                raise ValueError("If is_vol=True, data array must be specified")

            dmin, dmax, dmean, rms = data.min(), data.max(), data.mean(), data.std()
        else:  # use undefined values for image stacks
            dmin, dmax, dmean, rms = -1, -2, -3, -1

        vals = [
            nx,
            ny,
            nz,
            use_mode,  # mode = 2 for 32-bit float
            0,
            0,
            0,  # nxstart, nystart, nzstart
            nx,
            ny,
            nz,  # mx, my, mz
            Apix * nx,
            Apix * ny,
            Apix * nz,  # cella
            90.0,
            90.0,
            90.0,  # cellb
            1,
            2,
            3,  # mapc, mapr, maps
            dmin,
            dmax,
            dmean,
            ispg,
            0,  # exthd_size
            0,  # creatid
            20140,  # nversion
            0,
            0,  # nint, nreal
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            xorg,
            yorg,
            zorg,
            b"MAP ",
            cls.MACHST_FOR_ENDIANNESS["<" if sys.byteorder == "little" else ">"],
            rms,  # rms
            0,  # nlabl
            b"\x00" * 800,  # labels
        ]

        return cls(vals)

    def write(self, fh):
        STRUCT = struct.Struct(self.FSTR)
        buf = STRUCT.pack(*list(self.fields.values()))
        fh.write(buf)
        fh.write(self.extended_header)

    @property
    def apix(self) -> float:
        return round(self.fields["xlen"] / self.fields["nx"], 6)

    @apix.setter
    def apix(self, value: float) -> None:
        self.fields["xlen"] = self.fields["nx"] * value
        self.fields["ylen"] = self.fields["ny"] * value
        self.fields["zlen"] = self.fields["nz"] * value

    @property
    def origin(self) -> tuple[float, float, float]:
        return self.fields["xorg"], self.fields["yorg"], self.fields["zorg"]

    @origin.setter
    def origin(self, value: tuple[float, float, float]) -> None:
        self.fields["xorg"] = value[0]
        self.fields["yorg"] = value[1]
        self.fields["zorg"] = value[2]


def parse_mrc(fname: str) -> Tuple[np.ndarray, MRCHeader]:
    """Read in the array of data values and the header data stored in a .mrc(s) file."""
    header = MRCHeader.parse(fname)

    # get the number of bytes in extended header
    extbytes = header.fields["next"]
    start = 1024 + extbytes  # start of image data

    dtype = header.dtype
    nz, ny, nx = header.fields["nz"], header.fields["ny"], header.fields["nx"]

    with open(fname, "rb") as fh:
        fh.read(start)  # skip the header + extended header
        array = np.fromfile(fh, dtype=dtype).reshape((nz, ny, nx))

    return array, header


def get_mrc_header(
    array: Union[np.ndarray, torch.Tensor], is_vol: Optional[bool] = None, **header_args
) -> MRCHeader:
    """Create the default header corresponding to this image data array."""

    if is_vol is None:
        # If necessary, guess whether data is vol or image stack
        is_vol = len(set(array.shape)) == 1

    header = MRCHeader.make_default_header(
        nz=None,
        ny=None,
        nx=None,
        data=array,
        is_vol=is_vol,
        **header_args,
    )
    return header


def fix_mrc_header(header: MRCHeader) -> MRCHeader:
    """Fix older versions of MRCHeader with incorrect `cmap` and `stamp` fields."""
    header.fields["cmap"] = b"MAP "

    if header.ENDIANNESS == "=":
        endianness = {"little": "<", "big": ">"}[sys.byteorder]
    else:
        endianness = header.ENDIANNESS

    header.fields["stamp"] = header.MACHST_FOR_ENDIANNESS[endianness]

    return header


def write_mrc(
    filename: str,
    array: Union[np.ndarray, torch.Tensor],
    header: Optional[MRCHeader] = None,
    is_vol: Optional[bool] = None,
    transform_fn: Optional[Callable] = None,
    **header_args,
) -> None:
    """Save an image stack or volume to disk as an .mrc(s) file.

    Arguments
    ---------
    filename        Where the .mrc(s) will be saved.
    array           The image stack or volume to save to file.
    header          Optionally supply an MRCHeader instead of using the default one.
    is_vol          Don't infer whether this is a volume from the array itself.
    transform_fn    Apply this function to the array values before saving.
    header_args     Additional keyword arguments passed to `MRCHeader` if not using
                    your own header.

    """
    if header is None:
        header = get_mrc_header(array, is_vol, **header_args)
    else:
        if header_args:
            logger.warning(
                f"Passed header arguments {header_args} to `write_mrc` but these will "
                "not be used as header was also given!"
            )
        header = fix_mrc_header(header=header)

    if transform_fn is None:
        transform_fn = lambda chunk, indices: chunk  # noqa: E731

    new_dtype = np.dtype(header.dtype).newbyteorder(header.ENDIANNESS)  # type: ignore
    with open(filename, "wb") as f:
        header.write(f)
        indices = np.arange(array.shape[0])
        array = transform_fn(array, indices)

        if isinstance(array, torch.Tensor):
            array = np.array(array.cpu()).astype(new_dtype)

        assert isinstance(array, np.ndarray)
        f.write(array.tobytes())
