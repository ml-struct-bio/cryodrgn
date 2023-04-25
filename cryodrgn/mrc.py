import sys
import os
import struct
from collections import OrderedDict
from typing import Any, Optional, Tuple
import numpy as np
import cryodrgn.types as types
from cryodrgn import __version__

# See ref:
# MRC2014: Extensions to the MRC format header for electron cryo-microscopy and tomography
# And:
# https://www.ccpem.ac.uk/mrc_format/mrc2014.php

DTYPE_FOR_MODE = {
    0: np.int8,
    1: np.int16,
    2: np.float32,
    3: "2h",  # complex number from 2 shorts
    4: np.complex64,
    6: np.uint16,
    12: np.float16,
    16: "3B",
}  # RBG values
MODE_FOR_DTYPE = {vv: kk for kk, vv in DTYPE_FOR_MODE.items()}

MACHST_OFFSET = 212
MACHST_FOR_ENDIANNESS = {"<": b"\x44\x44\x00\x00", ">": b"\x11\x11\x00\x00"}
ENDIANNESS_FOR_MACHST = {v: k for k, v in MACHST_FOR_ENDIANNESS.items()}


class MRCHeader:
    """MRC header class"""

    # Class attribute denoting endianness (machst field in MRC file format spec)
    #   '<' = little-endian; '>' = big-endian; '=' = machine-native
    # Note that instances of this class may override it on '.parse()' classmethod invocation,
    # depending on what they actually find in the MRC file header at MACHST_OFFSET.
    endianness = "="

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

    def __init__(self, header_values, extended_header=b"", endianness="="):
        self.fields = OrderedDict(zip(self.FIELDS, header_values))
        self.extended_header = extended_header
        self.D = self.fields["nx"]
        assert endianness in ("=", "<", ">")
        self.endianness = endianness
        self.dtype = np.dtype(DTYPE_FOR_MODE[self.fields["mode"]]).newbyteorder(
            endianness
        )

    def __str__(self):
        return f"Header: {self.fields}\nExtended header: {self.extended_header}"

    @classmethod
    def parse(cls, fname):
        with open(fname, "rb") as f:
            f.seek(MACHST_OFFSET)
            # Look for valid machst field; assume '=' if invalid
            machst = f.read(4)
            endianness = ENDIANNESS_FOR_MACHST.get(machst, "=")

            f.seek(0)
            STRUCT = struct.Struct(
                endianness + cls.FSTR
            )  # prepend endianness specifier to python struct specification
            header = cls(STRUCT.unpack(f.read(1024)), endianness=endianness)

            # Older versions of MRCHeader in cryoDRGN had incorrect cmap and stamp fields.
            # Fix these before proceeding.
            header.fields["cmap"] = b"MAP "
            if header.endianness == "=":
                endianness = {"little": "<", "big": ">"}[sys.byteorder]
            else:
                endianness = header.endianness
            header.fields["stamp"] = MACHST_FOR_ENDIANNESS[endianness]

            extbytes = header.fields["next"]
            extended_header = f.read(extbytes)
            header.extended_header = extended_header
        return header

    @classmethod
    def make_default_header(
        cls, data, is_vol=True, Apix=1.0, xorg=0.0, yorg=0.0, zorg=0.0
    ):
        nz, ny, nx = data.shape
        ispg = 1 if is_vol else 0
        if is_vol:
            dmin, dmax, dmean, rms = data.min(), data.max(), data.mean(), data.std()
        else:  # use undefined values for image stacks
            dmin, dmax, dmean, rms = -1, -2, -3, -1
        vals = [
            nx,
            ny,
            nz,
            2,  # mode = 2 for 32-bit float
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
            MACHST_FOR_ENDIANNESS["<" if sys.byteorder == "little" else ">"],
            rms,  # rms
            1,  # nlabl
            ("cryoDRGN " + __version__[:80]).ljust(80, " ").encode("ascii")
            + b"\x00" * 720,  # Use one 80-char label, leave remaining 9 blank
        ]
        return cls(vals)

    def write(self, fh):
        STRUCT = struct.Struct(self.endianness + self.FSTR)
        buf = STRUCT.pack(*list(self.fields.values()))
        fh.write(buf)
        fh.write(self.extended_header)

    def get_apix(self):
        return self.fields["xlen"] / self.fields["nx"]

    def update_apix(self, Apix):
        self.fields["xlen"] = self.fields["nx"] * Apix
        self.fields["ylen"] = self.fields["ny"] * Apix
        self.fields["zlen"] = self.fields["nz"] * Apix

    def get_origin(self):
        return self.fields["xorg"], self.fields["yorg"], self.fields["zorg"]

    def update_origin(self, xorg, yorg, zorg):
        self.fields["xorg"] = xorg
        self.fields["yorg"] = yorg
        self.fields["zorg"] = zorg


def parse_header(fname: str) -> MRCHeader:
    return MRCHeader.parse(fname)


class LazyImage:
    """On-the-fly image loading"""

    def __init__(self, fname: str, shape: Tuple[int, int], dtype: Any, offset: int):
        self.fname = fname
        self.shape = shape
        self.dtype = dtype
        self.offset = offset

    def get(self) -> np.ndarray:
        with open(self.fname) as f:
            f.seek(self.offset)
            image = np.fromfile(
                f, dtype=self.dtype, count=np.product(self.shape)
            ).reshape(self.shape)
        return image


def parse_mrc_list(txtfile: str, lazy: bool = False) -> types.ImageArray:
    lines = open(txtfile, "r").readlines()

    def abspath(f):
        if os.path.isabs(f):
            return f
        base = os.path.dirname(os.path.abspath(txtfile))
        return os.path.join(base, f)

    lines = [abspath(x) for x in lines]
    if not lazy:
        arrays = []
        for line in lines:
            array = parse_mrc(line.strip(), lazy=False)[0]
            arrays.append(array)
        particles = np.vstack(arrays)
    else:
        particles = [img for x in lines for img in parse_mrc(x.strip(), lazy=True)[0]]
    return particles


def parse_mrc(fname: str, lazy: bool = False) -> Tuple[types.ImageArray, MRCHeader]:  # type: ignore
    # parse the header
    header = MRCHeader.parse(fname)

    # get the number of bytes in extended header
    extbytes = header.fields["next"]
    start = 1024 + extbytes  # start of image data

    dtype = header.dtype
    nz, ny, nx = header.fields["nz"], header.fields["ny"], header.fields["nx"]

    # load all in one block
    if not lazy:
        with open(fname, "rb") as fh:
            fh.read(start)  # skip the header + extended header
            array = np.fromfile(fh, dtype=dtype).reshape((nz, ny, nx))

    # or list of LazyImages
    else:
        stride = dtype.itemsize * ny * nx
        array = [
            LazyImage(fname, (ny, nx), dtype, start + i * stride) for i in range(nz)
        ]
    return array, header  # type: ignore


def write(
    fname,
    array: np.ndarray,
    header: Optional[MRCHeader] = None,
    Apix: float = 1.0,
    xorg: float = 0.0,
    yorg: float = 0.0,
    zorg: float = 0.0,
    is_vol: Optional[bool] = None,
):
    # get a default header
    if header is None:
        if is_vol is None:
            is_vol = (
                True if len(set(array.shape)) == 1 else False
            )  # Guess whether data is vol or image stack
        header = MRCHeader.make_default_header(array, is_vol, Apix, xorg, yorg, zorg)

    # write the header
    f = open(fname, "wb")
    header.write(f)

    f.write(array.astype(header.dtype).tobytes())
