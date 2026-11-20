"""Compact "columns" wire format for large Plotly scatter figures.

The scatter endpoints normally return a whole Plotly figure as JSON: per-point
coordinates, one hex colour string per point, and a ``customdata`` row per point
that includes a preformatted hover string. For ~1e5 particles that is 5-8 MB per
colour or axis change, most of it text that repeats a handful of values.

With ``wire="columns"`` the same figure is sent as:

* a small JSON header holding the figure (layout, trace attributes) in which every
  per-point array is replaced by ``{"$col": i}``, plus a description of each column;
* the columns themselves as raw little-endian buffers after the header.

Columns are stored once (the colour values used for ``marker.color`` and for the
hover text share one buffer), integer and categorical columns use the smallest
integer type, repeated strings are dictionary-encoded, and the continuous hover
text is rebuilt in the browser from its values. ``static/js/cryo_figure_columns.js``
turns the response back into the figure object the page code already expects, so
only the fetch call changes on each page.

Columns that do not depend on the colour covariate (coordinates, particle index,
table row) are tagged ``base``. A response carries a ``base_key``; when the client
sends it back as ``have`` and it still matches, those columns are left out and the
browser reuses the ones it holds. A colour change then transfers only the colour.

Layout of a response body::

    b"CDFC1\\n"  | uint32 header_len | header JSON (utf-8) | pad to 8 | buffers...

Buffer offsets in the header are relative to the start of the buffer section and
are multiples of 8, so the browser can view them in place as typed arrays.

:func:`decode_figure_columns` is the Python reference decoder (used by tests).
"""

from __future__ import annotations

import hashlib
import json
import struct
from typing import Any, Sequence

import numpy as np
import pandas as pd

MAGIC = b"CDFC1\n"
MIMETYPE = "application/vnd.cryodrgn.figure-columns"

_DICT_MAX_CATEGORIES = 65535
_NUMERIC_DICT_MAX_CATEGORIES = 256


class G5Text:
    """Per-point hover line ``prefix + '%.5g' % value`` (``missing`` when not finite).

    Sent as its values plus the two strings; the browser formats the text.
    """

    __slots__ = ("prefix", "missing", "values")

    def __init__(self, prefix: str, missing: str, values: np.ndarray) -> None:
        self.prefix = str(prefix)
        self.missing = str(missing)
        self.values = np.asarray(values, dtype=np.float64)

    def to_list(self) -> list[str]:
        """The strings this column stands for (same rule as the browser)."""
        finite = np.isfinite(self.values)
        formatted = np.char.mod("%.5g", self.values)
        return np.where(
            finite, np.char.add(self.prefix, formatted), self.missing
        ).tolist()


class NullableFloats:
    """Float column whose non-finite entries are ``null`` in the figure."""

    __slots__ = ("values",)

    def __init__(self, values: np.ndarray) -> None:
        self.values = np.asarray(values, dtype=np.float64)


class Rows:
    """``customdata``: one row per point, assembled from named member columns.

    ``base`` lists the member names that do not depend on the colour covariate.
    """

    __slots__ = ("columns", "names", "base")

    def __init__(
        self,
        columns: Sequence[Any],
        names: Sequence[str],
        base: Sequence[str] = (),
    ) -> None:
        if len(columns) != len(names):
            raise ValueError("Rows needs one name per column.")
        self.columns = list(columns)
        self.names = [str(n) for n in names]
        self.base = frozenset(base)


def base_key(*parts: Any) -> str:
    """Short digest identifying the colour-independent columns of a figure."""
    h = hashlib.sha1(repr(parts).encode("utf-8")).hexdigest()
    return h[:20]


# ---------------------------------------------------------------------------
# Encoder
# ---------------------------------------------------------------------------


class _Writer:
    def __init__(self, have_base: bool) -> None:
        self.have_base = have_base
        self.cols: list[dict[str, Any]] = []
        self.buffers: list[bytes] = []
        self.offset = 0
        self._dedupe: dict[tuple, int] = {}
        self._buffer_offsets: dict[tuple, int] = {}

    def _add_buffer(self, arr: np.ndarray) -> dict[str, Any]:
        arr = np.ascontiguousarray(arr)
        if arr.dtype.byteorder == ">":
            arr = arr.astype(arr.dtype.newbyteorder("<"))
        raw = arr.tobytes()
        # Identical buffers are stored once: e.g. the category codes behind the
        # per-point hex colour, filter key and hover text of a discrete covariate.
        bkey = (arr.dtype.str, hashlib.sha1(raw).digest())
        off = self._buffer_offsets.get(bkey)
        if off is None:
            pad = (-len(raw)) % 8
            off = self.offset
            self.buffers.append(raw + b"\0" * pad)
            self.offset += len(raw) + pad
            self._buffer_offsets[bkey] = off
        return {
            "dtype": arr.dtype.kind + str(arr.dtype.itemsize),
            "offset": off,
            "n": int(arr.shape[0]),
        }

    def _push(self, spec: dict[str, Any]) -> int:
        self.cols.append(spec)
        return len(self.cols) - 1

    @staticmethod
    def _int_dtype(lo: int, hi: int) -> np.dtype | None:
        for dt in (np.int8, np.int16, np.int32):
            info = np.iinfo(dt)
            if info.min <= lo and hi <= info.max:
                return np.dtype(dt)
        return None

    def _numeric(
        self,
        arr: np.ndarray,
        *,
        coords: bool,
        try_dict: bool,
    ) -> dict[str, Any]:
        n = int(arr.shape[0])
        if arr.dtype.kind == "b":
            arr = arr.astype(np.int8)
        if arr.dtype.kind in "iu":
            if n and int(arr[0]) == 0 and int(arr[-1]) == n - 1:
                if np.array_equal(arr, np.arange(n, dtype=arr.dtype)):
                    return {"kind": "arange", "n": n}
            lo = int(arr.min()) if n else 0
            hi = int(arr.max()) if n else 0
            dt = self._int_dtype(lo, hi)
            if dt is None:
                arr = arr.astype(np.float64)
            else:
                arr = arr.astype(dt, copy=False)
        elif arr.dtype.kind == "f":
            if try_dict and n > 64:
                probe = arr[:: max(1, n // 1024)]
                if np.unique(probe).size <= _NUMERIC_DICT_MAX_CATEGORIES:
                    uniq, inv = np.unique(arr, return_inverse=True)
                    if (
                        uniq.size <= _NUMERIC_DICT_MAX_CATEGORIES
                        and not np.isnan(uniq).any()
                    ):
                        spec = self._add_buffer(inv.astype(np.uint8))
                        spec.update(kind="dict", table=[float(u) for u in uniq])
                        return spec
            if arr.dtype.itemsize > 4:
                a32 = arr.astype(np.float32)
                if coords or np.array_equal(
                    a32.astype(np.float64), arr, equal_nan=True
                ):
                    arr = a32
        else:
            raise TypeError(f"Unsupported numeric dtype {arr.dtype!r}.")
        spec = self._add_buffer(arr)
        spec["kind"] = "num"
        return spec

    def _objects(self, values: Any) -> dict[str, Any]:
        ser = pd.Series(values, dtype=object)
        n = int(len(ser))
        codes, uniques = pd.factorize(ser, use_na_sentinel=False)
        if len(uniques) == 1:
            v = uniques[0]
            return {"kind": "const", "n": n, "value": None if pd.isna(v) else v}
        if len(uniques) > _DICT_MAX_CATEGORIES:
            return {"kind": "json", "values": [None if pd.isna(v) else v for v in ser]}
        table = [
            None if (not isinstance(v, str) and pd.isna(v)) else v for v in uniques
        ]
        table = [t.item() if isinstance(t, np.generic) else t for t in table]
        dt = np.uint8 if len(uniques) <= 256 else np.uint16
        spec = self._add_buffer(np.asarray(codes).astype(dt))
        spec.update(kind="dict", table=table)
        return spec

    def column(
        self,
        value: Any,
        *,
        name: str | None = None,
        base: bool = False,
        coords: bool = False,
    ) -> int:
        """Register one per-point array; returns its column index."""
        if base and self.have_base and name:
            return self._push({"kind": "cached", "name": name})

        if isinstance(value, G5Text):
            src = self.column(value.values)
            return self._push(
                {
                    "kind": "g5",
                    "src": src,
                    "prefix": value.prefix,
                    "missing": value.missing,
                }
            )

        if isinstance(value, NullableFloats):
            value = value.values
        if isinstance(value, (pd.Series, pd.Index)):
            value = value.to_numpy()
        if not isinstance(value, np.ndarray):
            value = np.asarray(value, dtype=object)

        if value.dtype.kind in "biuf":
            key = (
                "num",
                value.dtype.str,
                int(value.shape[0]),
                hashlib.sha1(np.ascontiguousarray(value).tobytes()).digest(),
                bool(coords),
            )
            hit = self._dedupe.get(key)
            if hit is not None and not base:
                return hit
            spec = self._numeric(
                value,
                coords=coords,
                try_dict=not coords and not base,
            )
        else:
            key = None
            spec = self._objects(value)
        if base and name:
            spec["name"] = name
            spec["base"] = True
        idx = self._push(spec)
        if key is not None and not base:
            self._dedupe[key] = idx
        return idx


def _is_point_array(value: Any, n_points: int) -> bool:
    if isinstance(value, (np.ndarray, pd.Series)):
        return value.ndim == 1 and len(value) == n_points
    if isinstance(value, list):
        return len(value) == n_points and n_points > 16
    return False


def encode_figure_columns(
    layout: dict[str, Any],
    trace: dict[str, Any],
    *,
    n_points: int,
    key: str,
    have: str | None = None,
    coord_keys: Sequence[str] = ("x", "y", "z"),
) -> bytes:
    """Serialize one-trace figure ``{"data": [trace], "layout": layout}``.

    ``trace`` is a plain dict of Plotly trace attributes. Values that are 1-D
    arrays of length ``n_points`` (directly or one level down, e.g. in ``marker``)
    become columns, as does a :class:`Rows` ``customdata``. ``coord_keys`` name the
    coordinate arrays: sent as float32 and tagged ``base``.
    """
    from plotly.utils import PlotlyJSONEncoder

    w = _Writer(have_base=bool(have) and have == key)

    def conv(obj: Any, name: str | None, depth: int) -> Any:
        if isinstance(obj, Rows):
            members = [
                w.column(c, name=f"cd.{nm}", base=nm in obj.base)
                for c, nm in zip(obj.columns, obj.names)
            ]
            nulls = [isinstance(c, NullableFloats) for c in obj.columns]
            return {
                "$col": w._push({"kind": "rows", "cols": members, "nan_null": nulls})
            }
        if _is_point_array(obj, n_points):
            is_coord = depth == 0 and name in coord_keys
            return {"$col": w.column(obj, name=name, base=is_coord, coords=is_coord)}
        if isinstance(obj, dict) and depth < 3:
            return {k: conv(v, k, depth + 1) for k, v in obj.items()}
        if isinstance(obj, (pd.Series, pd.Index)):
            return obj.to_numpy().tolist()
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.generic):
            return obj.item()
        return obj

    fig = {"data": [conv(trace, None, -1)], "layout": layout}
    header = {
        "v": 1,
        "n": int(n_points),
        "base_key": key,
        "base_omitted": bool(w.have_base),
        "fig": fig,
        "cols": w.cols,
    }
    hbytes = json.dumps(header, cls=PlotlyJSONEncoder, separators=(",", ":")).encode(
        "utf-8"
    )
    prefix = MAGIC + struct.pack("<I", len(hbytes)) + hbytes
    pad = (-len(prefix)) % 8
    return b"".join([prefix, b"\0" * pad, *w.buffers])


# ---------------------------------------------------------------------------
# Reference decoder (tests, debugging)
# ---------------------------------------------------------------------------


def _g5(value: float) -> str:
    return "%.5g" % value


def decode_figure_columns(
    body: bytes, cached: dict[str, Any] | None = None
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Inverse of :func:`encode_figure_columns`.

    Returns ``(figure, info)``; ``info`` has ``base_key`` and ``base`` (name →
    decoded base column, to pass back as ``cached`` for a later response that
    omits them). Numeric columns come back as numpy arrays; row and string
    columns as lists, mirroring what the browser builds.
    """
    if body[: len(MAGIC)] != MAGIC:
        raise ValueError("Not a figure-columns response.")
    (hlen,) = struct.unpack_from("<I", body, len(MAGIC))
    h0 = len(MAGIC) + 4
    header = json.loads(body[h0 : h0 + hlen].decode("utf-8"))
    data0 = h0 + hlen
    data0 += (-data0) % 8
    specs = header["cols"]
    decoded: dict[int, Any] = {}
    base_out: dict[str, Any] = {}

    def buf(spec: dict[str, Any]) -> np.ndarray:
        dt = np.dtype("<" + spec["dtype"])
        return np.frombuffer(
            body, dtype=dt, count=spec["n"], offset=data0 + spec["offset"]
        )

    def nulled(values: Any) -> list[Any]:
        return [
            None if (isinstance(v, float) and not np.isfinite(v)) else v
            for v in (values.tolist() if isinstance(values, np.ndarray) else values)
        ]

    def col(i: int) -> Any:
        if i in decoded:
            return decoded[i]
        s = specs[i]
        kind = s["kind"]
        if kind == "cached":
            if not cached or s["name"] not in cached:
                raise KeyError(f"Base column {s['name']!r} is not cached.")
            out = cached[s["name"]]
            base_out[s["name"]] = out
        elif kind == "arange":
            out = np.arange(s["n"], dtype=np.int32)
        elif kind == "num":
            out = buf(s)
        elif kind == "dict":
            table = s["table"]
            out = [table[c] for c in buf(s)]
        elif kind == "const":
            out = [s["value"]] * s["n"]
        elif kind == "json":
            out = list(s["values"])
        elif kind == "g5":
            out = [
                s["missing"] if v is None else s["prefix"] + _g5(v)
                for v in nulled(col(s["src"]))
            ]
        elif kind == "rows":
            members = []
            for j, null in zip(s["cols"], s["nan_null"]):
                m = col(j)
                if null:
                    m = nulled(m)
                members.append(m.tolist() if isinstance(m, np.ndarray) else m)
            out = [list(r) for r in zip(*members)]
        else:
            raise ValueError(f"Unknown column kind {kind!r}.")
        if s.get("base") and s.get("name"):
            base_out[s["name"]] = out
        decoded[i] = out
        return out

    def hydrate(obj: Any) -> Any:
        if isinstance(obj, dict):
            if len(obj) == 1 and "$col" in obj:
                return col(int(obj["$col"]))
            return {k: hydrate(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [hydrate(v) for v in obj]
        return obj

    fig = hydrate(header["fig"])
    return fig, {
        "base_key": header["base_key"],
        "base_omitted": header["base_omitted"],
        "base": base_out,
        "n": header["n"],
    }
