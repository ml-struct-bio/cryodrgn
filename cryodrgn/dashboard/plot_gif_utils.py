"""Assemble PNG frame sequences into animated GIFs for the dashboard."""

from __future__ import annotations

import base64
from io import BytesIO
from typing import Sequence

# Cream-ish paper background (matches Plotly layout paper_bgcolor in scatter3d).
_GIF_FLATTEN_BG = (250, 248, 244)


def png_base64_frames_to_gif_bytes(
    frames_b64: Sequence[str],
    *,
    durations_ms: Sequence[int] | int = 100,
    loop: int = 0,
) -> bytes:
    """Decode PNG base64 strings and build an animated GIF.

    ``frames_b64`` entries may be raw base64 or ``data:image/png;base64,...`` URLs.
    """
    try:
        from PIL import Image
    except ImportError as exc:  # pragma: no cover — matplotlib normally pulls in Pillow
        raise RuntimeError("Pillow (PIL) is required for GIF export.") from exc

    if len(frames_b64) < 2:
        raise ValueError("At least two PNG frames are required.")
    if len(frames_b64) > 120:
        raise ValueError("Too many frames (max 120).")

    if isinstance(durations_ms, int):
        durs: list[int] = [max(1, int(durations_ms))] * len(frames_b64)
    else:
        durs = [max(1, int(x)) for x in durations_ms]
        if len(durs) != len(frames_b64):
            raise ValueError("durations_ms must match the number of frames.")

    images_rgb: list[Image.Image] = []
    for raw in frames_b64:
        s = raw.strip()
        if "," in s and s.lower().startswith("data:"):
            s = s.split(",", 1)[1]
        png_bytes = base64.b64decode(s, validate=True)
        im = Image.open(BytesIO(png_bytes))
        if im.mode == "RGBA":
            bg = Image.new("RGB", im.size, _GIF_FLATTEN_BG)
            bg.paste(im, mask=im.split()[3])
            im_rgb = bg
        else:
            im_rgb = im.convert("RGB")
        images_rgb.append(im_rgb)

    out = BytesIO()
    first, *rest = images_rgb
    palettes: list[Image.Image] = [
        first.convert("P", palette=Image.ADAPTIVE, colors=256)
    ]
    for im in rest:
        palettes.append(im.convert("P", palette=Image.ADAPTIVE, colors=256))

    palettes[0].save(
        out,
        format="GIF",
        save_all=True,
        append_images=palettes[1:],
        duration=durs,
        loop=int(loop),
        optimize=False,
        disposal=2,
    )
    return out.getvalue()
