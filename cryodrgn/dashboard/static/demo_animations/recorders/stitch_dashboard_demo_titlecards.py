#!/usr/bin/env python3
"""
Stitch dashboard demo GIFs into a single MP4.

This is the canonical entry point living under
`cryodrgn/dashboard/static/demo_animations/recorders/` alongside the other demo
assets. The implementation is currently delegated to the upstream script in
the repo-root `scripts/` directory.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path


def _load_upstream():
    repo_root = Path(__file__).resolve().parents[5]
    upstream = repo_root / "scripts" / "stitch_dashboard_demo_titlecards.py"
    if not upstream.exists():
        raise FileNotFoundError(f"Upstream stitch script missing: {upstream}")

    spec = importlib.util.spec_from_file_location(
        "stitch_dashboard_demo_titlecards_upstream", str(upstream)
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load upstream stitch script: {upstream}")

    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    mod = _load_upstream()
    return int(mod.main())


if __name__ == "__main__":
    raise SystemExit(main())
