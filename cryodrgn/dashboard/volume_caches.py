"""In-memory caches for decoded volumes and long-running volume job state."""

from __future__ import annotations

import base64
import secrets
import shutil
import threading
import time


class VolumeMrcCache:
    """TTL-backed registry of decoded ``.mrc`` directories for montage / trajectory jobs."""

    def __init__(
        self,
        *,
        max_entries: int = 32,
        ttl_s: float = 7200.0,
    ) -> None:
        self._max_entries = max_entries
        self._ttl_s = ttl_s
        self._lock = threading.Lock()
        self._entries: dict[str, dict[str, object]] = {}

    @property
    def lock(self) -> threading.Lock:
        return self._lock

    @property
    def entries(self) -> dict[str, dict[str, object]]:
        """Raw token → metadata map (for tests and diagnostics)."""
        return self._entries

    def evict_unlocked(self, token: str) -> None:
        meta = self._entries.pop(token, None)
        if meta and meta.get("mrc_dir"):
            shutil.rmtree(meta["mrc_dir"], ignore_errors=True)

    def prune_unlocked(self) -> None:
        now = time.monotonic()
        dead = [
            tok
            for tok, meta in self._entries.items()
            if now - float(meta["t0"]) > self._ttl_s
        ]
        for tok in dead:
            self.evict_unlocked(tok)
        while len(self._entries) >= self._max_entries:
            oldest = min(self._entries.items(), key=lambda kv: float(kv[1]["t0"]))[0]
            self.evict_unlocked(oldest)

    def register(
        self,
        mrc_dir: str,
        vol_files: list[str],
        rows: tuple[int, ...],
    ) -> str:
        with self._lock:
            self.prune_unlocked()
            token = secrets.token_urlsafe(24)
            self._entries[token] = {
                "mrc_dir": mrc_dir,
                "vol_files": list(vol_files),
                "rows": rows,
                "t0": time.monotonic(),
            }
            return token

    def get_meta(self, token: str) -> dict[str, object] | None:
        with self._lock:
            meta = self._entries.get(token)
            if not meta:
                return None
            if time.monotonic() - float(meta["t0"]) > self._ttl_s:
                self.evict_unlocked(token)
                return None
            return meta

    def require_meta(
        self,
        token: str,
        *,
        rows_expected: tuple[int, ...] | None = None,
    ) -> dict[str, object]:
        meta = self.get_meta(token)
        if not meta:
            raise ValueError("Unknown or expired volume cache id.")
        if rows_expected is not None and meta["rows"] != rows_expected:
            raise ValueError("Montage rows do not match cached volumes.")
        return meta


class VolumeJobStore:
    """Progress counters and partial PNG payloads for decode / render pipelines."""

    def __init__(
        self,
        *,
        progress_ttl_s: float = 600.0,
        partial_ttl_s: float = 600.0,
    ) -> None:
        self._progress_ttl_s = progress_ttl_s
        self._partial_ttl_s = partial_ttl_s
        self._progress_lock = threading.Lock()
        self._partial_lock = threading.Lock()
        self._progress: dict[str, dict[str, object]] = {}
        self._partial: dict[str, dict[str, object]] = {}

    def progress_register(
        self,
        token: str,
        total: int,
        workers: int,
        phase: str,
        *,
        rerender: bool = False,
    ) -> None:
        with self._progress_lock:
            self._progress[token] = {
                "total": max(0, int(total)),
                "done": 0,
                "workers": max(1, int(workers)),
                "phase": str(phase),
                "rerender": bool(rerender),
                "t0": time.monotonic(),
            }

    def progress_set_done(self, token: str, done: int) -> None:
        with self._progress_lock:
            entry = self._progress.get(token)
            if not entry:
                return
            total = int(entry["total"])
            entry["done"] = max(0, min(total, int(done)))

    def progress_snapshot(self, token: str) -> dict[str, object] | None:
        with self._progress_lock:
            entry = self._progress.get(token)
            if not entry:
                return None
            if time.monotonic() - float(entry["t0"]) > self._progress_ttl_s:
                self._progress.pop(token, None)
                return None
            total = int(entry["total"])
            done = int(entry["done"])
            pct = (100.0 * done / total) if total else 0.0
            phase = str(entry.get("phase") or "decode")
            workers = int(entry["workers"])
            display_total = entry.get("display_total")
            snap_total = int(display_total) if display_total is not None else total
            snap: dict[str, object] = {
                "total": snap_total,
                "done": done,
                "workers": workers,
                "percent": round(pct, 1),
                "phase": phase,
                "rerender": bool(entry.get("rerender")),
            }
            if phase == "decode":
                snap["n_gpus"] = workers
            elif phase == "pipeline":
                snap["n_gpus"] = int(entry.get("n_gpus") or workers)
                snap["n_cpus"] = int(entry.get("n_cpus") or 1)
                decode_done = int(entry.get("decode_done") or 0)
                render_done = int(entry.get("render_done") or done)
                if snap_total > 0:
                    decode_done = min(snap_total, decode_done)
                    render_done = min(snap_total, render_done)
                snap["decode_done"] = decode_done
                snap["render_done"] = render_done
            else:
                snap["n_cpus"] = workers
            return snap

    def progress_unregister(self, token: str) -> None:
        with self._progress_lock:
            self._progress.pop(token, None)

    def progress_register_pipeline(
        self,
        token: str,
        total: int,
        *,
        n_gpus: int,
        n_cpus: int,
        display_total: int | None = None,
    ) -> None:
        with self._progress_lock:
            internal_total = max(0, int(total))
            progress_total = (
                max(0, int(display_total))
                if display_total is not None
                else internal_total
            )
            entry: dict[str, object] = {
                "total": progress_total,
                "done": 0,
                "workers": max(1, int(n_gpus)),
                "phase": "pipeline",
                "rerender": False,
                "n_gpus": max(1, int(n_gpus)),
                "n_cpus": max(1, int(n_cpus)),
                "decode_done": 0,
                "render_done": 0,
                "t0": time.monotonic(),
            }
            if display_total is not None:
                entry["display_total"] = progress_total
                entry["internal_total"] = internal_total
            self._progress[token] = entry

    def progress_update_pipeline(
        self,
        token: str,
        *,
        decode_done: int | None = None,
        render_done: int | None = None,
    ) -> None:
        with self._progress_lock:
            entry = self._progress.get(token)
            if not entry:
                return
            total = int(entry["total"])
            if decode_done is not None:
                entry["decode_done"] = max(0, min(total, int(decode_done)))
            if render_done is not None:
                rd = max(0, min(total, int(render_done)))
                entry["render_done"] = rd
                entry["done"] = rd
            entry["phase"] = "pipeline"

    def partial_register(self, token: str, total: int) -> None:
        with self._partial_lock:
            self._partial[token] = {
                "total": max(0, int(total)),
                "decode_done": 0,
                "render_done": 0,
                "images": {},
                "view_matrix": None,
                "complete": False,
                "error": None,
                "t0": time.monotonic(),
            }

    def partial_set_image(self, token: str, index: int, png_bytes: bytes) -> None:
        b64 = base64.standard_b64encode(png_bytes).decode("ascii")
        with self._partial_lock:
            entry = self._partial.get(token)
            if not entry:
                return
            images = entry.get("images")
            if not isinstance(images, dict):
                images = {}
                entry["images"] = images
            images[int(index)] = b64
            entry["render_done"] = len(images)

    def partial_update_decode(self, token: str, decode_done: int) -> None:
        with self._partial_lock:
            entry = self._partial.get(token)
            if not entry:
                return
            total = int(entry["total"])
            entry["decode_done"] = max(0, min(total, int(decode_done)))

    def partial_mark_complete(self, token: str, view_matrix: str | None = None) -> None:
        with self._partial_lock:
            entry = self._partial.get(token)
            if not entry:
                return
            entry["complete"] = True
            if view_matrix:
                entry["view_matrix"] = view_matrix

    def partial_snapshot(self, token: str) -> dict[str, object] | None:
        with self._partial_lock:
            entry = self._partial.get(token)
            if not entry:
                return None
            if time.monotonic() - float(entry["t0"]) > self._partial_ttl_s:
                self._partial.pop(token, None)
                return None
            total = int(entry["total"])
            images_raw = entry.get("images")
            images_list: list[dict[str, object]] = []
            if isinstance(images_raw, dict):
                for idx in sorted(images_raw):
                    images_list.append({"index": int(idx), "b64": str(images_raw[idx])})
            return {
                "total": total,
                "decode_done": int(entry.get("decode_done") or 0),
                "render_done": int(entry.get("render_done") or 0),
                "complete": bool(entry.get("complete")),
                "images": images_list,
                "view_matrix": entry.get("view_matrix"),
            }

    def partial_unregister(self, token: str) -> None:
        with self._partial_lock:
            self._partial.pop(token, None)

    def progress_entry(self, token: str) -> dict[str, object] | None:
        """Mutable progress entry for in-pipeline updates (caller holds no lock)."""
        with self._progress_lock:
            return self._progress.get(token)
