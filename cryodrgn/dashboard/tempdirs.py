"""Temporary folders of a dashboard process.

Decoded volumes are written to disk so they can be re-rendered, turned into GIFs
or saved later, and stay there for up to two hours. Every such folder is created
under one per-process root, which is removed when the process exits; roots left
behind by a process that was killed are removed the next time a dashboard starts
on the same machine.
"""

from __future__ import annotations

import atexit
import logging
import os
import re
import shutil
import socket
import tempfile
import threading
import time

logger = logging.getLogger(__name__)

_ROOT_PREFIX = "cryodrgn_dashboard_"
# ``cryodrgn_dashboard_<host>_<pid>_<random>``
_ROOT_RE = re.compile(r"^cryodrgn_dashboard_(?P<host>[A-Za-z0-9.-]+)_(?P<pid>\d+)_")

# Releases up to 4.3.2 created these directly under the temp directory and never
# removed them. Their cache entries expire after two hours, so day-old ones are dead.
_LEGACY_PREFIXES = ("cryodrgn_trajectory_mrc_", "cryodrgn_explorer_mrc_")
_LEGACY_MIN_AGE_S = 24 * 3600.0

_lock = threading.Lock()
_root: str | None = None
_root_pid: int | None = None
_atexit_registered = False


def _host_tag() -> str:
    return re.sub(r"[^A-Za-z0-9.-]", "-", socket.gethostname() or "host")


def session_root() -> str:
    """Folder that holds every temporary folder of this process (created on demand)."""
    global _root, _root_pid, _atexit_registered
    with _lock:
        pid = os.getpid()
        if _root is None or _root_pid != pid or not os.path.isdir(_root):
            _root = tempfile.mkdtemp(prefix=f"{_ROOT_PREFIX}{_host_tag()}_{pid}_")
            _root_pid = pid
            if not _atexit_registered:
                atexit.register(remove_session_root)
                _atexit_registered = True

        return _root


def make_tempdir(prefix: str) -> str:
    """Like :func:`tempfile.mkdtemp`, inside :func:`session_root`."""
    return tempfile.mkdtemp(prefix=prefix, dir=session_root())


def remove_session_root() -> None:
    """Delete this process's temporary folders (registered with :mod:`atexit`)."""
    global _root, _root_pid
    with _lock:
        root, pid = _root, _root_pid
        # A forked child inherits the parent's root: only its creator removes it.
        if root is None or pid != os.getpid():
            return
        _root = _root_pid = None

    shutil.rmtree(root, ignore_errors=True)


def _pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return True

    return True


def _newest_mtime(path: str) -> float:
    newest = os.stat(path).st_mtime
    with os.scandir(path) as it:
        for entry in it:
            try:
                newest = max(newest, entry.stat(follow_symlinks=False).st_mtime)
            except OSError:
                continue

    return newest


def sweep_stale_tempdirs(tmp_dir: str | None = None) -> list[str]:
    """Remove dashboard folders this user left in the temp directory; return them.

    Two kinds are removed:

    * roots of dashboard processes on this host that no longer exist;
    * volume folders written by earlier releases, once they are a day old.

    Folders of other users, of other hosts (a shared ``TMPDIR``) and of running
    dashboards are left alone.
    """

    tmp_dir = tmp_dir or tempfile.gettempdir()
    uid = os.getuid()
    host = _host_tag()
    now = time.time()
    removed: list[str] = []

    try:
        entries = list(os.scandir(tmp_dir))
    except OSError:
        return removed

    for entry in entries:
        name = entry.name
        try:
            if not entry.is_dir(follow_symlinks=False):
                continue
            if entry.stat(follow_symlinks=False).st_uid != uid:
                continue
            m = _ROOT_RE.match(name)
            if m:
                stale = m.group("host") == host and not _pid_alive(int(m.group("pid")))
            elif name.startswith(_LEGACY_PREFIXES):
                stale = now - _newest_mtime(entry.path) > _LEGACY_MIN_AGE_S
            else:
                continue
        except OSError:
            continue
        if stale:
            shutil.rmtree(entry.path, ignore_errors=True)
            if not os.path.exists(entry.path):
                removed.append(entry.path)

    if removed:
        logger.info(
            "Removed %d leftover dashboard temp folder(s) under %s",
            len(removed),
            tmp_dir,
        )

    return removed
