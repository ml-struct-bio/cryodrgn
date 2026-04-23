"""Extract argparse ``help=`` strings from command modules without importing them.

``abinit`` and training commands import ``torch`` at module level; parsing the
``.py`` files with :mod:`ast` keeps the dashboard importable in minimal envs.
"""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Any


def _string_from_ast(node: ast.expr | None) -> str | None:
    if node is None:
        return None
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        left = _string_from_ast(node.left)
        right = _string_from_ast(node.right)
        if left is not None and right is not None:
            return left + right
        return None
    if isinstance(node, ast.JoinedStr):
        parts: list[str] = []
        for v in node.values:
            if isinstance(v, ast.Constant) and isinstance(v.value, str):
                parts.append(v.value)
            else:
                return None
        return "".join(parts)
    return None


def _flags_and_positional_from_add_argument(
    call: ast.Call,
) -> tuple[list[str], str | None]:
    """Return (option_strings, first_positional_name).

    Positional names (e.g. ``particles``) are stored in the map for required
    fields that are not introduced with ``-``.
    """
    opts: list[str] = []
    first_plain: str | None = None
    for arg in call.args:
        s = _string_from_ast(arg)
        if s is None:
            continue
        if s.startswith("-"):
            opts.append(s)
        elif first_plain is None:
            first_plain = s
    return opts, first_plain


def _help_from_add_argument(call: ast.Call) -> str | None:
    for kw in call.keywords:
        if kw.arg == "help":
            return _string_from_ast(kw.value)
    return None


def help_map_from_command_py(path: Path) -> dict[str, str]:
    """Map each option string (e.g. ``--load``, ``-n``) to its help text."""
    text = path.read_text(encoding="utf-8")
    tree = ast.parse(text, filename=str(path))
    m: dict[str, str] = {}
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        if not isinstance(fn, ast.Attribute) or fn.attr != "add_argument":
            continue
        opts, pos_name = _flags_and_positional_from_add_argument(node)
        h = _help_from_add_argument(node)
        if not h or (not opts and not pos_name):
            continue
        h_norm = " ".join(h.split())
        for opt in opts:
            m.setdefault(opt, h_norm)
        if pos_name:
            m.setdefault(pos_name, h_norm)
    return m


def attach_help_to_groups(
    help_map: dict[str, str],
    groups: list[dict[str, Any]],
) -> None:
    """Set ``help`` on each arg dict in *groups* when a CLI flag matches *help_map*."""
    for g in groups:
        for a in g["args"]:
            w = a.get("w")
            if w == "no_amp":
                t = help_map.get("--no-amp")
                if t:
                    a["help"] = t
                continue
            for c in a.get("cli") or []:
                t = help_map.get(c)
                if t:
                    a["help"] = t
                    break


def load_cli_help_maps() -> dict[str, dict[str, str]]:
    cmd_dir = Path(__file__).resolve().parent.parent / "commands"
    out: dict[str, dict[str, str]] = {}
    for name in ("abinit", "train_vae", "train_nn", "train_dec"):
        p = cmd_dir / f"{name}.py"
        if p.is_file():
            out[name] = help_map_from_command_py(p)
    return out
