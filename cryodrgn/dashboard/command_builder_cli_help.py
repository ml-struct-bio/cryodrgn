"""Extract argparse ``help=`` strings from command modules without importing them.

``abinit`` and training commands import ``torch`` at module level; parsing the
``.py`` files with :mod:`ast` keeps the dashboard importable in minimal envs.
"""

from __future__ import annotations

import ast
import re
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


def _literal_from_ast(node: ast.expr | None) -> Any:
    if node is None:
        return None
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        inner = _literal_from_ast(node.operand)
        if isinstance(inner, (int, float)):
            return -inner
    return None


def _default_from_add_argument(call: ast.Call) -> Any:
    for kw in call.keywords:
        if kw.arg == "default":
            return _literal_from_ast(kw.value)
    return None


def resolve_help_default_placeholders(help: str, default: Any) -> str:
    """Replace argparse ``%(default)s`` with a concrete value for UI tooltips."""
    if "%(default)s" not in help:
        return help
    if default is None:
        disp = "None"
    elif isinstance(default, bool):
        disp = str(default)
    else:
        disp = str(default)
    return help.replace("%(default)s", disp)


def _default_for_arg(arg: dict[str, Any], defaults_map: dict[str, Any]) -> Any:
    for c in arg.get("cli") or []:
        if c in defaults_map:
            return defaults_map[c]
    if "placeholder" in arg:
        return arg.get("placeholder")
    if "default" in arg:
        return arg.get("default")
    return None


def defaults_map_from_command_py(path: Path) -> dict[str, Any]:
    """Map each CLI flag to its argparse ``default=`` when present in the source."""
    text = path.read_text(encoding="utf-8")
    tree = ast.parse(text, filename=str(path))
    m: dict[str, Any] = {}
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        if not isinstance(fn, ast.Attribute) or fn.attr != "add_argument":
            continue
        opts, pos_name = _flags_and_positional_from_add_argument(node)
        default = _default_from_add_argument(node)
        if default is None and not opts and not pos_name:
            continue
        for opt in opts:
            if default is not None:
                m.setdefault(opt, default)
        if pos_name is not None and default is not None:
            m.setdefault(pos_name, default)
    return m


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


def display_name_from_arg_help(help: str | None) -> str | None:
    """Short human label for a flag, derived from argparse ``help=`` text."""
    if not help or not str(help).strip():
        return None
    s = " ".join(str(help).split())
    s = re.sub(r"\s*\(default:.*$", "", s, flags=re.I).strip()
    s = s.rstrip(".")
    if not s:
        return None

    def lower(x: str) -> str:
        return x.strip().lower()

    m = re.match(r"^flag for (.+)$", s, re.I)
    if m:
        return lower(m.group(1))

    m = re.match(r"^input (.+)$", s, re.I)
    if m:
        inner = lower(re.split(r"[.(]", m.group(1), maxsplit=1)[0])
        return inner if inner else None

    m = re.match(r"^dimension of (.+)$", s, re.I)
    if m:
        phrase = lower(re.sub(r";.*$", "", m.group(1)))
        phrase = re.sub(r"\s+required\.?$", "", phrase).strip()
        return phrase if phrase and len(phrase) <= 42 else None

    m = re.match(r"^number of total epochs(?: to train for)?$", s, re.I)
    if m:
        return "total epochs"

    m = re.match(r"^number of (.+)$", s, re.I)
    if m:
        phrase = lower(m.group(1))
        phrase = re.sub(r" to train for$", "", phrase)
        phrase = re.sub(r" \(.*$", "", phrase)
        if phrase and len(phrase) <= 42:
            return phrase

    m = re.match(r"^path to (?:a |an |the )?(.+)$", s, re.I)
    if m:
        inner = lower(m.group(1))
        inner = re.split(r"[.(]", inner, maxsplit=1)[0].strip()
        if inner.endswith(" to load"):
            inner = inner[: -len(" to load")].strip()
        if inner and len(inner) <= 38:
            return inner
        return None

    m = re.match(r"^training batch size used for (.+)$", s, re.I)
    if m:
        return lower(f"batch size ({m.group(1)})")

    m = re.match(r"^(.+?) for the (.+?)(?: optimizer)?$", s, re.I)
    if m:
        head = lower(m.group(1))
        tail = re.sub(r"\s+optimizer$", "", lower(m.group(2)), flags=re.I)
        if head in ("learning rate", "optimizer") and tail:
            label = f"{head} ({tail})"
            return label if len(label) <= 42 else head
        if len(head) <= 40:
            return head

    for sep in (" for a ", " for an "):
        if sep in s.lower():
            part = re.split(re.escape(sep), s, maxsplit=1, flags=re.I)[0].strip()
            if 3 <= len(part) <= 40:
                return lower(part)

    if re.search(r"\s+used for\s+", s, re.I):
        part = re.split(r"\s+used for\s+", s, maxsplit=1, flags=re.I)[0].strip()
        if 3 <= len(part) <= 42:
            return lower(part)

    skip = (
        "disable ",
        "activate ",
        "do not ",
        "bypass ",
        "increase ",
        "print ",
        "fix the ",
        "indicate that ",
        "when using",
    )
    if any(s.lower().startswith(p) for p in skip):
        return None

    first = re.split(r"[;,]", s)[0].strip()
    if len(first) <= 38:
        return lower(first)
    return None


_DISPLAY_NAME_OVERRIDES: dict[str, str] = {
    "--ind": "filtering .pkl",
}


def jinja_arg_display_name(source: dict[str, Any] | str | None) -> str | None:
    """Jinja filter: accept an arg dict (``help`` key) or a help string."""
    if isinstance(source, dict):
        cli = source.get("cli")
        if isinstance(cli, (list, tuple)):
            for flag in cli:
                override = _DISPLAY_NAME_OVERRIDES.get(str(flag))
                if override:
                    return override
        help_text = source.get("help")
    elif isinstance(source, str):
        help_text = source
    else:
        return None
    return display_name_from_arg_help(help_text)


def attach_help_to_groups(
    help_map: dict[str, str],
    groups: list[dict[str, Any]],
    *,
    defaults_map: dict[str, Any] | None = None,
) -> None:
    """Set ``help`` on each arg dict in *groups* when a CLI flag matches *help_map*."""
    dm = defaults_map or {}
    for g in groups:
        for a in g["args"]:
            w = a.get("w")
            if w == "no_amp":
                t = help_map.get("--no-amp")
                if t:
                    a["help"] = resolve_help_default_placeholders(
                        t, _default_for_arg(a, dm)
                    )
                continue
            for c in a.get("cli") or []:
                t = help_map.get(c)
                if t:
                    a["help"] = resolve_help_default_placeholders(
                        t, dm.get(c, _default_for_arg(a, dm))
                    )
                    break


def _command_module_paths() -> dict[str, Path]:
    cmd_dir = Path(__file__).resolve().parent.parent / "commands"
    names = (
        "abinit",
        "abinit_het_old",
        "abinit_homo_old",
        "train_vae",
        "train_nn",
        "train_dec",
        "backproject_voxel",
        "analyze",
        "analyze_landscape",
        "analyze_landscape_full",
    )
    return {name: cmd_dir / f"{name}.py" for name in names}


def load_cli_help_maps() -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for name, p in _command_module_paths().items():
        if p.is_file():
            out[name] = help_map_from_command_py(p)
    return out


def load_cli_defaults_maps() -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for name, p in _command_module_paths().items():
        if p.is_file():
            out[name] = defaults_map_from_command_py(p)
    return out


def resolved_help_for_flag(
    help_map: dict[str, str],
    defaults_map: dict[str, Any],
    flag: str,
) -> str | None:
    """Help text for a CLI flag with ``%(default)s`` filled in."""
    t = help_map.get(flag)
    if not t:
        return None
    return resolve_help_default_placeholders(t, defaults_map.get(flag))


def module_summary_from_command_py(path: Path) -> str:
    """First paragraph of a command module docstring (before example blocks)."""
    text = path.read_text(encoding="utf-8")
    tree = ast.parse(text, filename=str(path))
    doc = ast.get_docstring(tree) or ""
    lines: list[str] = []
    for line in doc.splitlines():
        stripped = line.strip()
        if stripped.startswith("Example usage") or stripped.startswith("---"):
            break
        if stripped:
            lines.append(stripped)
        elif lines:
            break
    if lines:
        return " ".join(lines)
    first = doc.strip().splitlines()
    return first[0].strip() if first else ""


def load_command_module_docstrings() -> dict[str, str]:
    """Map command keys to one-line module summaries for the command builder."""
    cmd_dir = Path(__file__).resolve().parent.parent / "commands"
    out: dict[str, str] = {}
    for name in (
        "abinit",
        "abinit_het_old",
        "abinit_homo_old",
        "train_vae",
        "train_nn",
        "train_dec",
        "backproject_voxel",
        "analyze",
        "analyze_landscape",
        "analyze_landscape_full",
    ):
        p = cmd_dir / f"{name}.py"
        if p.is_file():
            out[name] = module_summary_from_command_py(p)
    return out
