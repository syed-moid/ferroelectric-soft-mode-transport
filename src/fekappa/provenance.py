"""Provenance headers for the derived tables, and the locations of inputs kept outside the repository.

Every table written by `scripts/revision/` starts with three comment lines:

    # script: scripts/revision/<name>.py
    # inputs: <repo-relative paths, or zenodo:<archive>/<path>, or $FEKAPPA_RUNS/<dir>>
    # generated: <UTC date>; fekappa <version>; manifest_final.json sha256 <12 hex>

followed by any descriptive comment lines the script writes itself. All readers use `comment="#"`.
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import os
import sys
from pathlib import Path

from fekappa import __version__

_KEYS = ("# script:", "# inputs:", "# generated:")


def _repo_root(script: Path) -> Path:
    p = Path(script).resolve()
    for d in p.parents:
        if (d / "pyproject.toml").exists():
            return d
    return p.parent


def manifest_sha12(root: Path) -> str:
    m = root / "data/provenance/manifest_final.json"
    return hashlib.sha256(m.read_bytes()).hexdigest()[:12] if m.exists() else "absent"


def header_lines(script, inputs) -> list[str]:
    root = _repo_root(Path(script))
    try:
        rel = Path(script).resolve().relative_to(root)
    except ValueError:
        rel = Path(script).name
    date = _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%d")
    return [f"# script: {rel}",
            f"# inputs: {'; '.join(inputs)}",
            f"# generated: {date}; fekappa {__version__}; manifest_final.json sha256 {manifest_sha12(root)}"]


def stamp(path, script, inputs) -> None:
    """Prepend the provenance header to a table the script has just written (replacing an earlier one); keep the
    script's own comment lines after it. Works on plain and gzipped CSV files."""
    import gzip
    path = Path(path)
    gz = path.suffix == ".gz"
    text = (gzip.decompress(path.read_bytes()).decode() if gz else path.read_text())
    lines = text.splitlines(keepends=True)
    i = 0
    while i < len(lines) and lines[i].startswith("#"):
        i += 1
    own = [ln for ln in lines[:i] if not ln.startswith(_KEYS)]
    out = "".join(ln + "\n" for ln in header_lines(script, inputs)) + "".join(own) + "".join(lines[i:])
    if gz:
        path.write_bytes(gzip.compress(out.encode(), mtime=0))
    else:
        path.write_text(out)


def external_dir(var: str, what: str) -> Path:
    """Directory of inputs that live outside the repository (pre-queue runs, the DFT dataset), from an environment
    variable; exits with a message naming what is needed when it is unset."""
    v = os.environ.get(var)
    if not v:
        sys.exit(f"set {var} to {what}")
    return Path(v)
