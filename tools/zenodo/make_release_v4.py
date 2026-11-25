"""Build the data archives of release v4.0.0 (parts a and b) for the Zenodo record, in release/zenodo_v4/ (ignored).

    python tools/zenodo/make_release_v4.py --runs RUNS --fits FITS --processed PROC [--out release/zenodo_v4]
                                           [--only ARCHIVE ...]

RUNS       a directory with one entry per run of data/provenance/v4/run_map.csv, named by its run-map name (the run
           directory or a link to it); only the runs of the archives v4_runs_a and v4_runs_b are read.
FITS       the model B fit directories b1/, b2/ and h1_fc2/.
PROC       the output directory of scripts/v4/regenerate_v4.py (PROC/a, PROC/b).

Archives (tar.gz; members sorted, fixed mtime, uid/gid 0, no owner names, so that a rebuild gives the same bytes):
  v4_runs_a, v4_runs_b   the run directories of the run map's archive column, top directory = run name. Every file is
                         checked against data/provenance/v4/manifest_v4.json (sha256 as produced). Text files are
                         scrubbed: the home directory of the machine that wrote the file (HOME_PATH) becomes
                         "<home>". A file changed by the scrub must match "sha256_archived" of the manifest; no other
                         file may change. The input files of the ALAMODE SrTiO3 example (STO_anharm.xml, STO222.xml,
                         BORN, identified by their sha256) are not included; they are distributed with ALAMODE
                         (commit c5b015c).
  v4_fits_b              FITS/b1, FITS/b2, FITS/h1_fc2 under v4_fits_b/.
  v4_processed_a, _b     the files of PROC/a and PROC/b that the repository does not hold under data/processed/v4/a, /b
                         (per-mode outputs: family labels, census tables, mode tables), under v4_processed_{a,b}/.
Every archive holds FILES.sha256 (sha256 of each member as archived). The output directory gets SHA256SUMS (the
archives) and a copy of manifest_v4.json. Nothing of the files of the 2015 study is read or written.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import shutil
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MTIME = 1760000000                      # fixed member mtime, as for the v3 archives
LIMIT = 4 * 1024 ** 3                   # 4 GB per file
HOME_PATH = re.compile(rb"/(?:Users|home)/[^/\s<>\"']+")
EXAMPLE_SHA256 = {                      # ALAMODE c5b015c, example/SrTiO3: not redistributed
    "STO_anharm.xml": "5e787b95eb59fdc2a5cd93a3fe298d976f71d3a83905dce568cc021d86e71b79",
    "STO222.xml": "9739c36af71eb1ab306760e0a86a7ba09f7337c6220901272fea1b20400f9f0f",
    "BORN": "e283fbf4a6d75c7dbe00a8d8589c1a0a241ea9a4255bc1fc1b078e2c229a1bb9",
}
RUN_ARCHIVES = ("v4_runs_a", "v4_runs_b")
ALL = RUN_ARCHIVES + ("v4_fits_b", "v4_processed_a", "v4_processed_b")


def sha256(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def is_text(b: bytes) -> bool:
    return b"\0" not in b[:1 << 20]


def scrub(b: bytes) -> bytes:
    """Home path of the machine that wrote the file -> "<home>" (text files only)."""
    return HOME_PATH.sub(b"<home>", b) if is_text(b) else b


def write_archive(out: Path, name: str, members) -> tuple[Path, int]:
    """members: iterable of (archive path, bytes), in any order; written sorted by path."""
    members = sorted(members, key=lambda m: m[0])
    listing = "".join(f"{sha256(b)}  {p}\n" for p, b in members).encode()
    members.append((f"{name}/FILES.sha256", listing))
    tar_path, final = out / f"{name}.tar", out / f"{name}.tar.gz"
    with open(tar_path, "wb") as fh, tarfile.open(fileobj=fh, mode="w", format=tarfile.PAX_FORMAT) as tf:
        for arc, data in members:
            ti = tarfile.TarInfo(arc)
            ti.size, ti.mtime, ti.mode, ti.uid, ti.gid, ti.uname, ti.gname = len(data), MTIME, 0o644, 0, 0, "", ""
            tf.addfile(ti, io.BytesIO(data))
    with open(tar_path, "rb") as src, open(final, "wb") as dst, \
            gzip.GzipFile(filename="", mode="wb", fileobj=dst, mtime=0, compresslevel=9) as g:
        shutil.copyfileobj(src, g, 1 << 20)
    tar_path.unlink()
    if final.stat().st_size > LIMIT:
        sys.exit(f"stop: {final.name} exceeds 4 GB")
    return final, len(members)


def run_members(runs_dir: Path, manifest: dict, archive: str, log: list):
    for run in manifest["runs"]:
        if run["archive"] != archive:
            continue
        d = runs_dir / run["name"]
        if not d.is_dir():
            sys.exit(f"stop: run directory {d} missing")
        on_disk = sorted(str(p.relative_to(d)) for p in d.rglob("*") if p.is_file() or p.is_symlink())
        listed = sorted(f["file"] for f in run["files"])
        if on_disk != listed:
            sys.exit(f"stop: files of {run['name']} differ from the manifest")
        for f in run["files"]:
            if f.get("link"):
                log.append(("link not archived", run["name"], f["file"]))
                continue
            b = (d / f["file"]).read_bytes()
            if sha256(b) != f["sha256"]:
                sys.exit(f"stop: sha256 of {run['name']}/{f['file']} differs from the manifest")
            if EXAMPLE_SHA256.get(Path(f["file"]).name) == f["sha256"]:
                log.append(("example input not archived", run["name"], f["file"]))
                continue
            s = scrub(b)
            if s != b:
                if sha256(s) != f.get("sha256_archived"):
                    sys.exit(f"stop: scrubbed {run['name']}/{f['file']} differs from sha256_archived")
                log.append(("scrubbed", run["name"], f["file"]))
            elif "sha256_archived" in f:
                sys.exit(f"stop: {run['name']}/{f['file']} has sha256_archived but nothing to scrub")
            yield f"{run['name']}/{f['file']}", s


def dir_members(top: str, base: Path, rels, log: list):
    for rel in rels:
        b = (base / rel).read_bytes()
        s = scrub(b)
        if s != b:
            sys.exit(f"stop: {base / rel} contains a home path")
        yield f"{top}/{rel}", b


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", type=Path, required=True)
    ap.add_argument("--fits", type=Path, required=True)
    ap.add_argument("--processed", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=ROOT / "release/zenodo_v4")
    ap.add_argument("--only", nargs="+", choices=ALL)
    a = ap.parse_args(argv)
    manifest_path = ROOT / "data/provenance/v4/manifest_v4.json"
    manifest = json.loads(manifest_path.read_text())
    with open(ROOT / "data/provenance/v4/run_map.csv", newline="") as fh:
        rm = {r["run"]: r["archive"] for r in csv.DictReader(fh)}
    if any(rm.get(r["name"]) != r["archive"] for r in manifest["runs"]):
        sys.exit("stop: run map and manifest disagree")
    a.out.mkdir(parents=True, exist_ok=True)
    built, log = [], []
    for name in a.only or ALL:
        if name in RUN_ARCHIVES:
            members = list(run_members(a.runs, manifest, name, log))
        elif name == "v4_fits_b":
            rels = sorted(str(p.relative_to(a.fits)) for sub in ("b1", "b2", "h1_fc2")
                          for p in (a.fits / sub).rglob("*") if p.is_file())
            members = list(dir_members(name, a.fits, rels, log))
        else:
            part = name[-1]
            src, held = a.processed / part, ROOT / "data/processed/v4" / part
            rels = sorted(str(p.relative_to(src)) for p in src.rglob("*")
                          if p.is_file() and not (held / p.relative_to(src)).exists())
            members = list(dir_members(name, src, rels, log))
        path, n = write_archive(a.out, name, members)
        built.append(path)
        print(f"{name}: {n} members, {path.stat().st_size:,} bytes")
    shutil.copyfile(manifest_path, a.out / "manifest_v4.json")
    sums = sorted(set(built) | set(a.out.glob("v4_*.tar.gz")))
    (a.out / "SHA256SUMS").write_text("".join(
        f"{sha256_file(p)}  {p.name}\n" for p in sums + [a.out / "manifest_v4.json"]))
    for kind in sorted({k for k, _, _ in log}):
        print(f"{kind}: {sum(1 for k, _, _ in log if k == kind)}")


if __name__ == "__main__":
    main()
