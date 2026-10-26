"""Build the complete Zenodo upload for release v3.0.0 in release/zenodo_v3/ (gitignored).

One Zenodo version holds code and data together, so the paper cites one version DOI:
    ferroelectric-soft-mode-transport-v3.0.0.tar.gz   git archive of tag v3.0.0 (code, tables, patches, records, docs)
    <variant>.tar.gz                                  raw outputs of the queue jobs, one archive per variant
                                                      (make_zenodo_package.py: checked against the manifest sha256s)
    prequeue_runs.tar.gz                              the runs made before the queue existed ($FEKAPPA_RUNS)
    manifest_final.json                               scrubbed job manifest (identical to data/provenance/ at the tag)
    software_versions.json                            ALAMODE commit, patches, full binary hashes, package lock
    README.md                                         layout, software, how to regenerate every table and figure
    zenodo_metadata.json                              metadata to enter on the Zenodo upload form
    SHA256SUMS                                        sha256 of every file above

Host and account strings are replaced by neutral tokens in every text file (as for the manifest). The force-constant
inputs of the ALAMODE example (STO_anharm.xml, STO222.xml, BORN) and symbolic links are not redistributed.

    python tools/zenodo/make_release_v3.py [--toolchain DIR] [--skip-raw]
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
import subprocess
import sys
import tarfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_zenodo_package as zp  # noqa: E402

ROOT = zp.ROOT
TAG = "v3.0.0"
CODE = f"ferroelectric-soft-mode-transport-{TAG}.tar.gz"
LIMIT = 4 * 1024 ** 3  # 4 GB per file
PREQUEUE_SKIP_DIRS = {"queue", "repro300"}  # queue/: work copies of archived jobs; repro300: superseded test
EXAMPLE_INPUTS = {"STO_anharm.xml", "STO222.xml", "BORN"}
CONCEPT_DOI = "10.5281/zenodo.20559513"
REPO = "https://github.com/syed-moid/ferroelectric-soft-mode-transport"


def gz_tar(out: Path, name: str, members: list[tuple[str, bytes]]) -> Path:
    tar_path, final = out / f"{name}.tar", out / f"{name}.tar.gz"
    with open(tar_path, "wb") as fh, tarfile.open(fileobj=fh, mode="w", format=tarfile.PAX_FORMAT) as tf:
        for arc, data in members:
            ti = tarfile.TarInfo(arc)
            ti.size, ti.mtime, ti.mode, ti.uid, ti.gid, ti.uname, ti.gname = len(data), zp.MTIME, 0o644, 0, 0, "", ""
            tf.addfile(ti, io.BytesIO(data))
    with open(tar_path, "rb") as src, open(final, "wb") as dst, \
            gzip.GzipFile(filename="", mode="wb", fileobj=dst, mtime=0, compresslevel=9) as g:
        shutil.copyfileobj(src, g, 1 << 20)
    tar_path.unlink()
    return final


def prequeue(out: Path, runs: Path, rules) -> dict:
    members, skipped, n_sub = [], {"symlinks": 0, "example inputs": 0}, 0
    for p in sorted(runs.rglob("*")):
        rel = p.relative_to(runs)
        if rel.parts[0] in PREQUEUE_SKIP_DIRS or p.is_dir():
            continue
        if p.is_symlink():
            skipped["symlinks"] += 1
            continue
        if p.name in EXAMPLE_INPUTS:
            skipped["example inputs"] += 1
            continue
        data, k = zp.scrub_bytes(p.read_bytes(), rules)
        n_sub += k
        members.append((f"prequeue_runs/{rel}", data))
    f = gz_tar(out, "prequeue_runs", members)
    return {"file": f.name, "files": len(members), "uncompressed_bytes": sum(len(d) for _, d in members),
            "bytes": f.stat().st_size, "substitutions": n_sub, "skipped": skipped,
            "directories": sorted({m[0].split("/")[1] for m in members if m[0].count("/") > 1})}


def builds() -> list[dict]:
    """Full binary hashes from tools/alamode_patches/BUILD.md at the tag, joined with the build manifest."""
    md = subprocess.run(["git", "show", f"{TAG}:tools/alamode_patches/BUILD.md"], cwd=ROOT, capture_output=True,
                        text=True, check=True).stdout
    rows = []
    for ln in md.splitlines():
        m = re.match(r"\| (.+?) \| `([^`]+)` \| (.+?) \| (HOME|NODE) \| `([0-9a-f]{64})` \| (.+?) \|$", ln)
        if m:
            patch = re.match(r"`([^`]+)` \(`([0-9a-f]{64})`\)", m.group(3))
            rows.append({"variant": m.group(1), "tree": m.group(2), "machine": m.group(4), "anphon_sha256": m.group(5),
                         "patch": patch.group(1) if patch else None, "patch_sha256": patch.group(2) if patch else None,
                         "jobs_or_use": m.group(6)})
    return rows


def readme(files: list[dict], sw: dict) -> str:
    rows = "\n".join(f"| `{f['file']}` | {f['bytes'] / 1e6:,.1f} | {f['what']} |" for f in files)
    brows = "\n".join(f"| {b['variant']} | `{b['tree']}` | {b['patch'] or '—'} | {b['machine']} | `{b['anphon_sha256']}` |"
                      for b in sw["builds"])
    return f"""# Zone-resolved soft-mode damping and the thermal-conductivity decomposition of SrTiO₃: release {TAG}

Code, derived tables and raw calculation outputs of the revised manuscript (J. Phys.: Condens. Matter, submission
JPCM-126876), archived together as one version of the record {CONCEPT_DOI} (all versions). The code archive is the
tagged release {REPO}/releases/tag/{TAG}.

## Files

| file | MB | content |
|---|---|---|
{rows}

Checksums: `shasum -a 256 -c SHA256SUMS`. Unpack: `tar -xzf <file>.tar.gz`.

## How the parts fit

1. Unpack the code archive; install with Python 3.11: `pip install -e ".[dev]"` (or `uv sync --frozen --extra dev`).
2. From tracked tables alone (no raw data needed): `pytest`; `python scripts/revision/build_manuscript_numbers_v8.py`
   regenerates every number of the main text; `cd paper/figures && python make_v8_figures.py && python
   make_si_v4_figures.py` regenerates Figs. 1–4 and S1–S8.
3. To regenerate the tables from the raw outputs: unpack the variant archives inside the code tree with
   `mkdir -p data/raw/alamode_sto/revision && for f in pm_q12_rta pm_q12_evec scph_q12 scph_q16 pm_q8 pm_q12_4ph_audit;
   do tar -xzf $f.tar.gz --strip-components=1 -C data/raw/alamode_sto/revision; done`; unpack `prequeue_runs.tar.gz`
   anywhere and set `FEKAPPA_RUNS` to its `prequeue_runs/` directory; then run the scripts of `scripts/revision/` in the
   order of its README. `FEKAPPA_TOOLCHAIN` (the ALAMODE example inputs) is needed only for the input checksums of SI
   Table S1, and `FEKAPPA_DFPT` only for the harmonic cross-check.
4. Every queue job, with the sha256 of every file and the build id of its binary, is in `manifest_final.json`
   (identical to `data/provenance/manifest_final.json` in the code archive); `archive_file` names the archive holding it.

Host names are replaced by `HOME` / `NODE` (the two machines that ran the jobs) and account paths by `<home>` /
`<user>`; no other content is changed, so the recorded output sha256s hold for the queue archives. The force constants
are those distributed with the ALAMODE SrTiO₃ example (ALAMODE commit `{sw['alamode_commit'][:7]}`) and are not
redistributed; nor are symbolic links in the pre-queue runs. The two timing probes of the four-phonon runs are
not archived (listed as excluded in the manifest).

## Software

- ALAMODE commit `{sw['alamode_commit']}`. Patches, build recipe and hashes: `tools/alamode_patches/BUILD.md` in the
  code archive. Binary hashes depend on the machine (embedded rpath); each job's build id records the binary that ran it.
- Toolchain lock `env.explicit.txt` (conda-forge, osx-arm64): sha256 `{sw['package_lock_env_explicit_sha256']}`.
- Analysis code: Python package `fekappa` {TAG[1:]}, locked environment `uv.lock`.

| build | tree | patch | machine | anphon sha256 |
|---|---|---|---|---|
{brows}

## Licences

Code: MIT (`LICENSE`). Data, tables and documentation: CC BY 4.0 (`DATA_LICENSE.md`). Literature values are
derived extractions; cite the original sources.
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--toolchain", default=str(ROOT.parent / "_toolchain"))
    ap.add_argument("--skip-raw", action="store_true", help="reuse the variant archives already in the output dir")
    a = ap.parse_args()
    out = ROOT / "release/zenodo_v3"
    out.mkdir(parents=True, exist_ok=True)
    tool = Path(a.toolchain)
    manifest = json.loads((ROOT / "queue/manifest.json").read_text())
    rules, info = zp.build_scrub(manifest, tool)

    # 1. code archive of the tag
    subprocess.run(["git", "archive", "--format=tar.gz", f"--prefix=ferroelectric-soft-mode-transport-{TAG}/",
                    "-o", str(out / CODE), TAG], cwd=ROOT, check=True)
    files = [{"file": CODE, "what": f"code, tables, patches, recorded rules and docs: `git archive {TAG}`"}]

    # 2. queue archives (sha256-checked against the manifest by make_zenodo_package.write_archive)
    jobs_by_variant = {}
    for j in manifest["jobs"].values():
        v = zp.variant_of(j)
        if v:
            jobs_by_variant.setdefault(v, []).append(j)
    jobs_by_variant = {v: jobs_by_variant[v] for v in zp.VARIANT_TEXT if v in jobs_by_variant}
    results = {}
    for v, js in jobs_by_variant.items():
        f = out / f"{v}.tar.gz"
        if a.skip_raw and f.exists():
            results[v] = {"file": f.name, "jobs": len(js), "bytes": f.stat().st_size}
        else:
            results[v] = zp.write_archive(out, v, js, rules, "gz")
            print(f"wrote {f.name}: {results[v]['bytes'] / 1e6:.1f} MB, {results[v]['files']} files")
        files.append({"file": f.name, "what": f"{len(js)} job{'s' if len(js) > 1 else ''}: {zp.VARIANT_TEXT[v]}"})
    if not a.skip_raw:
        problems = zp.verify(out, results, rules, jobs_by_variant)
        if problems:
            raise SystemExit("verify:\n" + "\n".join(problems))

    # 3. pre-queue runs
    pq = prequeue(out, tool / "runs", rules)
    print(f"wrote {pq['file']}: {pq['bytes'] / 1e6:.1f} MB, {pq['files']} files, skipped {pq['skipped']}")
    files.append({"file": pq["file"], "what": "runs made before the job queue ($FEKAPPA_RUNS): " + ", ".join(pq["directories"])})

    # 4. manifest (must equal the tracked export at the tag)
    archive_of = {j["id"]: results[v]["file"] for v, js in jobs_by_variant.items() for j in js}
    sm = json.dumps(zp.scrub_manifest(manifest, rules, archive_of), indent=1, sort_keys=True) + "\n"
    tracked = subprocess.run(["git", "show", f"{TAG}:data/provenance/manifest_final.json"], cwd=ROOT,
                             capture_output=True, text=True, check=True).stdout
    if sm != tracked:
        raise SystemExit("manifest differs from data/provenance/manifest_final.json at the tag")
    (out / "manifest_final.json").write_text(sm)
    files.append({"file": "manifest_final.json", "what": "scrubbed job manifest: 170 jobs, sha256 of every file"})

    # 5. software versions, README, metadata
    sw = {"alamode_commit": zp.ALAMODE_COMMIT, "package_lock_env_explicit_sha256": zp.sha256_file(tool / "env.explicit.txt"),
          "builds": builds(), "code_tag": TAG,
          "code_commit": subprocess.run(["git", "rev-parse", f"{TAG}^{{commit}}"], cwd=ROOT, capture_output=True,
                                        text=True, check=True).stdout.strip()}
    (out / "software_versions.json").write_text(json.dumps(sw, indent=1) + "\n")
    files.append({"file": "software_versions.json", "what": "ALAMODE commit, full binary and patch hashes, lock hash"})
    notes = subprocess.run(["git", "show", f"{TAG}:RELEASE_NOTES.md"], cwd=ROOT, capture_output=True, text=True,
                           check=True).stdout.split("\n---\n", 1)[0]
    meta = {"title": f"Zone-resolved soft-mode damping and the thermal-conductivity decomposition of SrTiO3: code, "
                     f"tables and raw calculation outputs ({TAG})",
            "upload_type": "dataset", "version": TAG[1:],
            "creators": [{"name": "Moid, Syed A.", "affiliation": "Ronin Institute for Independent Scholarship 2.0"}],
            "license": {"code": "MIT", "data_and_docs": "CC-BY-4.0"},
            "related_identifiers": [
                {"relation": "isSupplementTo", "identifier": f"{REPO}/tree/{TAG}", "resource_type": "software"},
                {"relation": "isVersionOf", "identifier": CONCEPT_DOI, "scheme": "doi"}],
            "notes": "Supplement to the manuscript JPCM-126876 (J. Phys.: Condens. Matter, under review).",
            "keywords": ["SrTiO3", "thermal conductivity", "self-consistent phonon", "Wigner transport",
                         "phonon damping", "ALAMODE", "four-phonon scattering"],
            "description_markdown": notes}
    (out / "zenodo_metadata.json").write_text(json.dumps(meta, indent=1, ensure_ascii=False) + "\n")
    files.append({"file": "zenodo_metadata.json", "what": "metadata to enter on the Zenodo form"})
    for f in files:
        f["bytes"] = (out / f["file"]).stat().st_size
    files.insert(len(files), {"file": "README.md", "what": "this file", "bytes": 0})
    (out / "README.md").write_text(readme([f for f in files if f["file"] != "README.md"], sw))
    files[-1]["bytes"] = (out / "README.md").stat().st_size

    # 6. checks, checksums
    for name in ("README.md", "software_versions.json", "zenodo_metadata.json", "manifest_final.json"):
        txt = (out / name).read_text()
        hits = [p.pattern for p, _ in rules if p.search(txt)]
        if hits:
            raise SystemExit(f"{name}: scrub patterns present: {hits}")
    big = [f["file"] for f in files if f["bytes"] > LIMIT]
    if big:
        raise SystemExit(f"over 4 GB: {big}")
    (out / "SHA256SUMS").write_text("".join(f"{zp.sha256_file(out / f['file'])}  {f['file']}\n" for f in files))
    total = sum(f["bytes"] for f in files) + (out / "SHA256SUMS").stat().st_size
    print(f"{len(files) + 1} files, {total / 1e9:.2f} GB in {out.relative_to(ROOT)}")
    for f in files:
        print(f"  {f['bytes'] / 1e6:10.1f} MB  {f['file']}")


if __name__ == "__main__":
    main()
