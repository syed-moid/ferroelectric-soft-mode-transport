"""Build the Zenodo data package for the revision calculations (raw ALAMODE outputs + scrubbed job manifest).

Reads the job manifest (`queue/manifest.json`) and the raw outputs under `data/raw/alamode_sto/revision/`; the
originals are never modified. For each archive variant, every job directory is read, each file is checked against
the sha256 recorded in the manifest, host and account strings are replaced by neutral tokens in memory, and the
result is written to one deterministic tarball (sorted members, fixed mtime, uid/gid 0).

Output (`release/zenodo_v3/`, gitignored):
    <variant>.tar.zst (or .tar.gz when no zstd binary is on PATH)   one per variant
    manifest_scrubbed.json      job records: hosts -> HOME / NODE, build ids keep the anphon sha256 and lose the host
    software_versions.json      ALAMODE commit, patch and binary sha256, package-lock sha256
    README.md                   layout, variants, job counts, software versions
    SHA256SUMS                  for every file above

Scrub list: the hostnames found in the manifest, the account names of the two machines (from `--toolchain` path,
the import paths in the manifest and the current home directory), any home-directory prefix, plus optional
extra strings, one per line, in `release/scrub_extra.txt` (kept out of git so that no identifying string is
tracked).

Usage:
    python tools/zenodo/make_zenodo_package.py --dry-run       # listing and size estimate only
    python tools/zenodo/make_zenodo_package.py                 # build and verify
Options: --toolchain PATH (default ../_toolchain relative to the repository root), --out DIR.
"""
from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import io
import json
import re
import shutil
import subprocess
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW_REL = "data/raw/alamode_sto/revision"
ALAMODE_COMMIT = "c5b015caaf72c1a86a5250b6e34a86c0e09c7986"
HOME_PARENT = str(Path.home().parent)  # directory that holds the account home directories
HOME_RE = re.compile(re.escape(HOME_PARENT) + r"/([^/]+)/")
META_KEEP = ("created", "host", "production_q1", "production_temps", "production_missing_T_ge_tmin")
MTIME = 1760000000  # fixed member mtime (2025-10-09) so that archives are byte-reproducible

# Jobs that are not part of the record: timing probes for run design, not reported anywhere.
PROBES = {
    "m4ph_probe_q8_m4": "timing probe for the four-phonon pilot design (q1 = 8^3); not a reported result",
    "m4ph_pm_q12_T300_m6_e5_timing": "wall-time probe at 6^3 for the four-phonon runs; not a reported result",
}


def variant_of(job: dict) -> str | None:
    """Archive variant for one manifest job (None: no outputs, or excluded)."""
    jid, kind, v = job["id"], job["kind"], job["params"].get("variant")
    if jid in PROBES or not job.get("archive"):
        return None
    if jid == "scph_pm_q12":
        return "scph_q12"
    if jid == "scph_pm_q16" or v == "pm_q16":
        return "scph_q16"
    if jid == "scph_pm_q8" or v in ("pm_q8", "pm_q8_mix0p1"):
        return "pm_q8"
    if kind == "rta" and v == "pm_q12":
        return "pm_q12_rta"
    if kind == "evec" and v == "pm_q12":
        return "pm_q12_evec"
    if kind == "mode4ph":
        return "pm_q12_4ph_audit"
    raise ValueError(f"no variant rule for job {jid}")


VARIANT_TEXT = {
    "scph_q12": "Production SCPH solution, q1 = 12^3 with polarization mixing: input, log, renormalized "
                "dispersion and the temperature-dependent FC2 used by every transport run.",
    "scph_q16": "SCPH on q1 = 16^3 (q1-convergence check) and the four RTA runs on it (12^3 transport mesh; "
                "250, 300, 500, 800 K), including the 250 K boundary point.",
    "pm_q8": "SCPH on q1 = 8^3 (q1-sampling comparison; three RTA runs on it) and the q1 = 8^3 SCPH run "
                      "with a reduced mixing parameter (R-point offset check).",
    "pm_q12_rta": "Production Wigner-RTA runs on q1 = 12^3: kappa_P, kappa_C and the per-pair coherence elements "
                  "for every temperature and transport mesh.",
    "pm_q12_evec": "Eigenvectors on each transport mesh and temperature (sector labelling by eigenvector "
                   "character).",
    "pm_q12_4ph_audit": "Per-mode four-phonon runs (diagram c): the unmodified 4^3 pilot and 6^3 check, the rule build "
                        "(enlarged identity exclusion, 1e-9 Ry), the 1e-6 Ry tolerance build, the one-quartet trace and the "
                        "three validation builds (omitted sum, strict-index exclusion, full-precision dump), including the "
                        "descriptive and heat-carrier samples.",
}


# ----------------------------------------------------------------------------- scrubbing
def build_scrub(manifest: dict, toolchain: Path) -> tuple[list[tuple[re.Pattern, str]], dict]:
    home_host = manifest["meta"]["host"]
    hosts = {j["host"] for j in manifest["jobs"].values() if j.get("host")}
    hosts |= {h["host"] for j in manifest["jobs"].values() for h in j.get("history", []) if h.get("host")}
    hosts.add(home_host)
    host_token = {h: ("HOME" if h == home_host else "NODE") for h in hosts}
    for h in list(host_token):                      # short form without the mDNS suffix
        if h.endswith(".local"):
            host_token[h[:-6]] = host_token[h]
    accounts = {Path.home().name}
    for j in manifest["jobs"].values():
        d = (j.get("import") or {}).get("dir", "")
        m = HOME_RE.match(d)
        if m:
            accounts.add(m.group(1))
    m = HOME_RE.match(str(toolchain.resolve()))
    if m:
        accounts.add(m.group(1))
    extra_file = ROOT / "release/scrub_extra.txt"
    extra = [s.strip() for s in extra_file.read_text().splitlines()] if extra_file.exists() else []
    rules: list[tuple[re.Pattern, str]] = []
    # directory prefixes first, so that no part of a local path survives: run directories, the toolchain, this
    # repository, then any remaining home directory
    for path, token in ((toolchain.resolve() / "runs", "<runs>"), (toolchain.resolve(), "<toolchain>"), (ROOT, "<repo>")):
        rules.append((re.compile(re.escape(str(path))), token))
        try:  # the same directory in its home-relative (tilde) form
            rules.append((re.compile(re.escape(str(Path("~") / path.relative_to(Path.home())))), token))
        except ValueError:
            pass
    rules.append((re.compile(re.escape(HOME_PARENT) + r"/[^\s\"'<>]+"), "<home>"))
    for h in sorted(host_token, key=len, reverse=True):
        rules.append((re.compile(re.escape(h)), host_token[h]))
    for a in sorted(accounts | {s for s in extra if s}, key=len, reverse=True):
        rules.append((re.compile(re.escape(a)), "<user>"))
    info = {"hosts": host_token, "n_accounts": len(accounts), "n_extra": len([s for s in extra if s])}
    return rules, info


def scrub_bytes(data: bytes, rules) -> tuple[bytes, int]:
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return data, 0                               # binary: left as is, reported by verify() if it matches
    n = 0
    for pat, tok in rules:
        text, k = pat.subn(tok, text)
        n += k
    return text.encode("utf-8"), n


def scrub_obj(obj, rules):
    if isinstance(obj, str):
        return scrub_bytes(obj.encode(), rules)[0].decode()
    if isinstance(obj, list):
        return [scrub_obj(x, rules) for x in obj]
    if isinstance(obj, dict):
        return {k: scrub_obj(v, rules) for k, v in obj.items()}
    return obj


def scrub_manifest(manifest: dict, rules, archive_of: dict) -> dict:
    m = copy.deepcopy(manifest)
    m["meta"] = {k: v for k, v in m["meta"].items() if k in META_KEEP}  # drop internal notes, keep provenance
    for jid, j in m["jobs"].items():
        if isinstance(j.get("build"), str):          # "alamode@c5b015ca/anphon-sha256:<h>/<host>" -> drop the host
            parts = j["build"].split("/")
            if len(parts) == 3 and "sha256" in parts[1]:
                j["build"] = "/".join(parts[:2])
        j["archive_file"] = archive_of.get(jid)
        if jid in PROBES:
            j["excluded"] = PROBES[jid]
    m = scrub_obj(m, rules)
    m["meta"]["scrubbed"] = "host names -> HOME / NODE; local paths -> <runs> / <toolchain> / <repo> / <home>; account names -> <user>"
    return m


# ----------------------------------------------------------------------------- software versions
def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def software_versions(manifest: dict, toolchain: Path) -> dict:
    out = {"alamode_commit": ALAMODE_COMMIT,
           "anphon_sha256_in_build_ids": sorted({re.search(r"anphon-sha256:([0-9a-f]+)", j["build"]).group(1)
                                                 for j in manifest["jobs"].values()
                                                 if isinstance(j.get("build"), str) and "anphon-sha256" in j["build"]}),
           "patches": {}, "binaries_this_machine": {}}
    if not toolchain.is_dir():
        out["note"] = "toolchain directory not found; patch and lock hashes not recorded"
        return out
    lock = toolchain / "env.explicit.txt"
    if lock.exists():
        out["package_lock_env_explicit_sha256"] = sha256_file(lock)
    for tree, desc in (("alamode_4ph_audit", "identity-set exclusion in selfenergy_c (1e-9 Ry), vs c5b015c"),
                       ("alamode_absw", "|omega| bare-harmonic diagnostic build, vs c5b015c")):
        t = toolchain / tree
        if (t / ".git").exists():
            diff = subprocess.run(["git", "-C", str(t), "diff", ALAMODE_COMMIT], capture_output=True, check=True).stdout
            out["patches"][tree] = {"description": desc, "base": ALAMODE_COMMIT,
                                    "git_diff_sha256": hashlib.sha256(diff).hexdigest(),
                                    "lines_changed": sum(1 for ln in diff.decode().splitlines()
                                                         if ln[:1] in "+-" and ln[:3] not in ("+++", "---"))}
    a, b = toolchain / "alamode_4ph_audit/anphon/selfenergy.cpp", toolchain / "alamode_4ph_audit_tol1e6/anphon/selfenergy.cpp"
    if a.exists() and b.exists():
        diff = subprocess.run(["diff", "-u", "--label", "a/anphon/selfenergy.cpp", "--label", "b/anphon/selfenergy.cpp",
                               str(a), str(b)], capture_output=True).stdout
        out["patches"]["alamode_4ph_audit_tol1e6"] = {
            "description": "identity tolerance 1e-6 Ry instead of 1e-9 Ry, vs alamode_4ph_audit",
            "base": "alamode_4ph_audit", "diff_u_sha256": hashlib.sha256(diff).hexdigest()}
    for tree in ("alamode", "alamode_4ph_audit", "alamode_absw"):
        p = toolchain / tree / "_build/anphon/anphon"
        if p.exists():
            out["binaries_this_machine"][tree] = sha256_file(p)
    out["binaries_note"] = ("binary hashes differ between the two machines because absolute rpaths are embedded; "
                            "the build ids in the manifest carry the hash of the binary that ran each job")
    return out


# ----------------------------------------------------------------------------- archives
def have_zstd() -> bool:
    return shutil.which("zstd") is not None


def write_archive(out: Path, variant: str, jobs: list[dict], rules, compress: str) -> dict:
    raw = ROOT / RAW_REL
    tar_path = out / f"{variant}.tar"
    n_files = n_sub = raw_bytes = 0
    with open(tar_path, "wb") as fh, tarfile.open(fileobj=fh, mode="w", format=tarfile.PAX_FORMAT) as tf:
        for job in sorted(jobs, key=lambda j: j["id"]):
            d = ROOT / job["archive"]
            files = sorted(p for p in d.iterdir() if p.is_file())
            recorded = job.get("outputs", {})
            for p in files:
                data = p.read_bytes()
                if p.name in recorded and hashlib.sha256(data).hexdigest() != recorded[p.name]:
                    raise SystemExit(f"sha256 mismatch against the manifest: {p.relative_to(ROOT)}")
                data, k = scrub_bytes(data, rules)
                n_sub += k
                ti = tarfile.TarInfo(f"{variant}/{p.relative_to(raw)}")
                ti.size, ti.mtime, ti.mode = len(data), MTIME, 0o644
                ti.uid = ti.gid = 0
                ti.uname = ti.gname = ""
                tf.addfile(ti, io.BytesIO(data))
                n_files += 1
                raw_bytes += len(data)
            missing = set(recorded) - {p.name for p in files}
            if missing:
                raise SystemExit(f"{job['id']}: outputs in the manifest but not on disk: {sorted(missing)}")
    if compress == "zst":
        final = out / f"{variant}.tar.zst"
        subprocess.run(["zstd", "-19", "-T0", "-q", "-f", "--rm", str(tar_path), "-o", str(final)], check=True)
    else:
        final = out / f"{variant}.tar.gz"
        with open(tar_path, "rb") as src, open(final, "wb") as dst, \
                gzip.GzipFile(filename="", mode="wb", fileobj=dst, mtime=0, compresslevel=9) as gz:
            shutil.copyfileobj(src, gz, 1 << 20)
        tar_path.unlink()
    return {"file": final.name, "jobs": len(jobs), "files": n_files, "uncompressed_bytes": raw_bytes,
            "bytes": final.stat().st_size, "substitutions": n_sub}


def open_tar(p: Path):
    if p.suffix == ".zst":
        data = subprocess.run(["zstd", "-dc", str(p)], capture_output=True, check=True).stdout
        return tarfile.open(fileobj=io.BytesIO(data), mode="r:")
    return tarfile.open(p, mode="r:gz")


def verify(out: Path, results: dict, rules, jobs_by_variant: dict) -> list[str]:
    """Re-read every archive: member count, expected job directories, no scrub pattern left in any member."""
    problems = []
    for v, r in results.items():
        with open_tar(out / r["file"]) as tf:
            names = []
            for ti in tf:
                names.append(ti.name)
                data = tf.extractfile(ti).read()
                for pat, _ in rules:
                    if pat.search(data.decode("utf-8", "replace")):
                        problems.append(f"{r['file']}:{ti.name}: matches {pat.pattern}")
        if len(names) != r["files"]:
            problems.append(f"{r['file']}: {len(names)} members, expected {r['files']}")
        dirs = {n.split("/")[2] for n in names}
        want = {Path(j["archive"]).name for j in jobs_by_variant[v]}
        if dirs != want:
            problems.append(f"{r['file']}: job directories differ: {sorted(dirs ^ want)}")
    return problems


# ----------------------------------------------------------------------------- README
def readme(results: dict, jobs_by_variant: dict, sw: dict, compress: str) -> str:
    ext = "tar.zst" if compress == "zst" else "tar.gz"
    rows = "\n".join(f"| `{r['file']}` | {r['jobs']} | {r['files']} | {r['uncompressed_bytes'] / 1e6:,.1f} | "
                     f"{r['bytes'] / 1e6:,.1f} | {VARIANT_TEXT[v]} |" for v, r in results.items())
    pat = "\n".join(f"- `{k}`: {p['description']}; base `{p['base']}`; "
                    f"sha256 of the diff `{p.get('git_diff_sha256') or p.get('diff_u_sha256')}`"
                    for k, p in sw.get("patches", {}).items()) or "- not recorded on this machine"
    probes = "\n".join(f"- `{k}`: {v}" for k, v in PROBES.items())
    unpack = "zstd -dc <file>.tar.zst | tar -x" if compress == "zst" else "tar -xzf <file>.tar.gz"
    return f"""# SrTiO3 SCPH + Wigner transport: raw calculation outputs (revision data set)

Raw ALAMODE outputs behind the tables in `data/processed/` of the code repository. Every file here is listed,
with its sha256, in `manifest_scrubbed.json` (one record per job: parameters, build id, timings, outputs and
the metrics extracted at validation time). Host names are replaced by `HOME` / `NODE` (the two machines that ran
the jobs) and account paths by `<home>` / `<user>`; no other content is changed.

## Layout

Each archive unpacks to `<variant>/<source directory>/<job id>/`, where the source directory is the run variant
used in the code repository (`pm_q12`, `pm_q16`, `pm_q8`, `pm_q8_mix0p1`, `pm_q12_4ph_audit`,
`pm_q12_4ph_audit_tol1e6`). Unpack with `{unpack}`.

| archive | jobs | files | uncompressed (MB) | compressed (MB) | content |
|---|---|---|---|---|---|
{rows}

Typical files per job: `*.in` / `*.log` (exact anphon input and log), `dfc2.log` (temperature-dependent FC2
from the SCPH solution), `*.result`, `*.kl` (particle conductivity), `*.kl_coherent` and `*.kc_elem` (coherence
conductivity and its per-pair elements), `*.mesh.evec` (eigenvectors), `*.scph_bands` / `*.scph_dfc2` (SCPH
output), `*.Gamma.*` (per-mode four-phonon self-energy runs).

Force constants are those of the ALAMODE SrTiO3 example (`example/` of ALAMODE {ALAMODE_COMMIT[:7]}); they are not
redistributed here.

## Excluded jobs

{probes}

The job `d0_decision` is a bookkeeping record (no outputs); it is kept in the manifest.

## Software

- ALAMODE commit `{sw['alamode_commit']}`; anphon binaries by sha256 prefix in the build ids:
  {', '.join(f'`{h}`' for h in sw['anphon_sha256_in_build_ids'])}. {sw.get('binaries_note', '')}.
- Package lock (`env.explicit.txt`, conda-forge, 111 pinned packages): sha256
  `{sw.get('package_lock_env_explicit_sha256', 'not recorded')}`.
- Source patches (the patch files and build scripts are in the code repository, `tools/alamode_patches/`):
{pat}

Full hashes: `software_versions.json`. Archive checksums: `SHA256SUMS` (`shasum -a 256 -c SHA256SUMS`).
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--toolchain", default=str(ROOT.parent / "_toolchain"))
    ap.add_argument("--out", default=str(ROOT / "release/zenodo_v3"))
    a = ap.parse_args()
    manifest = json.loads((ROOT / "queue/manifest.json").read_text())
    toolchain, out = Path(a.toolchain), Path(a.out)
    jobs_by_variant: dict[str, list[dict]] = {}
    for j in manifest["jobs"].values():
        v = variant_of(j)
        if v:
            jobs_by_variant.setdefault(v, []).append(j)
    order = [v for v in VARIANT_TEXT if v in jobs_by_variant]
    jobs_by_variant = {v: jobs_by_variant[v] for v in order}
    rules, info = build_scrub(manifest, toolchain)
    compress = "gz"  # tar.gz for every archive: readable without extra tools
    for v, js in jobs_by_variant.items():
        size = sum(p.stat().st_size for j in js for p in (ROOT / j["archive"]).iterdir() if p.is_file())
        print(f"{v:16s} {len(js):4d} jobs  {size / 1e6:9.1f} MB uncompressed")
    print(f"excluded: {sorted(PROBES)}; scrub: {len(info['hosts'])} host strings, {info['n_accounts']} accounts, "
          f"{info['n_extra']} extra; compression: {compress}")
    if a.dry_run:
        return
    out.mkdir(parents=True, exist_ok=True)
    results = {}
    for v, js in jobs_by_variant.items():
        results[v] = write_archive(out, v, js, rules, compress)
        print(f"wrote {results[v]['file']}: {results[v]['bytes'] / 1e6:.1f} MB, {results[v]['substitutions']} substitutions")
    archive_of = {j["id"]: results[v]["file"] for v, js in jobs_by_variant.items() for j in js}
    sm = scrub_manifest(manifest, rules, archive_of)
    (out / "manifest_scrubbed.json").write_text(json.dumps(sm, indent=1, sort_keys=True) + "\n")
    sw = software_versions(manifest, toolchain)
    (out / "software_versions.json").write_text(json.dumps(scrub_obj(sw, rules), indent=1) + "\n")
    (out / "README.md").write_text(readme(results, jobs_by_variant, sw, compress))
    for name in ("manifest_scrubbed.json", "software_versions.json", "README.md"):
        txt = (out / name).read_text()
        hits = [p.pattern for p, _ in rules if p.search(txt)]
        if hits:
            raise SystemExit(f"{name}: scrub patterns still present: {hits}")
    files = sorted(r["file"] for r in results.values()) + ["README.md", "manifest_scrubbed.json", "software_versions.json"]
    (out / "SHA256SUMS").write_text("".join(f"{sha256_file(out / f)}  {f}\n" for f in files))
    problems = verify(out, results, rules, jobs_by_variant)
    print("verify:", "OK" if not problems else "\n".join(problems))
    total = sum((out / f).stat().st_size for f in files)
    print(f"total {total / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
