#!/usr/bin/env python3
"""Checkpointed, resumable local queue for the ALAMODE (anphon) revision runs.

Subcommands
    run      run jobs until all are done, STOP / STOP_NOW is set, or the --until time is reached (default)
    init     create the manifest with the D0 jobs (and register already-computed units)
    status   regenerate queue/status.md and print it
    requeue  put a failed job back to pending (attempt counter reset)

Design
* A job is one atomic anphon run (SCPH sweep or one RTA point). A job is DONE only after its
  outputs exist, parse and pass validation; anything else is re-run from scratch, in a freshly
  emptied work directory (partial outputs are purged before every start).
* queue/manifest.json is the single source of truth; every read-modify-write holds an exclusive
  flock on queue/manifest.lock and writes via temp-file + rename.
* A running job holds queue/locks/<job>.lock with host, pid and a heartbeat refreshed every
  HEARTBEAT_S. On start the runner re-queues running jobs whose runner pid is dead (same host)
  or whose heartbeat is older than STALE_S (any host).
* Automatic retries: a job that fails (non-zero exit or validation) is retried at most
  MAX_RETRIES times, then marked failed. A job killed by STOP_NOW or at the --until time goes back to
  pending without consuming an attempt.
* queue/STOP stops gracefully after the current job; queue/STOP_NOW kills the running job and
  re-queues it. Both files are consumed (deleted, logged) by the next runner start.
* --until "YYYY-MM-DDTHH:MM:SS <IANA zone>": no job is started whose expected duration would cross that
  time; at that time a running job is killed and re-queued.
* Optional per-job params (absent = unchanged behaviour):
    mpi_np        run anphon under mpirun -np N with OMP_NUM_THREADS = 1;
    anphon_tree   use TOOL/<tree>/_build/anphon/anphon (an audit build) instead of the production
                  binary; the job's build id is then prefixed AUDIT:<tree> so it differs visibly;
    mixalpha      SCPH jobs only: MIXALPHA written into scph.in instead of the template's 0.2 (checked
                  against the value echoed in scph.log);
    measure_mem   wrap anphon in /usr/bin/time -l (time_l.txt, archived; peak RSS in metrics);
    rss_kill_gb   watchdog: kill the job (failed, no retry) if the job's process group exceeds
                  this resident set size; the sampled peak is recorded in metrics.

Standard library only (runs under /usr/bin/python3, 3.9+).
"""

import argparse
import fcntl
import hashlib
import json
import os
import shutil
import signal
import socket
import statistics
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

try:
    from zoneinfo import ZoneInfo
except ImportError:  # pragma: no cover
    ZoneInfo = None

REPO = Path(__file__).resolve().parent.parent
QDIR = Path(os.environ.get("QUEUE_DIR", str(REPO / "queue")))
MANIFEST = QDIR / "manifest.json"
MLOCK = QDIR / "manifest.lock"
LOCKS = QDIR / "locks"
STATUS = QDIR / "status.md"
PIDFILE = QDIR / "runner.pid"
STOP, STOP_NOW = QDIR / "STOP", QDIR / "STOP_NOW"

TOOL = Path(os.environ.get("QUEUE_TOOLCHAIN", "QUEUE_TOOLCHAIN-unset"))  # required; checked in main()
WORK = TOOL / "runs" / "queue"
EXAMPLE = TOOL / "sto_example"
ANPHON = TOOL / "alamode/_build/anphon/anphon"
DFC2 = TOOL / "alamode/_build/tools/dfc2"
ENV_BIN = TOOL / "env/bin"
ARCHIVE = REPO / "data/raw/alamode_sto/revision"
RTA_TEMPLATE = REPO / "data/raw/alamode_sto/RTA_scph_300K.in"
SCPH_TEMPLATE = EXAMPLE / "scph.in"
DECISION_MD = QDIR / "q1_decision.md"  # record of the q1 rule (docs/recorded_rules.md)

HOST = socket.gethostname()
HEARTBEAT_S = 60
STALE_S = 600
POLL_S = 10
MAX_RETRIES = 2
OMP_THREADS = int(os.environ.get("QUEUE_OMP_THREADS", "12"))
CLOSURE_TOL = 5e-3  # kc_elem pair-sum vs solver kappa_C (relative)
PRODUCTION_TMIN = 250

# default expected durations (s), replaced by measured medians as jobs complete
DEFAULT_EST = {("scph", 8): 2400, ("scph", 12): 9000, ("evec", 12): 60, ("evec", 16): 150, ("evec", 20): 300,
               ("rta", 8): 240, ("rta", 12): 1600, ("rta", 16): 6600, ("rta", 20): 25000}


# --------------------------------------------------------------------------- utilities
def now_iso():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def log(msg):
    print(f"[{now_iso()}] {msg}", flush=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def pid_alive(pid):
    try:
        os.kill(int(pid), 0)
        return True
    except (OSError, ValueError):
        return False


class Manifest:
    """Context manager: exclusive flock, load, and atomic save on exit (if modified)."""

    def __enter__(self):
        QDIR.mkdir(exist_ok=True)
        self._lf = open(MLOCK, "a+")
        fcntl.flock(self._lf, fcntl.LOCK_EX)
        self.data = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {"jobs": {}, "meta": {}}
        return self.data

    def __exit__(self, exc_type, exc, tb):
        if exc_type is None:
            tmp = MANIFEST.with_suffix(".json.tmp")
            tmp.write_text(json.dumps(self.data, indent=1, sort_keys=True))
            os.replace(tmp, MANIFEST)
        fcntl.flock(self._lf, fcntl.LOCK_UN)
        self._lf.close()
        return False


def build_id():
    rev = "unknown"
    try:
        rev = subprocess.run(["git", "-C", str(TOOL / "alamode"), "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True).stdout.strip() or rev
    except Exception:
        pass
    return f"alamode@{rev}/anphon-sha256:{sha256(ANPHON)[:12]}/{HOST}"


def parse_until(s):
    if not s:
        return None
    stamp, _, tzname = s.partition(" ")
    dt = datetime.fromisoformat(stamp)
    if tzname:
        dt = dt.replace(tzinfo=ZoneInfo(tzname))
    return dt.timestamp()


# --------------------------------------------------------------------------- job specs
def scph_job(q1, pm=True, imported_from=None, priority=10):
    jid = f"scph_{'pm' if pm else 'nopm'}_q{q1}"
    job = {"id": jid, "kind": "scph", "priority": priority, "deps": [],
           "params": {"q1": q1, "pm": pm, "tmin": 0, "tmax": 1000, "dt": 50}}
    if imported_from:
        job["import"] = imported_from
    return job


def rta_job(q1, T, mesh, pm=True, imported_from=None, priority=20, extra=None):
    variant = f"{'pm' if pm else 'nopm'}_q{q1}"
    jid = f"rta_{variant}_T{T}_m{mesh}"
    job = {"id": jid, "kind": "rta", "priority": priority, "deps": [f"scph_{variant}"],
           "params": {"q1": q1, "pm": pm, "T": T, "mesh": mesh, "variant": variant,
                      "extra": extra or {}}}
    if imported_from:
        job["import"] = imported_from
    return job


def new_entry(job):
    e = dict(job)
    e.update({"status": "pending", "attempts": 0, "history": [], "outputs": {}, "metrics": {}})
    return e


def prefix_of(job):
    return job.get("import", {}).get("prefix", job["id"])


def workdir_of(job):
    imp = job.get("import")
    return Path(imp["dir"]).expanduser() if imp else WORK / job["id"]


# --------------------------------------------------------------------------- inputs
def write_scph_inputs(job, wd):
    p = job["params"]
    s = SCPH_TEMPLATE.read_text()
    s = s.replace("PREFIX = STO_scph2-2", f"PREFIX = {job['id']}")
    s = s.replace("SELF_OFFDIAG = 0", f"SELF_OFFDIAG = {1 if p['pm'] else 0}")
    s = s.replace("KMESH_SCPH = 2 2 2", f"KMESH_SCPH = {p['q1']} {p['q1']} {p['q1']}")
    s = s.replace("TMIN =0; TMAX = 1000; DT = 50", f"TMIN = {p['tmin']}; TMAX = {p['tmax']}; DT = {p['dt']}")
    if p.get("mixalpha") is not None:
        s = s.replace("MIXALPHA = 0.2", f"MIXALPHA = {p['mixalpha']}")
        assert f"MIXALPHA = {p['mixalpha']}" in s, "MIXALPHA"
    for must in (f"PREFIX = {job['id']}", f"KMESH_SCPH = {p['q1']}", f"SELF_OFFDIAG = {1 if p['pm'] else 0}"):
        assert must in s, must
    (wd / "scph.in").write_text(s)
    (wd / "STO_anharm.xml").symlink_to(EXAMPLE / "STO_anharm.xml")
    shutil.copy(EXAMPLE / "BORN", wd / "BORN")
    return f'"{ANPHON}" scph.in > scph.log 2>&1'


def dfc2_file_for(scph_entry):
    return workdir_of(scph_entry) / f"{prefix_of(scph_entry)}.scph_dfc2"


def write_rta_inputs(job, wd, manifest):
    p = job["params"]
    scph = manifest["jobs"][job["deps"][0]]
    s = RTA_TEMPLATE.read_text()
    s = s.replace("PREFIX = STO_RTA_scph_300K", f"PREFIX = {job['id']}")
    s = s.replace("FC2XML = STO_renorm_300K.xml", "FC2XML = renorm.xml")
    s = s.replace("TMIN = 300; TMAX = 300", f"TMIN = {p['T']}; TMAX = {p['T']}")
    s = s.replace("\n8 8 8\n", f"\n{p['mesh']} {p['mesh']} {p['mesh']}\n")
    for k, v in (p.get("extra") or {}).items():
        s = s.replace("\n/\n", f"\n {k} = {v}\n/\n", 1)
    for must in (f"PREFIX = {job['id']}", f"TMIN = {p['T']}", f"\n{p['mesh']} {p['mesh']} {p['mesh']}\n",
                 "KAPPA_COHERENT = 2", "NONANALYTIC = 3"):
        assert must in s, must
    (wd / "rta.in").write_text(s)
    (wd / "STO_anharm.xml").symlink_to(EXAMPLE / "STO_anharm.xml")
    shutil.copy(EXAMPLE / "BORN", wd / "BORN")
    shutil.copy(EXAMPLE / "STO222.xml", wd / "STO222.xml")
    return (f'"{DFC2}" STO222.xml renorm.xml "{dfc2_file_for(scph)}" {p["T"]} > dfc2.log 2>&1 && '
            f'"{ANPHON}" rta.in > rta.log 2>&1')


def write_mode4ph_inputs(job, wd, manifest):
    """Four-phonon per-mode pilot: MODE = RTA + KS_INPUT + QUARTIC = 2 (diagram c), FC2 at T_fc2."""
    p = job["params"]
    scph = manifest["jobs"][job["deps"][0]]
    s = RTA_TEMPLATE.read_text()
    s = s.replace("PREFIX = STO_RTA_scph_300K", f"PREFIX = {job['id']}")
    s = s.replace("FC2XML = STO_renorm_300K.xml", "FC2XML = renorm.xml")
    s = s.replace("TMIN = 300; TMAX = 300; DT = 50", f"TMIN = {p['tmin']}; TMAX = {p['tmax']}; DT = {p['dt']}")
    s = s.replace("\n8 8 8\n", f"\n{p['mesh']} {p['mesh']} {p['mesh']}\n")
    s = s.replace(" KAPPA_COHERENT = 2\n", " KS_INPUT = modes.txt\n QUARTIC = 2\n SELF_ENERGY = 1\n")
    s = s.replace("\n/\n", f"\n EPSILON = {p['epsilon']}\n/\n", 1)
    for must in ("QUARTIC = 2", "SELF_ENERGY = 1", "KS_INPUT = modes.txt", f"EPSILON = {p['epsilon']}", f"TMIN = {p['tmin']}"):
        assert must in s, must
    (wd / "rta.in").write_text(s)
    (wd / "modes.txt").write_text(f"{len(p['modes'])}\n" + "".join(f"{a} {b} {c} {n}\n" for a, b, c, n in p["modes"]))
    (wd / "STO_anharm.xml").symlink_to(EXAMPLE / "STO_anharm.xml")
    shutil.copy(EXAMPLE / "BORN", wd / "BORN")
    shutil.copy(EXAMPLE / "STO222.xml", wd / "STO222.xml")
    return (f'"{DFC2}" STO222.xml renorm.xml "{dfc2_file_for(scph)}" {p["T_fc2"]} > dfc2.log 2>&1 && '
            f'"{ANPHON}" rta.in > rta.log 2>&1')


def evec_job(q1, T, mesh, pm=True, priority=29):
    variant = f"{'pm' if pm else 'nopm'}_q{q1}"
    return {"id": f"evec_{variant}_T{T}_m{mesh}", "kind": "evec", "priority": priority,
            "deps": [f"scph_{variant}"], "params": {"q1": q1, "pm": pm, "T": T, "mesh": mesh, "variant": variant}}


def write_evec_inputs(job, wd, manifest):
    """MODE = phonons on the transport mesh with the same FC2XML / NAC as the RTA run; PRINTEVEC = 1."""
    p = job["params"]
    scph = manifest["jobs"][job["deps"][0]]
    s = (f"&general\n PREFIX = {job['id']}\n MODE = phonons\n NKD = 3; KD = Sr Ti O\n FCSXML = renorm.xml\n"
         f" NONANALYTIC = 3; BORNINFO = BORN\n/\n&cell\n7.363\n1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n/\n"
         f"&kpoint\n2\n{p['mesh']} {p['mesh']} {p['mesh']}\n/\n&analysis\n PRINTEVEC = 1\n/\n")
    (wd / "ph.in").write_text(s)
    shutil.copy(EXAMPLE / "BORN", wd / "BORN")
    shutil.copy(EXAMPLE / "STO222.xml", wd / "STO222.xml")
    return (f'"{DFC2}" STO222.xml renorm.xml "{dfc2_file_for(scph)}" {p["T"]} > dfc2.log 2>&1 && '
            f'"{ANPHON}" ph.in > ph.log 2>&1')


def validate_evec(job, wd):
    f = wd / f"{job['id']}.mesh.evec"
    if not f.exists():
        raise ValidationError(f"missing {f.name}")
    nk = sum(1 for l in open(f) if l.startswith("## kpoint"))
    if nk != job["params"]["mesh"] ** 3:
        raise ValidationError(f"{f.name}: {nk} k-points, expected {job['params']['mesh'] ** 3}")
    return {"nk": nk}, [f, wd / "ph.in", wd / "ph.log", wd / "dfc2.log"]


# --------------------------------------------------------------------------- validation
class ValidationError(Exception):
    pass


def _diag(path):
    rows = [l.split() for l in Path(path).read_text().splitlines() if l.strip() and not l.startswith("#")]
    if len(rows) != 1:
        raise ValidationError(f"{path.name}: expected one temperature row, got {len(rows)}")
    v = [float(x) for x in rows[0]]
    return v


def validate_rta(job, wd):
    pre = wd / prefix_of(job)
    files = {ext: Path(f"{pre}.{ext}") for ext in ("kl", "kl_coherent", "result", "kc_elem")}
    for ext, f in files.items():
        if not f.exists() or f.stat().st_size == 0:
            raise ValidationError(f"missing/empty {f.name}")
    kl = _diag(files["kl"])
    if abs(kl[0] - job["params"]["T"]) > 1e-6:
        raise ValidationError(f".kl temperature {kl[0]} != {job['params']['T']}")
    kP = [kl[1], kl[5], kl[9]]
    if max(kP) - min(kP) > 1e-3 * max(kP) or not all(x > 0 for x in kP):
        raise ValidationError(f".kl not isotropic/positive: {kP}")
    klc = _diag(files["kl_coherent"])
    kC = sum(klc[1:4]) / 3
    lines = files["result"].read_text().splitlines()
    i = lines.index("#KPOINT")
    grid = [int(x) for x in lines[i + 1].split()]
    nk, n_irred = grid[0] * grid[1] * grid[2], int(lines[i + 2])
    if grid != [job["params"]["mesh"]] * 3:
        raise ValidationError(f".result mesh {grid}")
    n_gamma = sum(1 for l in lines if l.startswith("#GAMMA_EACH"))
    if n_gamma != 15 * n_irred:
        raise ValidationError(f".result has {n_gamma} GAMMA_EACH blocks, expected {15 * n_irred}")
    if any(("nan" in l.lower()) for l in lines if not l.startswith("#")):
        raise ValidationError(".result contains NaN")
    pair = 0.0
    with open(files["kc_elem"]) as fh:
        for l in fh:
            if l.startswith("#"):
                continue
            c = l.split()
            if len(c) < 10:
                continue
            if c[1] == c[2] and c[3] != c[4]:
                pair += float(c[8])
    pair = pair / 3.0 / nk
    closure = abs(pair - kC) / kC if kC > 0 else abs(pair)
    if closure > CLOSURE_TOL:
        raise ValidationError(f"kc_elem pair-sum {pair:.5f} vs kappa_C {kC:.5f} (rel {closure:.2e})")
    metrics = {"kappa_P": sum(kP) / 3, "kappa_C": kC, "kc_pairsum": pair, "closure_rel": closure,
               "n_irred": n_irred, "nk": nk}
    outs = [files[e] for e in files] + [wd / x for x in ("rta.in", "rta.log", "dfc2.log") if (wd / x).exists()]
    return metrics, outs


def validate_scph(job, wd):
    pre = wd / prefix_of(job)
    log_path = wd / "scph.log"
    dfc2 = Path(f"{pre}.scph_dfc2")
    bands = Path(f"{pre}.scph_bands")
    for f in (log_path, dfc2, bands):
        if not f.exists() or f.stat().st_size == 0:
            raise ValidationError(f"missing/empty {f.name}")
    text = log_path.read_text()
    if "Job finished" not in text:
        raise ValidationError("scph.log has no 'Job finished'")
    conv, notconv = [], []
    for l in text.splitlines():
        if "convergence achieved" in l and "Temp =" in l:
            conv.append(float(l.split("=")[1].split(":")[0]))
        elif "not converged" in l and "Temp =" in l:
            notconv.append(float(l.split("=")[1].split(":")[0]))
    crit = {}
    for key in ("TOL_SCPH", "MIXALPHA", "MAXITER", "KMESH_SCPH", "KMESH_INTERPOLATE", "SELF_OFFDIAG", "IALGO", "LOWER_TEMP"):
        for l in text.splitlines():
            if l.strip().startswith(key):
                crit[key] = l.split("=", 1)[1].strip()
                break
    p = job["params"]
    if crit.get("KMESH_SCPH", "").split()[:1] != [str(p["q1"])]:
        raise ValidationError(f"KMESH_SCPH in log = {crit.get('KMESH_SCPH')}")
    if p.get("mixalpha") is not None and float(crit.get("MIXALPHA", "nan").split()[0]) != float(p["mixalpha"]):
        raise ValidationError(f"MIXALPHA in log = {crit.get('MIXALPHA')}")
    above = sorted(t for t in conv if t >= PRODUCTION_TMIN)
    boundary = None
    for t in sorted(conv + notconv, reverse=True):  # contiguous converged range from the top
        if t in notconv:
            break
        boundary = t
    metrics = {"converged_T": sorted(conv), "not_converged_T": sorted(notconv), "criteria": crit,
               "lowest_contiguous_converged_T": boundary, "converged_T_ge_tmin": above}
    outs = [log_path, dfc2, bands, wd / "scph.in"]
    return metrics, outs


def validate_mode4ph(job, wd):
    p = job["params"]
    nT = round((p["tmax"] - p["tmin"]) / p["dt"]) + 1
    if "Job finished" not in (wd / "rta.log").read_text(errors="replace"):
        raise ValidationError("rta.log has no 'Job finished'")
    res, outs = {}, []
    for i, mode in enumerate(p["modes"], 1):
        f = wd / f"{job['id']}.Gamma.{i}"
        if not f.exists():
            raise ValidationError(f"missing {f.name}")
        lines = f.read_text().splitlines()
        omega = float([l for l in lines if l.startswith("# Frequency")][0].split("=")[1])
        rows = [l.split() for l in lines if l.strip() and not l.startswith("#")]
        if len(rows) != nT:
            raise ValidationError(f"{f.name}: {len(rows)} temperature rows, expected {nT}")
        vals = [[float(r[0]), float(r[1]), float(r[2])] for r in rows]
        if any(v != v for r in vals for v in r):
            raise ValidationError(f"{f.name}: NaN")
        res[f"{mode[0]},{mode[1]},{mode[2]},{mode[3]}"] = {"omega_cm1": omega,
            "T": [v[0] for v in vals], "Gamma3_hwhm": [v[1] / 2 for v in vals], "Gamma4_hwhm": [v[2] / 2 for v in vals]}
        outs.append(f)
    outs += [wd / x for x in ("rta.in", "rta.log", "modes.txt", "dfc2.log")]
    metrics = {"modes": res}
    audit = [l.split("AUDIT 4ph(c) identity set:")[1].strip() for l in
             (wd / "rta.log").read_text(errors="replace").splitlines() if "AUDIT 4ph(c) identity set:" in l]
    if p.get("anphon_tree"):  # the audit build must have reported its exclusion counters for every mode
        if len(audit) != len(p["modes"]):
            raise ValidationError(f"audit build: {len(audit)} AUDIT lines for {len(p['modes'])} modes")
        metrics["audit_4ph_c"] = audit
    return metrics, outs


VALIDATORS = {"rta": validate_rta, "scph": validate_scph, "mode4ph": validate_mode4ph, "evec": validate_evec}


# --------------------------------------------------------------------------- D0 decision
def _soft_freqs(scph_entry, T):
    """Gamma-TO1, M, R frequencies (cm^-1) at temperature T from the .scph_bands path G-X-M-G-R-M."""
    f = workdir_of(scph_entry) / f"{prefix_of(scph_entry)}.scph_bands"
    rows = [list(map(float, l.split())) for l in f.read_text().splitlines() if l.strip() and not l.startswith("#")]
    head = [l for l in f.read_text().splitlines() if l.startswith("#")][1].lstrip("#").split()
    kx = {"Gamma": float(head[0]), "M": float(head[2]), "R": float(head[4])}
    rows = [r for r in rows if abs(r[0] - T) < 1e-6]
    out = {}
    for lab, x in kx.items():
        r = min(rows, key=lambda r: abs(r[1] - x))[2:]
        out[lab] = r[3] if lab == "Gamma" else r[0]
    return out


def run_decision(entry, m):
    jobs = m["jobs"]
    rows, ok = [], True
    for T in (250, 300):
        a, b = jobs[f"rta_pm_q8_T{T}_m12"]["metrics"], jobs[f"rta_pm_q12_T{T}_m12"]["metrics"]
        fa, fb = _soft_freqs(jobs["scph_pm_q8"], T), _soft_freqs(jobs["scph_pm_q12"], T)
        dP = 100 * (b["kappa_P"] - a["kappa_P"]) / a["kappa_P"]
        dC = 100 * (b["kappa_C"] - a["kappa_C"]) / a["kappa_C"]
        dF = {k: 100 * (fb[k] - fa[k]) / fa[k] for k in fa}
        okT = abs(dP) < 3 and abs(dF["R"]) < 5 and abs(dF["M"]) < 5
        ok &= okT
        rows.append((T, a, b, dP, dC, fa, fb, dF, okT))
    q1 = 8 if ok else 12
    s8, s12 = jobs["scph_pm_q8"]["metrics"], jobs["scph_pm_q12"]["metrics"]
    md = [f"# D0 — SCPH q1-grid rule (applied automatically)\n",
          f"Generated by `queue/run_queue.py` at {now_iso()}.\n",
          "Rule: adopt q1 = 8³ if, at both 250 K and 300 K, going 8³ → 12³ moves κ_P by < 3% "
          "AND the R and M frequencies by < 5%; otherwise adopt q1 = 12³. "
          "Transport mesh fixed at 12³; PM on; tetrahedron; NONANALYTIC = 3; KAPPA_COHERENT = 2.\n",
          f"## Adopted: **q1 = {q1}³**\n",
          "| T (K) | q1 | κ_P | κ_C | Ω Γ-TO1 | Ω M | Ω R (cm⁻¹) |", "|---|---|---|---|---|---|---|"]
    for T, a, b, dP, dC, fa, fb, dF, okT in rows:
        md.append(f"| {T} | 8³ | {a['kappa_P']:.4f} | {a['kappa_C']:.4f} | {fa['Gamma']:.2f} | {fa['M']:.2f} | {fa['R']:.2f} |")
        md.append(f"| {T} | 12³ | {b['kappa_P']:.4f} | {b['kappa_C']:.4f} | {fb['Gamma']:.2f} | {fb['M']:.2f} | {fb['R']:.2f} |")
        md.append(f"| {T} | Δ (%) | {dP:+.2f} | {dC:+.2f} | {dF['Gamma']:+.2f} | {dF['M']:+.2f} | {dF['R']:+.2f} |")
    md.append("")
    for T, *_r, okT in rows:
        md.append(f"- {T} K: criteria {'met' if okT else 'NOT met'}.")
    md += ["", "## SCPH convergence record (criteria as recorded in scph.log)\n"]
    for name, s in (("8³", s8), ("12³", s12)):
        md.append(f"- q1 = {name}: criteria {s['criteria']}; converged T = {s['converged_T']}; "
                  f"not converged = {s['not_converged_T']}; lowest contiguous converged T = "
                  f"{s['lowest_contiguous_converged_T']} K.")
    md.append("")
    if not ok:
        md.append("Criteria not met (q1 8³→12³ moved κ_P by > 3% or R/M by > 5%); q1 = 12³ adopted.")
    DECISION_MD.write_text("\n".join(md) + "\n")
    # enqueue production for the adopted grid
    scph = jobs[f"scph_pm_q{q1}"]["metrics"]
    temps = [int(t) for t in scph["converged_T_ge_tmin"]]
    missing = [t for t in range(PRODUCTION_TMIN, 1001, 50) if t not in temps]
    added = []
    for T in temps:
        for mesh, pr in ((12, 30), (16, 40)):
            j = rta_job(q1, T, mesh, priority=pr)
            if j["id"] not in jobs:
                jobs[j["id"]] = new_entry(j); added.append(j["id"])
    for T, pr in ((300, 50), (800, 60)):
        for j in (rta_job(q1, T, 20, priority=pr), evec_job(q1, T, 20, priority=pr)):
            if j["id"] not in jobs:
                jobs[j["id"]] = new_entry(j); added.append(j["id"])
    for T in temps:
        for j in (evec_job(q1, T, 12), evec_job(q1, T, 16), evec_job(q1, T, 8, priority=65),
                  rta_job(q1, T, 8, priority=65)):
            if j["id"] not in jobs:
                jobs[j["id"]] = new_entry(j); added.append(j["id"])
    m["meta"]["production_q1"] = q1
    m["meta"]["production_temps"] = temps
    m["meta"]["production_missing_T_ge_tmin"] = missing
    entry["metrics"] = {"q1": q1, "criteria_met": ok, "enqueued": added, "missing_T": missing}
    return f"adopted q1={q1}; enqueued {len(added)} production jobs; missing T>={PRODUCTION_TMIN}: {missing}"


# --------------------------------------------------------------------------- status
def _ekey(j):
    return (j["kind"], j["params"].get("mesh") if j["kind"] in ("rta", "mode4ph", "evec") else j["params"].get("q1"))


def expected_s(job, jobs):
    key = _ekey(job)
    if job.get("est_s") and not any(j["status"] == "done" and _ekey(j) == key and not j.get("import") for j in jobs.values()):
        return job["est_s"]
    if job["kind"] == "decision":
        return 1
    done = [j["duration_s"] for j in jobs.values()
            if j["status"] == "done" and j.get("duration_s") and not j.get("import")
            and _ekey(j) == key]
    return statistics.median(done) if done else DEFAULT_EST.get(key, 3600)


def write_status(m, note=""):
    jobs = m["jobs"]
    by = {}
    for j in jobs.values():
        by.setdefault(j["status"], []).append(j)
    pend = sorted(by.get("pending", []), key=lambda j: (j["priority"], j["id"]))
    eta = sum(expected_s(j, jobs) for j in pend)
    L = [f"# Queue status — {now_iso()} ({HOST})", "", note, "",
         " | ".join(f"{k}: {len(v)}" for k, v in sorted(by.items())),
         f"Pending work (sum of expected durations): {eta / 3600:.1f} h", ""]
    if m["meta"].get("production_q1"):
        L.append(f"Production q1 = {m['meta']['production_q1']}³; temperatures {m['meta']['production_temps']}; "
                 f"missing ≥{PRODUCTION_TMIN} K: {m['meta']['production_missing_T_ge_tmin']}")
        L.append("")
    L += ["| job | status | attempts | duration | κ_P | κ_C | note |", "|---|---|---|---|---|---|---|"]
    for j in sorted(jobs.values(), key=lambda j: (j["priority"], j["id"])):
        mt = j.get("metrics", {})
        dur = f"{j['duration_s'] / 60:.1f} min" if j.get("duration_s") else ""
        kp = f"{mt['kappa_P']:.4f}" if "kappa_P" in mt else ""
        kc = f"{mt['kappa_C']:.4f}" if "kappa_C" in mt else ""
        nt = "imported" if j.get("import") else ""
        if j["kind"] == "scph" and mt:
            nt += f" lowest contiguous converged T = {mt.get('lowest_contiguous_converged_T')} K"
        L.append(f"| {j['id']} | {j['status']} | {j['attempts']} | {dur} | {kp} | {kc} | {nt} |")
    for j in by.get("failed", []):
        L += ["", f"## FAILED {j['id']}", "```", *j.get("last_log", []), "```"]
    STATUS.write_text("\n".join(L) + "\n")


# --------------------------------------------------------------------------- execution
def archive(job, outs):
    variant = job["params"].get("variant") or f"{'pm' if job['params']['pm'] else 'nopm'}_q{job['params']['q1']}"
    dest = ARCHIVE / variant / job["id"]
    dest.mkdir(parents=True, exist_ok=True)
    for f in outs:
        shutil.copy2(f, dest / Path(f).name)
    return str(dest.relative_to(REPO))


def finish_ok(job_id, metrics, outs, start, build):
    with Manifest() as m:
        j = m["jobs"][job_id]
        j.update({"status": "done", "end": now_iso(), "duration_s": time.time() - start, "metrics": metrics,
                  "outputs": {Path(f).name: sha256(f) for f in outs}, "build": build, "host": HOST})
        j["archive"] = archive(j, outs)
        j["history"].append({"t": now_iso(), "event": "done", "host": HOST})
        write_status(m, f"last: {job_id} done")


def tail(path, n=50):
    try:
        return Path(path).read_text(errors="replace").splitlines()[-n:]
    except Exception:
        return []


def finish_fail(job_id, reason, wd, final=False):
    log(f"failed {job_id}: {reason}")
    with Manifest() as m:
        j = m["jobs"][job_id]
        j["attempts"] += 1
        j["history"].append({"t": now_iso(), "event": "failed", "reason": reason, "host": HOST})
        logname = {"scph": "scph.log", "evec": "ph.log"}.get(j["kind"], "rta.log")
        j["last_log"] = [f"reason: {reason}"] + tail(wd / logname)
        j["status"] = "failed" if final or j["attempts"] > MAX_RETRIES else "pending"
        write_status(m, f"last: {job_id} failed ({reason})")
    (LOCKS / f"{job_id}.lock").unlink(missing_ok=True)


def requeue_killed(job_id, why):
    with Manifest() as m:
        j = m["jobs"][job_id]
        j["status"] = "pending"
        j["history"].append({"t": now_iso(), "event": f"killed:{why}", "host": HOST})
        write_status(m, f"{job_id} killed ({why}) and re-queued")
    (LOCKS / f"{job_id}.lock").unlink(missing_ok=True)


def kill_group(proc):
    try:
        os.killpg(proc.pid, signal.SIGTERM)
        for _ in range(20):
            if proc.poll() is not None:
                return
            time.sleep(0.5)
        os.killpg(proc.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass


def group_rss_kb(pgid):
    """Resident set size (KiB) summed over the process group pgid (the job's session)."""
    try:
        out = subprocess.run(["ps", "-A", "-o", "pgid=,rss="], capture_output=True, text=True).stdout
    except Exception:
        return 0
    return sum(int(r) for g, r in (l.split() for l in out.splitlines() if l.strip()) if int(g) == pgid)


def launcher_for(job):
    """The anphon invocation for this job (replaces the quoted production ANPHON path in the command)."""
    p = job["params"]
    binary = TOOL / p["anphon_tree"] / "_build/anphon/anphon" if p.get("anphon_tree") else ANPHON
    tok = f'"{binary}"'
    np_ = int(p.get("mpi_np") or 0)
    if np_ > 1:
        tok = f'"{ENV_BIN / "mpirun"}" -np {np_} {tok}'
    if p.get("measure_mem"):
        tok = f"/usr/bin/time -l -o time_l.txt {tok}"
    return tok, binary


def time_l_max_rss_gb(path):
    for line in Path(path).read_text().splitlines():
        if "maximum resident set size" in line:
            return int(line.split()[0]) / 2 ** 30  # macOS reports bytes
    return None


def pick(m, until):
    jobs = m["jobs"]
    cands = sorted((j for j in jobs.values() if j["status"] == "pending"), key=lambda j: (j["priority"], j["id"]))
    blocked_by_until = []
    for j in cands:
        if not all(jobs.get(d, {}).get("status") == "done" for d in j["deps"]):
            continue
        if until and time.time() + expected_s(j, jobs) > until:
            blocked_by_until.append(j["id"])
            continue
        return j, blocked_by_until
    return None, blocked_by_until


def run_one(job_id, until, build):
    with Manifest() as m:
        j = m["jobs"][job_id]
        if j["status"] != "pending":
            return "skip"
        j.update({"status": "running", "start": now_iso(), "host": HOST})
        j["history"].append({"t": now_iso(), "event": "start", "host": HOST, "attempt": j["attempts"] + 1})
        job = json.loads(json.dumps(j))
        snap = json.loads(json.dumps(m))
        write_status(m, f"running {job_id}")
    if job["kind"] == "decision":
        with Manifest() as m:
            msg = run_decision(m["jobs"][job_id], m)
            m["jobs"][job_id].update({"status": "done", "end": now_iso(), "duration_s": 0.0})
            write_status(m, msg)
        log(msg)
        return "done"
    wd = workdir_of(job)
    start = time.time()
    lock = LOCKS / f"{job_id}.lock"
    LOCKS.mkdir(exist_ok=True)
    if job.get("import"):  # register an already-computed unit: validate in place, no run
        try:
            metrics, outs = VALIDATORS[job["kind"]](job, wd)
        except Exception as exc:
            finish_fail(job_id, f"import validation: {exc}", wd)
            return "failed"
        finish_ok(job_id, metrics, outs, start, f"imported:{job['import'].get('build', 'same toolchain')}")
        return "done"
    if wd.exists():  # purge partial outputs of any earlier attempt
        shutil.rmtree(wd)
    wd.mkdir(parents=True)
    cmd = (write_scph_inputs(job, wd) if job["kind"] == "scph" else
           write_mode4ph_inputs(job, wd, snap) if job["kind"] == "mode4ph" else
           write_evec_inputs(job, wd, snap) if job["kind"] == "evec" else write_rta_inputs(job, wd, snap))
    env = dict(os.environ, OMP_NUM_THREADS=str(OMP_THREADS), PATH=f"{ENV_BIN}:{os.environ.get('PATH', '')}")
    np_ = int(job["params"].get("mpi_np") or 0)
    if np_ > 1:  # MPI-parallel anphon (e.g. the four-phonon sum is distributed over ranks only)
        env["OMP_NUM_THREADS"] = "1"
    tok, binary = launcher_for(job)
    if job["params"].get("anphon_tree"):
        if not binary.exists():
            finish_fail(job_id, f"audit binary missing: {binary}", wd, final=True)
            return "failed"
        build = f"AUDIT:{job['params']['anphon_tree']}/anphon-sha256:{sha256(binary)[:12]}/{HOST}"
    cmd = cmd.replace(f'"{ANPHON}"', tok)
    rss_kill = float(job["params"].get("rss_kill_gb") or 0)
    peak_kb = 0
    proc = subprocess.Popen(["/bin/bash", "-c", cmd], cwd=wd, env=env, start_new_session=True)
    log(f"started {job_id} (pid {proc.pid}) in {wd}" + (f" [mpirun -np {np_}]" if np_ > 1 else "")
        + (f" [{build}]" if job["params"].get("anphon_tree") else "")
        + (f" [rss watchdog {rss_kill:g} GB]" if rss_kill else ""))
    last_hb = 0
    while True:
        if time.time() - last_hb > HEARTBEAT_S:
            lock.write_text(json.dumps({"host": HOST, "runner_pid": os.getpid(), "child_pid": proc.pid,
                                        "heartbeat": time.time(), "job": job_id}))
            last_hb = time.time()
        rc = proc.poll()
        if rc is not None:
            break
        if rss_kill or job["params"].get("measure_mem"):
            peak_kb = max(peak_kb, group_rss_kb(proc.pid))
            if rss_kill and peak_kb > rss_kill * 2 ** 20:
                log(f"watchdog: {job_id} resident set {peak_kb / 2 ** 20:.1f} GB > {rss_kill:g} GB, killing")
                kill_group(proc)
                finish_fail(job_id, f"watchdog: RSS {peak_kb / 2 ** 20:.1f} GB > {rss_kill:g} GB", wd, final=True)
                return "failed"
        if STOP_NOW.exists():
            log(f"STOP_NOW: killing {job_id}")
            kill_group(proc)
            requeue_killed(job_id, "stop_now")
            return "stop_now"
        if until and time.time() >= until:
            log(f"--until time reached: killing {job_id}")
            kill_group(proc)
            requeue_killed(job_id, "until")
            return "until"
        time.sleep(POLL_S)
    if rc != 0:
        finish_fail(job_id, f"exit status {rc}", wd)
        return "failed"
    try:
        metrics, outs = VALIDATORS[job["kind"]](job, wd)
    except Exception as exc:
        finish_fail(job_id, f"validation: {exc}", wd)
        return "failed"
    if peak_kb:
        metrics["peak_rss_sampled_gb"] = peak_kb / 2 ** 20
    if (wd / "time_l.txt").exists():
        outs = list(outs) + [wd / "time_l.txt"]
        rss = time_l_max_rss_gb(wd / "time_l.txt")
        if rss is not None:
            metrics["max_rss_time_gb"] = rss
    finish_ok(job_id, metrics, outs, start, build)
    lock.unlink(missing_ok=True)
    log(f"done {job_id} in {(time.time() - start) / 60:.1f} min: "
        + ", ".join(f"{k}={v:.4g}" for k, v in metrics.items() if isinstance(v, float)))
    return "done"


def recover():
    """Re-queue jobs left 'running' by a dead runner (same host) or with a stale heartbeat."""
    with Manifest() as m:
        for j in m["jobs"].values():
            if j["status"] != "running":
                continue
            lock = LOCKS / f"{j['id']}.lock"
            info = json.loads(lock.read_text()) if lock.exists() else {}
            stale = time.time() - info.get("heartbeat", 0) > STALE_S
            dead = info.get("host", HOST) == HOST and not pid_alive(info.get("runner_pid", -1))
            if stale or dead:
                j["status"] = "pending"
                j["history"].append({"t": now_iso(), "event": "recovered:" + ("stale" if stale else "runner-dead")})
                lock.unlink(missing_ok=True)
                log(f"recovered {j['id']} -> pending ({'stale heartbeat' if stale else 'runner dead'})")
        write_status(m, "runner start")


def cmd_run(args):
    if PIDFILE.exists() and pid_alive(PIDFILE.read_text().strip()):
        sys.exit(f"another runner is active (pid {PIDFILE.read_text().strip()})")
    PIDFILE.write_text(str(os.getpid()))
    try:
        for f in (STOP, STOP_NOW):
            if f.exists():
                log(f"consuming stale {f.name} from a previous session")
                f.unlink()
        until = parse_until(args.until)
        if until:
            log(f"until {datetime.fromtimestamp(until, timezone.utc).isoformat()}")
        build = build_id()
        log(f"build {build}; OMP_NUM_THREADS={OMP_THREADS}")
        recover()
        code_mtime = Path(__file__).stat().st_mtime
        while True:
            if Path(__file__).stat().st_mtime != code_mtime:
                log("run_queue.py changed on disk: re-executing the runner with the new code")
                PIDFILE.unlink(missing_ok=True)
                os.execv(sys.executable, [sys.executable, __file__] + sys.argv[1:])
            if STOP.exists():
                log("STOP file present: stopping before the next job"); STOP.unlink(); break
            with Manifest() as m:
                job, blocked = pick(m, until)
                if job is None:
                    pend = [j["id"] for j in m["jobs"].values() if j["status"] == "pending"]
                    note = ("queue complete" if not pend else
                            f"no runnable job fits before the --until time ({len(blocked)} blocked)" if blocked else
                            f"{len(pend)} pending jobs blocked by unfinished/failed dependencies")
                    write_status(m, note)
            if job is None:
                log(note); break
            res = run_one(job["id"], until, build)
            if res in ("stop_now", "until"):
                STOP_NOW.unlink(missing_ok=True)
                break
    finally:
        PIDFILE.unlink(missing_ok=True)
        with Manifest() as m:
            write_status(m, "runner exited")
        log("runner exit")


# --------------------------------------------------------------------------- admin
def cmd_init(args):
    if MANIFEST.exists() and not args.force:
        sys.exit("manifest exists (use --force to rebuild)")
    runs = TOOL / "runs"
    jobs = [
        scph_job(8, imported_from={"dir": str(runs / "scph_on_k8"), "prefix": "scph_on_k8",
                                   "build": "2026-10-07 runs/scph_on_k8, same toolchain"}, priority=1),
        rta_job(8, 300, 12, imported_from={"dir": str(runs / "rta_on_k8_m12"), "prefix": "rta_on_k8_m12",
                                          "build": "2026-10-07 runs/rta_on_k8_m12, same toolchain"}, priority=2),
        rta_job(8, 300, 16, imported_from={"dir": str(runs / "rta_on_k8_m16"), "prefix": "rta_on_k8_m16",
                                          "build": "2026-10-07 runs/rta_on_k8_m16, same toolchain"}, priority=40),
        rta_job(8, 250, 12, priority=3),
        scph_job(12, priority=4),
        rta_job(12, 300, 12, priority=5),
        rta_job(12, 250, 12, priority=6),
    ]
    dec = {"id": "d0_decision", "kind": "decision", "priority": 7, "params": {},
           "deps": ["rta_pm_q8_T250_m12", "rta_pm_q8_T300_m12", "rta_pm_q12_T250_m12", "rta_pm_q12_T300_m12",
                    "scph_pm_q8", "scph_pm_q12"]}
    with Manifest() as m:
        m["jobs"] = {j["id"]: new_entry(j) for j in jobs + [dec]}
        m["meta"] = {"created": now_iso(), "host": HOST}
        write_status(m, "initialised")
    print(STATUS.read_text())


def cmd_status(args):
    with Manifest() as m:
        write_status(m, "status requested")
    print(STATUS.read_text())


def cmd_add(args):
    """Add jobs from a JSON list of full job specs (id, kind, priority, deps, params[, est_s, import])."""
    specs = json.loads(Path(args.spec).read_text())
    with Manifest() as m:
        added = []
        for j in specs:
            if j["id"] in m["jobs"]:
                continue
            for d in j["deps"]:
                if d not in m["jobs"]:
                    raise SystemExit(f"{j['id']}: unknown dependency {d}")
            m["jobs"][j["id"]] = new_entry(j); added.append(j["id"])
        write_status(m, f"added {len(added)} jobs from {args.spec}")
    print("added:", added)


def cmd_requeue(args):
    with Manifest() as m:
        j = m["jobs"][args.job]
        j.update({"status": "pending", "attempts": 0})
        j["history"].append({"t": now_iso(), "event": "manual-requeue"})
        write_status(m, f"{args.job} requeued")


def main():
    if "QUEUE_TOOLCHAIN" not in os.environ:
        sys.exit("set QUEUE_TOOLCHAIN to the toolchain directory (ALAMODE builds, env/, sto_example/, runs/)")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    r = sub.add_parser("run"); r.add_argument("--until"); r.add_argument("--resume", action="store_true")
    i = sub.add_parser("init"); i.add_argument("--force", action="store_true")
    sub.add_parser("status")
    q = sub.add_parser("requeue"); q.add_argument("job")
    a = sub.add_parser("add"); a.add_argument("spec")
    args = ap.parse_args()
    if args.cmd in (None, "run"):
        if args.cmd is None:
            args.until = None
        cmd_run(args)
    else:
        {"init": cmd_init, "status": cmd_status, "requeue": cmd_requeue, "add": cmd_add}[args.cmd](args)


if __name__ == "__main__":
    main()
