"""Unit tests for queue/run_queue.py (validation, manifest atomicity, recovery, retries).

Run: QUEUE_DIR=<tmp> /usr/bin/python3 queue/test_queue.py   (uses a scratch queue dir; reads the
existing 2026-10-07 run directories under $QUEUE_TOOLCHAIN/runs read-only, copying them
before corrupting anything).
"""

import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

TMP = Path(tempfile.mkdtemp(prefix="qtest_"))
os.environ["QUEUE_DIR"] = str(TMP / "queue")
sys.path.insert(0, str(Path(__file__).resolve().parent))
if "QUEUE_TOOLCHAIN" not in os.environ:
    sys.exit("set QUEUE_TOOLCHAIN to the toolchain directory that holds runs/ (the stored test runs)")
import run_queue as rq  # noqa: E402

RUNS = rq.TOOL / "runs"
results = []


def check(name, cond, detail=""):
    results.append((name, bool(cond)))
    print(f"{'PASS' if cond else 'FAIL'}  {name}  {detail}")


def rta_imp(d, prefix, T=300, mesh=12):
    j = rq.rta_job(8, T, mesh, imported_from={"dir": str(d), "prefix": prefix})
    return rq.new_entry(j)


# 1. a good RTA unit validates and reproduces the recorded numbers
good = RUNS / "rta_on_k8_m12"
m, outs = rq.validate_rta(rta_imp(good, "rta_on_k8_m12"), good)
check("rta validates", abs(m["kappa_P"] - 9.2384) < 1e-4 and abs(m["kappa_C"] - 0.2790) < 1e-4,
      f"kP={m['kappa_P']:.4f} kC={m['kappa_C']:.4f} closure={m['closure_rel']:.2e} n_irred={m['n_irred']}")

# 2. corrupted copies fail
for what in ("no_kl", "truncated_kc", "truncated_result", "wrong_T"):
    d = TMP / what
    shutil.copytree(good, d, symlinks=True)
    pre = d / "rta_on_k8_m12"
    if what == "no_kl":
        Path(f"{pre}.kl").unlink()
    elif what == "truncated_kc":
        lines = Path(f"{pre}.kc_elem").read_text().splitlines()
        Path(f"{pre}.kc_elem").write_text("\n".join(lines[: len(lines) // 2]) + "\n")
    elif what == "truncated_result":
        lines = Path(f"{pre}.result").read_text().splitlines()
        Path(f"{pre}.result").write_text("\n".join(lines[: len(lines) * 2 // 3]) + "\n")
    try:
        rq.validate_rta(rta_imp(d, "rta_on_k8_m12", T=250 if what == "wrong_T" else 300), d)
        check(f"corrupt[{what}] rejected", False)
    except Exception as exc:
        check(f"corrupt[{what}] rejected", True, f"({type(exc).__name__}: {str(exc)[:70]})")

# 3. SCPH validation records criteria and the convergence boundary
sj = rq.new_entry(rq.scph_job(8, imported_from={"dir": str(RUNS / "scph_on_k8"), "prefix": "scph_on_k8"}))
sm, _ = rq.validate_scph(sj, RUNS / "scph_on_k8")
check("scph validates", sm["lowest_contiguous_converged_T"] == 250.0 and sm["criteria"]["TOL_SCPH"] == "1e-10",
      f"boundary={sm['lowest_contiguous_converged_T']} notconv={sm['not_converged_T']} crit={sm['criteria']}")

# 4. input writers substitute every field (and refuse silently-unchanged templates)
wd = TMP / "inp_rta"; wd.mkdir()
man = {"jobs": {"scph_pm_q8": sj}}
cmd = rq.write_rta_inputs(rq.new_entry(rq.rta_job(8, 450, 16)), wd, man)
s = (wd / "rta.in").read_text()
check("rta.in written", "TMIN = 450; TMAX = 450" in s and "\n16 16 16\n" in s and "FC2XML = renorm.xml" in s
      and "scph_on_k8.scph_dfc2" in cmd and " 450 " in cmd)
wd2 = TMP / "inp_scph"; wd2.mkdir()
rq.write_scph_inputs(rq.new_entry(rq.scph_job(12)), wd2)
s2 = (wd2 / "scph.in").read_text()
check("scph.in written", "KMESH_SCPH = 12 12 12" in s2 and "SELF_OFFDIAG = 1" in s2 and "KMESH_INTERPOLATE = 2 2 2" in s2)
check("scph.in default MIXALPHA unchanged", "MIXALPHA = 0.2" in s2)
wd3 = TMP / "inp_scph_mix"; wd3.mkdir()
jm = rq.new_entry(rq.scph_job(8)); jm["params"]["mixalpha"] = 0.1
rq.write_scph_inputs(jm, wd3)
check("scph.in mixalpha", "MIXALPHA = 0.1" in (wd3 / "scph.in").read_text())
jv = rq.new_entry(rq.scph_job(8, imported_from={"dir": str(RUNS / "scph_on_k8"), "prefix": "scph_on_k8"}))
jv["params"]["mixalpha"] = 0.1
try:
    rq.validate_scph(jv, RUNS / "scph_on_k8"); check("scph mixalpha mismatch rejected", False)
except rq.ValidationError as exc:
    check("scph mixalpha mismatch rejected", "MIXALPHA" in str(exc), str(exc))

# 5. manifest atomic write + recovery of a job whose runner died / heartbeat stale
rq.QDIR.mkdir(parents=True); rq.LOCKS.mkdir()
with rq.Manifest() as mm:
    a = rq.new_entry(rq.rta_job(8, 500, 12)); a["status"] = "running"
    b = rq.new_entry(rq.rta_job(8, 550, 12)); b["status"] = "running"
    c = rq.new_entry(rq.rta_job(8, 600, 12)); c["status"] = "running"
    mm["jobs"] = {a["id"]: a, b["id"]: b, c["id"]: c}
(rq.LOCKS / f"{a['id']}.lock").write_text(json.dumps({"host": rq.HOST, "runner_pid": 999999, "heartbeat": __import__('time').time()}))
(rq.LOCKS / f"{b['id']}.lock").write_text(json.dumps({"host": "other-mac", "runner_pid": 1, "heartbeat": 0}))
(rq.LOCKS / f"{c['id']}.lock").write_text(json.dumps({"host": "other-mac", "runner_pid": 1, "heartbeat": __import__('time').time()}))
rq.recover()
st = json.loads(rq.MANIFEST.read_text())["jobs"]
check("recover: dead runner -> pending", st[a["id"]]["status"] == "pending")
check("recover: stale heartbeat -> pending", st[b["id"]]["status"] == "pending")
check("recover: live other-host job untouched", st[c["id"]]["status"] == "running")
check("no temp manifest left", not rq.MANIFEST.with_suffix(".json.tmp").exists())

# 6. retries: two failures -> pending, third -> failed
with rq.Manifest() as mm:
    mm["jobs"][a["id"]]["status"] = "running"
for k in range(3):
    rq.finish_fail(a["id"], f"test failure {k}", TMP)
st = json.loads(rq.MANIFEST.read_text())["jobs"][a["id"]]
check("retry policy (max 2 retries)", st["status"] == "failed" and st["attempts"] == 3)

# 7. --until parsing
ts = rq.parse_until("2026-10-12T12:00:00 UTC")
check("--until parse", abs(ts - 1791806400) < 1, f"{ts}")

# 8. launcher composition (optional per-job params; absent = the production binary, unchanged)
base = {"params": {}}
tok, b = rq.launcher_for(base)
check("launcher default unchanged", tok == f'"{rq.ANPHON}"' and b == rq.ANPHON, tok)
tok, _ = rq.launcher_for({"params": {"mpi_np": 16}})
check("launcher mpi", tok == f'"{rq.ENV_BIN / "mpirun"}" -np 16 "{rq.ANPHON}"', tok)
tok, b = rq.launcher_for({"params": {"anphon_tree": "alamode_4ph_audit", "mpi_np": 16}})
check("launcher audit tree", b == rq.TOOL / "alamode_4ph_audit/_build/anphon/anphon" and str(b) in tok
      and str(rq.ANPHON) not in tok, tok)
tok, _ = rq.launcher_for({"params": {"measure_mem": True}})
check("launcher time -l", tok == f'/usr/bin/time -l -o time_l.txt "{rq.ANPHON}"', tok)
(TMP / "time_l.txt").write_text("        1.00 real\n  3221225472  maximum resident set size\n")
check("time -l parse", abs(rq.time_l_max_rss_gb(TMP / "time_l.txt") - 3.0) < 1e-9)

# 9. watchdog RSS sampling sees a ~300 MB process group
import subprocess as sp, time as _t
pr = sp.Popen([sys.executable, "-c", "import time; a = bytearray(300 << 20); a[::4096] = b'x' * len(a[::4096]); time.sleep(20)"],
              start_new_session=True)
_t.sleep(3)
kb = rq.group_rss_kb(pr.pid)
rq.kill_group(pr)
check("watchdog rss sampling", kb > 250 * 1024, f"{kb / 1024:.0f} MB")

shutil.rmtree(TMP)
n_fail = sum(1 for _, ok in results if not ok)
print(f"\n{len(results) - n_fail}/{len(results)} passed")
sys.exit(1 if n_fail else 0)
