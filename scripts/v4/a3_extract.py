#!/usr/bin/env python3
"""A3: Omega and Omega^2 at Gamma, R, M (300 K) from the A3 SCPH runs; stop rule; spreads (reading table).

Extraction rule as scripts/v3/r_point_check.py and scripts/v4/a2_compare.py: path G-X-M-G-R-M of
.scph_bands; Gamma = 4th branch at k = 0, X/M/R = lowest branch. Omega^2 = sign(Omega) * Omega^2 (anphon writes
imaginary frequencies as negative numbers). A value is used only if 'convergence achieved' at 300 K in scph.log.

Usage: a3_extract.py ARCHIVE_ROOT A345_Q8_REFERENCE_DIR
  ARCHIVE_ROOT: pulled node archive (searched recursively for scph_rf_<fit>_pm_q<n> directories)
  A345_Q8_REFERENCE_DIR: scph_pm_q8_rb (A2 rerun = distributed file; reference for a345 at 8^3)
Output: OUT_CSV (benchmark table A3_frequencies.csv) and a summary on stdout. Verifies sha256 against the manifest
of each queue if queue manifests are found under ARCHIVE_ROOT/../queues.
"""
import csv
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np

FITS = ["a198", "a262", "a345", "a455", "a600"]
PROBES = ["a345_pic1", "rr1053"]
REF_2015 = {2: 69.0, 8: 37.0}  # Omega_R(300 K), 2015 paper, q1 = 2^3 and 8^3
T = 300.0


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read_run(d):
    d = Path(d)
    f = next(d.glob("*.scph_bands"))
    txt = f.read_text().splitlines()
    head = [float(x) for x in txt[1].lstrip("#").split()]
    kx = {"Gamma": head[0], "X": head[1], "M": head[2], "R": head[4]}
    a = np.loadtxt(f)
    conv = sorted({float(l.split("=")[1].split(":")[0]) for l in (d / "scph.log").read_text().splitlines()
                   if "convergence achieved" in l})
    out = {"file": f, "conv": conv, "lowest_conv": min(conv) if conv else None, "conv300": T in conv}
    b = a[a[:, 0] == T]
    for lab in ("Gamma", "R", "M"):
        row = b[np.argmin(np.abs(b[:, 1] - kx[lab]))][2:]
        w = row[3] if lab == "Gamma" else row[0]
        out[lab] = w
        out[lab + "2"] = np.sign(w) * w * w
    return out


def main():
    root, ref_dir = Path(sys.argv[1]), Path(sys.argv[2])
    runs = {}
    for d in sorted(root.rglob("scph_rf_*_pm_q*")):
        if not d.is_dir() or not list(d.glob("*.scph_bands")):
            continue
        m = re.fullmatch(r"scph_rf_(.+)_pm_q(\d+)", d.name)
        runs[(m.group(1), int(m.group(2)))] = read_run(d)
    runs[("a345", 8)] = read_run(ref_dir) | {"ref": True}

    # sha256 against queue manifests, if pulled
    qroot = root.parent / "queues"
    for man in sorted(qroot.glob("queue_a3_*/manifest.json")):
        for jid, j in json.load(open(man))["jobs"].items():
            if j.get("status") != "done":
                continue
            for rel, h in (j.get("outputs") or {}).items():
                cand = [c for c in root.rglob(Path(rel).name) if jid in str(c)]
                ok = bool(cand) and sha(cand[0]) == h
                if not ok:
                    print(f"sha256 MISMATCH/absent: {jid} {Path(rel).name}")
        print(f"sha256 checked: {man.parent.name}")

    rows = []
    for (fit, q), r in sorted(runs.items(), key=lambda kv: ((FITS + PROBES).index(kv[0][0])
                                                                if kv[0][0] in FITS + PROBES else 99, kv[0][1])):
        src = "A2 rerun scph_pm_q8_rb (distributed file)" if r.get("ref") else f"scph_rf_{fit}_pm_q{q}"
        rows.append([fit, f"{q}^3", r["conv300"], r["lowest_conv"], f"{r['Gamma']:.5f}", f"{r['R']:.5f}",
                     f"{r['M']:.5f}", f"{r['Gamma2']:.2f}", f"{r['R2']:.2f}", f"{r['M2']:.2f}",
                     f"{r['R'] - REF_2015.get(q, np.nan):.5f}" if q in REF_2015 else "", src, sha(r["file"])[:16]])
        print(f"{fit:10s} q{q}: conv300={r['conv300']} lowest={r['lowest_conv']}  "
              f"G {r['Gamma']:9.4f}  R {r['R']:9.4f}  M {r['M']:9.4f}")

    out = Path(__file__).resolve().parents[1] / "benchmark/A3_frequencies.csv"
    with open(out, "w", newline="") as f:
        f.write("# A3 SCPH (polarization mixing), 300 K. Omega in cm^-1; Omega^2 = sign*Omega^2 in cm^-2. "
                "Script scripts/v4/a3_extract.py. a345 at 8^3 = A2 rerun (distributed file).\n")
        w = csv.writer(f)
        w.writerow(["fit", "q1", "converged_300K", "lowest_converged_T", "Omega_Gamma", "Omega_R", "Omega_M",
                    "Omega2_Gamma", "Omega2_R", "Omega2_M", "Omega_R_minus_2015", "source", "bands_sha256_16"])
        w.writerows(rows)
    print(f"wrote {out}")

    # stop rule
    if ("a345_pic1", 8) in runs:
        p, a = runs[("a345_pic1", 8)], runs[("a345", 8)]
        d = p["R"] - a["R"]
        print(f"STOP RULE: Omega_R(8^3) a345_pic1 {p['R']:.4f} - a345 {a['R']:.4f} = {d:+.4f} cm^-1 -> "
              f"{'STOP (> 1 cm^-1)' if abs(d) > 1 or not p['conv300'] else 'continue (<= 1 cm^-1)'}"
              f"{'' if p['conv300'] else ' [300 K not converged for a345_pic1]'}")

    # spreads over the five accepted fits (reading table)
    for q in (2, 8):
        have = [f for f in FITS if (f, q) in runs and runs[(f, q)]["conv300"]]
        if len(have) == len(FITS):
            v = [runs[(f, q)]["R"] for f in FITS]
            print(f"spread Omega_R {q}^3: {max(v) - min(v):.4f} cm^-1 (min {min(v):.4f}, max {max(v):.4f}); "
                  f"closest to {REF_2015[q]}: {min(abs(x - REF_2015[q]) for x in v):.4f}")
        else:
            print(f"spread {q}^3: incomplete ({len(have)}/5 fits converged at 300 K)")


if __name__ == "__main__":
    main()
