#!/usr/bin/env python3
"""A2: compare the rerun of the three imported benchmark runs with the imported values (gate G-A2).

Frequencies: extraction rule of scripts/v3/r_point_check.py (path G-X-M-G-R-M of .scph_bands; Gamma = 4th
branch at k = 0, X/M/R = lowest branch); every temperature with 'convergence achieved' in both runs. Also the
largest difference over all branches and k points at those temperatures.
Transport: kappa_P and kappa_C from the .kl / .kl_coherent files (300 K row) when present, else from the manifests.
Output: OUT_CSV (benchmark table A2_rerun.csv) (provenance header) and a summary on stdout.

Usage: a2_compare.py IMPORTED_SCPH_DIR RERUN_SCPH_DIR IMPORTED_MANIFEST RERUN_MANIFEST INCOMING_ARCHIVE_DIR
       [IMPORTED_RTA_M12_DIR IMPORTED_RTA_M16_DIR]
"""
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

TAB_300 = {"Gamma": 136.0, "R": 37.0, "M": 85.0}  # 2015 Table I, q1 = 8^3 (arXiv:1506.01781v1 p. 7)
KAPPA_FIG = 8.7  # 2015 Fig. 7, read from figure


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def bands(d):
    f = next(Path(d).glob("*.scph_bands"))
    txt = f.read_text().splitlines()
    head = [float(x) for x in txt[1].lstrip("#").split()]
    kx = {"Gamma": head[0], "X": head[1], "M": head[2], "R": head[4]}
    a = np.loadtxt(f)
    conv = {float(l.split("=")[1].split(":")[0]) for l in (Path(d) / "scph.log").read_text().splitlines()
            if "convergence achieved" in l}
    return f, a, kx, conv


def soft(a, kx, T):
    b = a[a[:, 0] == T]
    out = {}
    for lab, x in kx.items():
        row = b[np.argmin(np.abs(b[:, 1] - x))][2:]
        out[lab] = row[3] if lab == "Gamma" else row[0]
    return out


def kappa_from_files(d):
    d = Path(d)
    out = {}
    for f in d.glob("*.kl"):
        for l in f.read_text().splitlines():
            s = l.split()
            if s and not l.startswith("#") and abs(float(s[0]) - 300) < 1e-6:
                v = [float(x) for x in s[1:10]]
                out["kappa_P"] = (v[0] + v[4] + v[8]) / 3
    for f in d.glob("*.kl_coherent"):
        for l in f.read_text().splitlines():
            s = l.split()
            if s and not l.startswith("#") and abs(float(s[0]) - 300) < 1e-6:
                v = [float(x) for x in s[1:4]]  # xx, yy, zz
                out["kappa_C"] = sum(v) / 3
    return out


def main():
    imp_scph, rr_scph, imp_man, rr_man, arch = sys.argv[1:6]
    imp_rta = sys.argv[6:8]
    mi, mr = json.load(open(imp_man)), json.load(open(rr_man))
    rows, ok = [], True
    build = {j: mr["jobs"][j].get("build", "") for j in mr["jobs"]}
    # sha256 of pulled outputs against the node manifest
    for jid, j in mr["jobs"].items():
        if j.get("status") != "done":
            continue
        for rel, h in (j.get("outputs") or {}).items():
            p = Path(arch) / Path(rel).name if not (Path(arch) / rel).exists() else Path(arch) / rel
            cand = list(Path(arch).rglob(Path(rel).name))
            cand = [c for c in cand if jid in str(c)]
            got = sha(cand[0]) if cand else None
            print(f"sha256 {jid} {Path(rel).name}: {'OK' if got == h else 'MISMATCH/absent'}")
            ok &= got == h

    fi, ai, kxi, ci = bands(imp_scph)
    fr, ar, kxr, cr = bands(rr_scph)
    Ts = sorted(ci & cr)
    print(f"converged T imported {sorted(ci)}\nconverged T rerun    {sorted(cr)}")
    same_conv = sorted(ci) == sorted(cr)
    ok &= same_conv
    worst = 0.0
    for T in Ts:
        si, sr = soft(ai, kxi, T), soft(ar, kxr, T)
        for lab in ("Gamma", "M", "R"):
            d = sr[lab] - si[lab]
            worst = max(worst, abs(d))
            rows.append([f"Omega_{lab}_T{int(T)}_cm-1", si[lab], sr[lab], d, build.get("scph_pm_q8_rb", ""),
                         "scph_pm_q8 -> scph_pm_q8_rb", sha(fr)[:16]])
    full = max(np.abs(ai[ai[:, 0] == T][:, 2:] - ar[ar[:, 0] == T][:, 2:]).max() for T in Ts)
    print(f"largest |dOmega| at Gamma/M/R over converged T: {worst:.2e} cm^-1; all branches/k: {full:.2e} cm^-1")
    ok &= worst <= 1e-3
    s300 = soft(ar, kxr, 300.0)
    for lab in ("Gamma", "R", "M"):
        rows.append([f"Omega_{lab}_T300_vs_2015_TableI_cm-1", TAB_300[lab], s300[lab], s300[lab] - TAB_300[lab],
                     build.get("scph_pm_q8_rb", ""), "scph_pm_q8_rb vs Table I (stated)", ""])

    pairs = [("rta_pm_q8_T300_m12", "rta_pm_q8_rb_T300_m12", imp_rta[0] if imp_rta else None),
             ("rta_pm_q8_T300_m16", "rta_pm_q8_rb_T300_m16", imp_rta[1] if len(imp_rta) > 1 else None)]
    for ji, jr, idir in pairs:
        if mr["jobs"].get(jr, {}).get("status") != "done":
            print(f"{jr}: not done yet")
            continue
        ki = kappa_from_files(idir) if idir else {}
        rdir = [c for c in Path(arch).rglob(jr) if c.is_dir()]
        kr = kappa_from_files(rdir[0]) if rdir else {}
        for q in ("kappa_P", "kappa_C"):
            vi = ki.get(q, mi["jobs"][ji]["metrics"][q])
            vr = kr.get(q, mr["jobs"][jr]["metrics"][q])
            rel = (vr - vi) / vi
            ok &= abs(rel) <= 1e-3
            rows.append([f"{q}_{jr.split('_')[-1]}_W/mK", vi, vr, vr - vi, build.get(jr, ""), f"{ji} -> {jr}",
                         ""])
            print(f"{jr} {q}: imported {vi:.6g} rerun {vr:.6g} rel {rel:+.2e}")
            if q == "kappa_P" and jr.endswith("m12"):
                rows.append(["kappa_P_m12_vs_2015_Fig7_W/mK", KAPPA_FIG, vr, vr - KAPPA_FIG, build.get(jr, ""),
                             f"{jr} vs 8.7 (read from figure)", ""])
    out = Path(__file__).resolve().parents[1] / "benchmark/A2_rerun.csv"
    with open(out, "w", newline="") as f:
        f.write(f"# A2 rerun vs imported runs; script scripts/v4/a2_compare.py; imported {imp_scph}; rerun {rr_scph}\n")
        w = csv.writer(f)
        w.writerow(["quantity", "imported_or_reference", "rerun", "difference", "build_id", "job_id", "output_sha256_16"])
        w.writerows(rows)
    print(f"same lowest converged T: {same_conv}; wrote {out}")
    print("criterion so far:", "MET" if ok else "NOT MET")


if __name__ == "__main__":
    main()
