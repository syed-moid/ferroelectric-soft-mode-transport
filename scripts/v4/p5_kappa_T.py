#!/usr/bin/env python3
"""P5: kappa(T) of case e (archived fit, NONANALYTIC = 2, Born data of the 2015 study, q1 8^3,
12^3 transport) against the numerical kappa(T) curve of the 2015 study (two columns: T, kappa). Reads pulled run directories; writes a CSV and prints a table.

Per run: kappa_P = xx entry of the .kl row at T; kappa_C = xx entry of .kl_coherent; trace/3 and the largest
off-diagonal entry are recorded as well. Frequencies: from the "##Phonon Frequency" block of .result; the number of
modes with Omega < 0 or |Omega| < 1e-6 cm^-1 is counted, with the three acoustic modes at Gamma (point 1, branches
1-3) excluded. The 300 K row is the existing run e-kappa (p15_e_q8_T300_m12, P2). No setting is chosen here.

Usage: p5_kappa_T.py INCOMING KAPPA_T_2015 OUT_CSV
  INCOMING: directory with the run directories p15_e_q8_T{T}_m12; KAPPA_T_2015: the 2015 curve (T, kappa).
"""
import csv
import sys
from pathlib import Path

INC, REF, OUT = (Path(a) for a in sys.argv[1:4])
TEMPS = (200, 300, 400, 500, 600, 700, 800, 900, 1000)


def kl_row(path, tk):
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        v = [float(x) for x in line.split()]
        if abs(v[0] - tk) < 1e-6:
            return v[1:]
    raise ValueError(f"{path}: no row at {tk} K")


def freq_check(result):
    lines = result.read_text().splitlines()
    i = lines.index("##Phonon Frequency") + 2
    n_neg = n_zero = n = 0
    omin = None
    while not lines[i].startswith("##END"):
        k, b, w = lines[i].split()
        k, b, w = int(k), int(b), float(w)
        if not (k == 1 and b <= 3):
            n += 1
            n_neg += w < 0
            n_zero += abs(w) < 1e-6
            omin = w if omin is None or w < omin else omin
        i += 1
    return n, n_neg, n_zero, omin


def main():
    ref = {}
    for line in REF.read_text().split("\n"):
        if line.strip():
            t, k = line.split()
            ref[round(float(t))] = float(k)
    rows = []
    for tk in TEMPS:
        job = f"p15_e_q8_T{tk}_m12"
        d = INC / job
        row = {"T_K": tk, "job": job, "kappa_2015": ref[tk]}
        kl, klc, res = d / f"{job}.kl", d / f"{job}.kl_coherent", d / f"{job}.result"
        if not (kl.exists() and klc.exists() and res.exists()):
            row["status"] = "not available"
            rows.append(row)
            continue
        v = kl_row(kl, tk)
        c = kl_row(klc, tk)
        n, n_neg, n_zero, omin = freq_check(res)
        kp = v[0]
        row.update(kappa_P_xx=kp, kappa_P_trace3=(v[0] + v[4] + v[8]) / 3,
                   kappa_P_offdiag_max=max(abs(v[j]) for j in (1, 2, 3, 5, 6, 7)),
                   kappa_C_xx=c[0], diff_abs=kp - ref[tk], diff_rel_pct=100 * (kp - ref[tk]) / ref[tk],
                   modes_checked=n, n_negative=n_neg, n_zero=n_zero, omega_min=omin,
                   status="existing run e-kappa (P2)" if tk == 300 else "new (P5)")
        rows.append(row)
    keys = ["T_K", "job", "status", "kappa_P_xx", "kappa_2015", "diff_abs", "diff_rel_pct", "kappa_C_xx",
            "kappa_P_trace3", "kappa_P_offdiag_max", "modes_checked", "n_negative", "n_zero", "omega_min"]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as f:
        f.write(f"# script: scripts/v4/p5_kappa_T.py; 2015 column from {REF.name} (2015 study)\n")
        w = csv.DictWriter(f, keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        if "kappa_P_xx" in r:
            print(f"{r['T_K']:5d} {r['kappa_P_xx']:8.4f} {r['kappa_2015']:9.6f} {r['diff_abs']:+8.4f} "
                  f"{r['diff_rel_pct']:+7.3f} %  kC {r['kappa_C_xx']:.4f}  neg {r['n_negative']} zero {r['n_zero']} "
                  f"min {r['omega_min']:.3f}  {r['status']}")
        else:
            print(f"{r['T_K']:5d} not available")
    print("written:", OUT)


if __name__ == "__main__":
    main()
