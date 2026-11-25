#!/usr/bin/env python3
"""Soft-mode table from SCPH runs, with the extraction rule of scripts/v3/r_point_check.py.

The rule: path G-X-M-G-R-M of .scph_bands; Gamma = 4th branch at k = 0; X, M, R = lowest branch.
Only temperatures with 'convergence achieved' in scph.log count as converged.

Mean-field fit per run: Omega_R^2 = A (T - T_c) over the converged temperatures in 300-600 K.
Omega^2 is stored as sign*Omega^2, since anphon prints imaginary modes as negative.

Usage: scph_softmodes.py OUT.csv LABEL=SCPH_DIR [LABEL=SCPH_DIR ...]
Prints the lowest contiguous converged T, the 300 K values and the fit per run.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def load(d):
    f = next(Path(d).glob("*.scph_bands"))
    txt = f.read_text().splitlines()
    head = [float(x) for x in txt[1].lstrip("#").split()]
    kx = {"Gamma": head[0], "X": head[1], "M": head[2], "R": head[4]}
    a = np.loadtxt(f)
    log = (Path(d) / "scph.log").read_text().splitlines()
    conv = {float(l.split("=")[1].split(":")[0]) for l in log if "convergence achieved" in l and "Temp =" in l}
    notc = {float(l.split("=")[1].split(":")[0]) for l in log if "not converged" in l and "Temp =" in l}
    rows = []
    for T in sorted(set(a[:, 0])):
        b = a[a[:, 0] == T]
        r = {"T_K": T, "converged": T in conv}
        for lab, x in kx.items():
            row = b[np.argmin(np.abs(b[:, 1] - x))][2:]
            r[lab] = row[3] if lab == "Gamma" else row[0]
        rows.append(r)
    df = pd.DataFrame(rows)
    boundary = None
    for T in sorted(conv | notc, reverse=True):
        if T in notc:
            break
        boundary = T
    return df, boundary


def main():
    out, specs = sys.argv[1], [a.split("=", 1) for a in sys.argv[2:]]
    frames = []
    for lab, d in specs:
        df, low = load(d)
        for k in ("Gamma", "X", "M", "R"):
            df[f"{k}2"] = np.sign(df[k]) * df[k] ** 2
        df.insert(0, "run", lab)
        df["lowest_contiguous_converged_T"] = low
        fit = df[df.converged & (df.T_K >= 300) & (df.T_K <= 600)]
        if len(fit) >= 2:
            A, B = np.polyfit(fit.T_K, fit.R2, 1)
            df["mf_A_cm2_per_K"], df["mf_Tc_K"], df["mf_npts"] = A, -B / A, len(fit)
        else:
            df["mf_A_cm2_per_K"], df["mf_Tc_K"], df["mf_npts"] = np.nan, np.nan, len(fit)
        frames.append(df)
        r = df[df.T_K == 300].iloc[0]
        print(f"{lab}: lowest contiguous converged T = {low} K; 300 K converged={r.converged}; "
              f"Gamma {r.Gamma:.2f} R {r.R:.2f} M {r.M:.2f} cm^-1; "
              f"mean-field A={df.mf_A_cm2_per_K.iloc[0]:.3f} cm^-2/K, Tc={df.mf_Tc_K.iloc[0]:.1f} K "
              f"({int(df.mf_npts.iloc[0])} pts)")
    pd.concat(frames).to_csv(out, index=False, float_format="%.6g")


if __name__ == "__main__":
    main()
