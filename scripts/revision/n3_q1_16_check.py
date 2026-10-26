"""N3: SCPH grid check q1 = 12^3 -> 16^3 (PM; same T grid and criteria).
Frequencies at Gamma-TO1, X, M, R (lowest) at every converged T (same extraction as d_q1_freq_shift.py), and the
D0 rule applied to 12^3 -> 16^3 (kappa_P < 3 %, R and M < 5 %) using RTA on the 12^3 transport mesh at 250/300/500/800 K.
Source: data/raw/alamode_sto/revision/pm_q12/scph_pm_q12, pm_q16/scph_pm_q16, and the rta_pm_q{12,16}_T*_m12 manifest entries.
Output: data/processed/revision/n3_q1_16_check.csv
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
REV = ROOT / "data/raw/alamode_sto/revision"
D = {12: REV / "pm_q12/scph_pm_q12", 16: REV / "pm_q16/scph_pm_q16"}
man = json.loads((ROOT / "data/provenance/manifest_final.json").read_text())["jobs"]


def load(q):
    f = D[q] / f"scph_pm_q{q}.scph_bands"
    txt = f.read_text().splitlines()
    head = [float(x) for x in txt[1].lstrip("#").split()]
    kx = {"Gamma": head[0], "X": head[1], "M": head[2], "R": head[4]}
    a = np.loadtxt(f)
    conv = {float(l.split("=")[1].split(":")[0]) for l in (D[q] / "scph.log").read_text().splitlines() if "convergence achieved" in l}
    out = {}
    for T in sorted(set(a[:, 0])):
        if T not in conv:
            continue
        b = a[a[:, 0] == T]
        for lab, x in kx.items():
            row = b[np.argmin(np.abs(b[:, 1] - x))][2:]
            out[(T, lab)] = row[3] if lab == "Gamma" else row[0]
    return out, sorted(conv)


a12, c12 = load(12)
a16, c16 = load(16)
print("converged T q12:", c12[0], "-", c12[-1], f"({len(c12)})", " q16:", c16[0], "-", c16[-1], f"({len(c16)})")
rows = []
for (T, p) in sorted(a12):
    if (T, p) in a16:
        rows.append({"T_K": int(T), "quantity": f"omega_{p}", "q1_12": a12[(T, p)], "q1_16": a16[(T, p)]})
for T in (250, 300, 500, 800):
    for k in ("kappa_P", "kappa_C"):
        rows.append({"T_K": T, "quantity": k, "q1_12": man[f"rta_pm_q12_T{T}_m12"]["metrics"][k],
                     "q1_16": man[f"rta_pm_q16_T{T}_m12"]["metrics"][k]})
df = pd.DataFrame(rows)
df["shift_pct"] = 100 * (df.q1_16 - df.q1_12) / df.q1_12
df.to_csv(ROOT / "data/processed/revision/n3_q1_16_check.csv", index=False, float_format="%.5g")
print(df[df.quantity.str.startswith("omega")].pivot(index="T_K", columns="quantity", values="shift_pct").round(2).to_string())
print(df[df.quantity.str.startswith("omega_R") | df.quantity.str.startswith("omega_M")].pivot(index="T_K", columns="quantity", values="q1_16").round(2).head(4).to_string())
print(df[df.quantity.str.startswith("kappa")].round(4).to_string(index=False))


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"n3_q1_16_check.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['data/provenance/manifest_final.json', 'zenodo:scph_q12.tar.gz (scph_pm_q12)', 'zenodo:scph_q16.tar.gz (scph_pm_q16)'])
