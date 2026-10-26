"""Largest three-phonon damping ratio η₃ = Γ/Ω from the per-mode widths of the 8³ (300–1000 K) and 24³ (300 K) RTA
runs, which have no eigenvector data and are therefore not part of the damping map (Sec. 4).
Inputs: data/raw/alamode_sto/revision/pm_q12/rta_pm_q12_T*_m{8,24}/*.result (zenodo:pm_q12_rta.tar.gz).
Output: data/processed/revision/eta_max_outside_map.csv (one row per run, modes with Ω > 1e-3 cm⁻¹).
"""
from pathlib import Path

import pandas as pd

from fekappa.alamode_io import read_result
from fekappa.provenance import stamp

ROOT = Path(__file__).resolve().parents[2]
rows = []
for m in (8, 24):
    for d in sorted((ROOT / "data/raw/alamode_sto/revision/pm_q12").glob(f"rta_pm_q12_T*_m{m}")):
        T = int(d.name.split("_T")[1].split("_")[0])
        if T < 300:
            continue
        r = read_result(d / f"{d.name}.result")
        r = r[r.omega_cm1 > 1e-3]
        rows.append({"mesh": m, "T_K": T, "job": d.name, "eta_max": float((r.gamma_hwhm_cm1 / r.omega_cm1).max())})
out = ROOT / "data/processed/revision/eta_max_outside_map.csv"
pd.DataFrame(rows).to_csv(out, index=False, float_format="%.6g")
stamp(out, __file__, ["zenodo:pm_q12_rta.tar.gz (rta_pm_q12_T*_m8, rta_pm_q12_T300_m24 .result)"])
print(pd.DataFrame(rows).groupby("mesh").eta_max.max())
