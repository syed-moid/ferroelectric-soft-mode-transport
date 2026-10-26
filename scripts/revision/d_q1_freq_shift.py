"""SCPH frequencies at Gamma-TO1, X, M, R (lowest) for q1 = 8^3 and 12^3 at every converged T,
and the 8^3 -> 12^3 shift. Source: .scph_bands of $FEKAPPA_RUNS/scph_on_k8 (q1 8^3, PM) and of the production SCPH job scph_pm_q12.
Path G-X-M-G-R-M; Gamma = lowest optic at k=0 approached along [100]; X, M, R = lowest branch.
Output: data/processed/revision/d_q1_freq_shift.csv
"""
from pathlib import Path
from fekappa.provenance import external_dir
import numpy as np, pandas as pd
R = external_dir("FEKAPPA_RUNS", "the pre-queue run directory (Zenodo prequeue_runs)")
Q12 = Path(__file__).resolve().parents[2] / "data/raw/alamode_sto/revision/pm_q12/scph_pm_q12"  # zenodo:scph_q12
F = {8: R / "scph_on_k8/scph_on_k8.scph_bands", 12: Q12 / "scph_pm_q12.scph_bands"}
LOG = {8: R / "scph_on_k8/scph.log", 12: Q12 / "scph.log"}
def load(q):
    txt = F[q].read_text().splitlines()
    head = [float(x) for x in txt[1].lstrip("#").split()]
    kx = {"Gamma": head[0], "X": head[1], "M": head[2], "R": head[4]}
    a = np.loadtxt(F[q]); conv = {float(l.split("=")[1].split(":")[0]) for l in LOG[q].read_text().splitlines() if "convergence achieved" in l}
    out = {}
    for T in sorted(set(a[:, 0])):
        if T not in conv: continue
        b = a[a[:, 0] == T]
        for lab, x in kx.items():
            row = b[np.argmin(np.abs(b[:, 1] - x))][2:]
            out[(T, lab)] = row[3] if lab == "Gamma" else row[0]
    return out
a8, a12 = load(8), load(12)
rows = [{"T_K": T, "point": p, "omega_q1_8_cm1": a8[(T, p)], "omega_q1_12_cm1": a12[(T, p)],
         "shift_pct": 100 * (a12[(T, p)] - a8[(T, p)]) / a8[(T, p)]} for (T, p) in sorted(a8) if (T, p) in a12]
df = pd.DataFrame(rows); df.to_csv(Path(__file__).resolve().parents[2] / "data/processed/revision/d_q1_freq_shift.csv", index=False, float_format="%.4g")
print(df.pivot(index="T_K", columns="point", values="shift_pct").round(2).to_string())
print(df[df.point.isin(["R", "M"])].pivot(index="T_K", columns="point", values="omega_q1_12_cm1").round(2).T.to_string())


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"d_q1_freq_shift.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['$FEKAPPA_RUNS/scph_on_k8 (zenodo:prequeue_runs.tar.gz)', 'zenodo:scph_q12.tar.gz (scph_pm_q12)'])
