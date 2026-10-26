"""Harmonic transfer cross-check: ALAMODE example FC2 against an independent PBEsol DFPT dispersion.

ALAMODE side (this repo's toolchain, `MODE = phonons`, NONANALYTIC = 3, same path):
  STO222.xml      harmonic-only fit (2x2x2) — the base the dfc2 renormalization is added to (production FC2)
  STO_anharm.xml  harmonic block of the anharmonic fit — the bare-control FC2 and the SCPH reference D_HA
  cell 7.363 bohr = 3.8963 Angstrom (PBEsol, ALAMODE example)
DFPT side (SV repo, read-only): dft/qe/SrTiO3/dispersion_pbesol/SrTiO3_pbesol.freq — QE ph.x 4x4x4
  q-grid, q2r (zasr crystal), matdyn band path G-X-M-G-R-X, 40 points/segment, asr crystal; relaxed
  PBEsol cell a = 3.8930 Angstrom (0.085% smaller). Gamma, X, M, R are exact DFPT grid points.
At Gamma both sides approach along [100] (non-analytic term included -> LO split).
Branches are compared in ascending order (imaginary modes as negative numbers).
Output: data/processed/revision/a2_harmonic_xcheck.csv
"""

from pathlib import Path

from fekappa.provenance import external_dir
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RUN = external_dir("FEKAPPA_RUNS", "the pre-queue run directory (Zenodo prequeue_runs)") / "harmonic_xcheck"
DFPT = external_dir("FEKAPPA_DFPT", "the root of the lattice-shear-viscosity repository (DFT dataset)") / "dft/qe/SrTiO3/dispersion_pbesol/SrTiO3_pbesol.freq"
POINTS = {"Gamma": (0.0, 0), "X": (0.426673, 40), "M": (0.853346, 80), "R": (2.195772, 160)}


def read_alamode(path):
    a = np.loadtxt(path)
    out = {}
    for lab, (k, _) in POINTS.items():
        out[lab] = np.sort(a[np.argmin(np.abs(a[:, 0] - k))][1:])
    return out


def read_qe(path):
    tok = path.read_text().split("/", 1)[1].split()
    nb = 15
    vals = np.array(tok, dtype=float).reshape(-1, 3 + nb)
    return {lab: np.sort(vals[i, 3:]) for lab, (_, i) in POINTS.items()}


qe = read_qe(DFPT)
rows = []
for src in ("STO222", "STO_anharm"):
    al = read_alamode(RUN / f"harm_{src}.bands")
    for lab in POINTS:
        for b in range(15):
            w_al, w_qe = al[lab][b], qe[lab][b]
            rows.append({"alamode_fc2": src, "point": lab, "branch": b + 1, "omega_alamode_cm1": w_al,
                         "omega_dfpt_pbesol_cm1": w_qe, "diff_cm1": w_al - w_qe,
                         "diff_pct": 100 * (w_al - w_qe) / abs(w_qe) if abs(w_qe) > 1 else np.nan})
df = pd.DataFrame(rows)
df.to_csv(ROOT / "data/processed/revision/a2_harmonic_xcheck.csv", index=False, float_format="%.3f")
pd.set_option("display.width", 200)
for src in ("STO222", "STO_anharm"):
    d = df[df.alamode_fc2 == src]
    print(f"\n== {src}")
    print(d.pivot(index="branch", columns="point", values="omega_alamode_cm1").round(1).join(
        d.pivot(index="branch", columns="point", values="omega_dfpt_pbesol_cm1").round(1), rsuffix="_DFPT").to_string())
    st = d[d.omega_dfpt_pbesol_cm1 > 50]
    print("stable branches (DFPT > 50 cm-1): median |diff| %.1f cm-1, max |diff| %.1f cm-1, median |diff%%| %.2f, max |diff%%| %.2f"
          % (st.diff_cm1.abs().median(), st.diff_cm1.abs().max(), st.diff_pct.abs().median(), st.diff_pct.abs().max()))
s1 = read_alamode(RUN / "harm_STO222.bands"); s2 = read_alamode(RUN / "harm_STO_anharm.bands")
print("\nSTO222 vs STO_anharm harmonic block, max |diff| per point:",
      {k: round(float(np.max(np.abs(s1[k] - s2[k]))), 3) for k in POINTS})


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"a2_harmonic_xcheck.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['$FEKAPPA_RUNS/harmonic_xcheck (zenodo:prequeue_runs.tar.gz)', '$FEKAPPA_DFPT (PBEsol DFPT dispersion, lattice-shear-viscosity repository)'])
