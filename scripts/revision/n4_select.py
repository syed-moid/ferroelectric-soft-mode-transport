"""N4 heat-carrier selection.
Two independent descriptive samples for the diagram-(c) audit build (external q must lie on the k-mesh):
  8^3 sample: q on the 8^3 mesh (8^3 is a subset of 16^3), ranked from the 16^3 production table;
  6^3 sample: q on the 6^3 mesh (6^3 is a subset of 12^3), ranked from the 12^3 production table.
Rule (each sample): 300 K; per-mode kappa_P contribution = kappa_mode_W_mK (already multiplicity-weighted: its sum
over the table is kappa_P of that mesh); classes by the character rule (sector_char A = acoustic, soft = low optic,
hard = hard optic); the 2 largest per class, at most one mode per q-point within a sample; ties broken by lower
q_index then lower branch. The modes are identified by (q, branch index) and run unchanged at 800 K.
Output: data/processed/revision/n4_selection.csv"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
m = pd.read_csv(ROOT / "data/processed/sto_production/mode_table_T.csv.gz", comment="#")
CLS = {"A": "acoustic", "soft": "low optic", "hard": "hard optic"}


def pick(table_mesh, job_mesh):
    t = m[(m.mesh == table_mesh) & (m.T_K == 300)].copy()
    q = t[["q_frac_x", "q_frac_y", "q_frac_z"]].to_numpy()
    t = t[np.all(np.abs(q * job_mesh - np.round(q * job_mesh)) < 1e-3, axis=1)]  # CSV stores 6 digits (1/6 = 0.166667)
    kP = m[(m.mesh == table_mesh) & (m.T_K == 300)].kappa_mode_W_mK.sum()
    t = t.sort_values(["kappa_mode_W_mK", "q_index", "branch"], ascending=[False, True, True])
    out, used = [], set()
    for c in CLS:
        n = 0
        for _, r in t[t.sector_char == c].iterrows():
            if n == 2:
                break
            if r.q_index in used:
                continue
            used.add(r.q_index); n += 1
            out.append({"sample": f"{job_mesh}^3", "table_mesh": table_mesh, "class": CLS[c], "q_index": int(r.q_index),
                        "q_frac_x": r.q_frac_x, "q_frac_y": r.q_frac_y, "q_frac_z": r.q_frac_z, "branch": int(r.branch),
                        "omega_cm1_300": r.omega_cm1, "gamma3_prod_300": r.gamma_hwhm_cm1, "eta3_300": r.eta,
                        "kappa_mode_300": r.kappa_mode_W_mK, "kappa_share_pct": 100 * r.kappa_mode_W_mK / kP,
                        "multiplicity": int(r.multiplicity)})
    return out


sel = pd.DataFrame(pick(16, 8) + pick(12, 6))
for T in (800,):  # same (q, branch) at 800 K: production values and character for reference
    t = m[m.T_K == T].set_index(["mesh", "q_index", "branch"])
    for col, src in (("omega_cm1_800", "omega_cm1"), ("gamma3_prod_800", "gamma_hwhm_cm1"), ("class_char_800", "sector_char")):
        sel[col] = [t.loc[(r.table_mesh, r.q_index, r.branch), src] for r in sel.itertuples()]
sel.to_csv(ROOT / "data/processed/revision/n4_selection.csv", index=False, float_format="%.6g")
pd.set_option("display.width", 250)
print(sel.to_string(index=False))


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"n4_selection.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['data/processed/sto_production/mode_table_T.csv.gz'])
