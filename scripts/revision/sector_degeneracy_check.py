"""Sector rule audit: ties and exactly degenerate subspaces under the character rule.

The per-eigenvector rule (production until 22-Oct; column sector_char_eig) labels each e_j(q) by argmax_S c_S(q, j), with
c_S(q, j) = sum_{i in S(Gamma)} |<e_i(Gamma)|e_j(q)>|^2 and S = A (Gamma branches 1-3), soft (4-6), hard (7-15), the Gamma
eigenvectors taken from the same (T, mesh) run. Inside an exactly degenerate subspace at q the individual e_j are fixed only up
to a unitary rotation (LAPACK's choice), so a per-eigenvector label is basis-dependent there. Basis-invariant variant: sum c_S
over the subspace and give every member the label of the largest sum; this is the production rule from 22-Oct
(prod_analysis.py, column sector_char), and the script asserts that it reproduces it. Columns *_char_* below refer to the
per-eigenvector rule, *_inv_* to the subspace rule.

Per (T, mesh): number of degenerate subspaces, number whose members carry different per-eigenvector labels, the kappa_P
share carried by those, ties (top two c_S within 1e-6), and the hard/hard share of kappa_C under both variants.
Output: data/processed/revision/sector_degeneracy.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd

from fekappa.kappa_decomposition import read_kc_elem

ROOT = Path(__file__).resolve().parents[2]
ARCH = ROOT / "data/raw/alamode_sto/revision/pm_q12"
OUT = ROOT / "data/processed/revision/sector_degeneracy.csv"
TOL = 1e-3   # cm^-1; frequencies in the mode table carry 3-4 decimals
LAB = np.array(["A", "soft", "hard"])

modes = pd.read_csv(ROOT / "data/processed/sto_production/mode_table_T.csv.gz", comment="#")
rows = []
for (T, mesh), df in modes.groupby(["T_K", "mesh"]):
    df = df.sort_values(["q_index", "omega_cm1"]).copy()
    c = df[["c_A", "c_soft", "c_hard"]].to_numpy()
    s = np.sort(c, axis=1)
    ties = int(((s[:, 2] - s[:, 1]) < 1e-6).sum())
    grp = (df.groupby("q_index").omega_cm1.transform(lambda w: np.concatenate([[0], np.cumsum(np.diff(w.to_numpy()) > TOL)])))
    df["sub"] = df.q_index.astype(str) + ":" + grp.astype(str)
    size = df.groupby("sub").branch.transform("size")
    csum = df.groupby("sub")[["c_A", "c_soft", "c_hard"]].transform("sum").to_numpy()
    df["sector_inv"] = LAB[np.argmax(csum, axis=1)]
    nlab = df.groupby("sub").sector_char_eig.transform("nunique")
    deg = size > 1
    mixed = deg & (nlab > 1)
    kP = df.kappa_mode_W_mK.sum()
    out = {"T_K": T, "mesh": mesh, "n_subspaces_deg": int(df[deg]["sub"].nunique()),
           "n_subspaces_mixed": int(df[mixed]["sub"].nunique()), "kappaP_share_mixed_pct": 100 * df[mixed].kappa_mode_W_mK.sum() / kP,
           "n_ties": ties, "n_modes_relabelled": int((df.sector_inv != df.sector_char_eig).sum())}
    assert (df.sector_inv == df.sector_char).all(), (T, mesh)
    job = f"rta_pm_q12_T{T}_m{mesh}"
    f = ARCH / job / f"{job}.kc_elem"
    if T >= 300 and f.exists():
        kc = read_kc_elem(f)
        kc = kc[(kc.i == kc.j) & (kc.ibranch != kc.jbranch)]
        for var, col in (("char", "sector_char_eig"), ("inv", "sector_inv")):
            lab = df.set_index(["q_index", "branch"])[col]
            s1 = lab.loc[list(zip(kc.ik, kc.ibranch))].to_numpy()
            s2 = lab.loc[list(zip(kc.ik, kc.jbranch))].to_numpy()
            hh = kc.kc_real.to_numpy()[(s1 == "hard") & (s2 == "hard")].sum() / kc.kc_real.sum()
            out[f"hh_share_{var}_pct"] = 100 * hh
        lab = df.set_index(["q_index", "branch"])
        for var, col in (("char", "sector_char_eig"), ("inv", "sector_inv")):
            sh = df.groupby(col).kappa_mode_W_mK.sum() / kP
            out[f"kP_hard_{var}_pct"] = 100 * sh.get("hard", 0.0)
    rows.append(out)
res = pd.DataFrame(rows)
with open(OUT, "w") as fh:
    fh.write("# character-rule audit: degenerate subspaces, ties, basis-invariant variant (scripts/sector_degeneracy_check.py)\n")
    res.to_csv(fh, index=False, float_format="%.4f")
pd.set_option("display.width", 220)
print(res.round(3).to_string(index=False))


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"sector_degeneracy.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['data/processed/sto_production/mode_table_T.csv.gz', 'zenodo:pm_q12_evec.tar.gz', 'zenodo:pm_q12_rta.tar.gz (*.result)'])
