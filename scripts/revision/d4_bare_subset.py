"""D4: bare-harmonic diagnostic — three framings, with a CHARACTER-matched common subset.

U = the 11 irreducible bare-unstable modes (36 of 7680 mesh modes, 8^3), T-independent.
For each SCPH mode j at an irreducible q hosting U modes:
  w_j = sum_{b in U_q} |<e_b(bare) | e_j(SCPH)>|^2   (weight of the bare-unstable subspace in mode j)
Character-matched subsets of the SCPH kappa_P:
  hard : drop, at each q, the |U_q| SCPH modes with the largest w_j
  frac : sum_j (1 - w_j) kappa_j
(The bare side is unchanged: anphon assigns U zero kappa and removes U as scattering partners.)
Index subset (q, branch-order) and |omega| substitution from c2_*.csv for comparison.
SCPH variants: production (own SCPH, PM, q1 12^3) at 250-1000 K; superseded q1 8^3 at 300 K;
submitted (tutorial dfc2) at 50-1000 K.
Output: data/processed/revision/d4_bare_framings.csv, d4_character_match_300K.csv
"""
from fekappa.provenance import external_dir
import sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
from fekappa.evec_io import irred_evecs
from fekappa.kappa_decomposition import mode_kappa_table  # noqa: E402

ROOT = HERE.parents[1]
RAW = ROOT / "data/raw/alamode_sto"
RUNS = external_dir("FEKAPPA_RUNS", "the pre-queue run directory (Zenodo prequeue_runs)")
OUT = ROOT / "data/processed/revision"

bare_res = RAW / "STO_RTA_production.result"
bare = mode_kappa_table(bare_res)
U = bare[(bare.T_K == bare.T_K.min()) & (bare.omega_cm1 < 0)][["q_index", "branch"]]
Uq = U.groupby("q_index")["branch"].apply(list).to_dict()
Eb, _ = irred_evecs(RUNS / "evec_d4/ev_bare.mesh.evec", bare_res)


def framings(scph_result, scph_evec, T, label):
    s = mode_kappa_table(scph_result)
    s = s[s.T_K == T].copy()
    Es, _ = irred_evecs(scph_evec, scph_result)
    s["w"] = 0.0
    match = []
    for q, bl in Uq.items():
        P = sum(np.abs(np.einsum("c,jc->j", Eb[q - 1, b - 1].conj(), Es[q - 1])) ** 2 for b in bl)
        for j in range(15):
            s.loc[(s.q_index == q) & (s.branch == j + 1), "w"] = P[j]
        top = np.argsort(-P)[: len(bl)]
        for j in top:
            r = s[(s.q_index == q) & (s.branch == j + 1)].iloc[0]
            match.append({"variant": label, "T_K": T, "q_index": q, "bare_branches": bl, "scph_branch": j + 1,
                          "scph_omega_cm1": r.omega_cm1, "w": P[j], "kappa_mode": r.kappa_mode_W_mK})
        s.loc[(s.q_index == q) & (s.branch.isin([j + 1 for j in top])), "hard_drop"] = True
    s["hard_drop"] = s["hard_drop"].fillna(False).astype(bool) if "hard_drop" in s else False
    idx = s.apply(lambda r: (r.q_index in Uq) and (r.branch in Uq[r.q_index]), axis=1)
    kb = bare[bare.T_K == T].kappa_mode_W_mK.sum()
    full = s.kappa_mode_W_mK.sum()
    row = {"variant": label, "T_K": T, "bare_kappaP": kb, "scph_full": full,
           "scph_subset_index": s.loc[~idx, "kappa_mode_W_mK"].sum(),
           "scph_subset_char_hard": s.loc[~s.hard_drop, "kappa_mode_W_mK"].sum(),
           "scph_subset_char_frac": ((1 - s.w) * s.kappa_mode_W_mK).sum()}
    for k in ("scph_full", "scph_subset_index", "scph_subset_char_hard", "scph_subset_char_frac"):
        row["ratio_" + k.replace("scph_", "")] = row[k] / kb
    return row, match


rows, matches = [], []
r, mt = framings(RUNS / "rta_on_k8_m8/rta_on_k8_m8.result", RUNS / "evec_d4/ev_pmq8_300K.mesh.evec", 300, "superseded_pm_q8")
rows.append(r); matches += mt
PROD = ROOT / "data/raw/alamode_sto/revision/pm_q12"
for T in range(250, 1001, 50):
    r, mt = framings(PROD / f"rta_pm_q12_T{T}_m8/rta_pm_q12_T{T}_m8.result",
                     PROD / f"evec_pm_q12_T{T}_m8/evec_pm_q12_T{T}_m8.mesh.evec", T, "production_pm_q12")
    rows.append(r)
    if T == 300:
        matches += mt
for T in range(50, 1001, 50):
    r, mt = framings(RAW / f"STO_RTA_scph_{T}K.result", RUNS / f"evec8/ev_{T}K.mesh.evec", T, "submitted_q2_nopm")
    rows.append(r)
    if T == 300:
        matches += mt
df = pd.DataFrame(rows)
absw = pd.read_csv(OUT / "c2_absw_bare.csv", comment="#")[["T_K", "bare_absw"]]
df = df.merge(absw, on="T_K", how="left")
df["ratio_absw"] = df.scph_full / df.bare_absw
df.to_csv(OUT / "d4_bare_framings.csv", index=False, float_format="%.5g")
pd.DataFrame(matches).to_csv(OUT / "d4_character_match_300K.csv", index=False, float_format="%.4g")
pd.set_option("display.width", 250)
print(pd.DataFrame(matches).round(3).to_string(index=False))
print(df.round(3).to_string(index=False))


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"d4_bare_framings.csv", f"d4_character_match_300K.csv", f"c2_absw_bare.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['data/raw/alamode_sto (distributed example outputs; IMPORT_NOTE.md)', 'zenodo:pm_q12_rta.tar.gz (rta_pm_q12_T*_m8)', 'zenodo:pm_q12_evec.tar.gz (evec_pm_q12_T*_m8)', '$FEKAPPA_RUNS/{evec8,evec_d4,rta_on_k8_m8,bare_stock,bare_absw} (zenodo:prequeue_runs.tar.gz)'])
