#!/usr/bin/env python3
"""P6: where the conductivity difference at 300 K sits. Descriptive.

Mode sum of run e-kappa (anphon's form, as fekappa.kappa_decomposition.mode_kappa_table; group-averaged widths) and
the same sum with the lifetimes (step 1) and also the heat capacities (step 2) taken from the 2015 per-mode file (TAU_FILE), paired by
position as in scripts/p2_analysis.py (P2 §3). Mismatched group (point 25, branches 10-11): our values; also given
with those modes, and with all of point 25, left out of both sums. Splits of the differences by band of our
frequency (< 105, 105-145, 145-160, >= 160 cm^-1) and by acoustic / optical branch index. Repeat of step 1 with the
lifetimes paired by sorted frequency inside each wave-vector group. Output private (per-mode file values).

Usage: p6_kappa_split.py RTA_DIR TAU_FILE KL_ARCHIVED OUT_JSON
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from p0_archive_files import blocks_from_tau, deg_pattern  # noqa: E402
from p2_analysis import kl_row, table  # noqa: E402

from fekappa.alamode_io import tau_ps_from_gamma_hwhm_cm1  # noqa: E402
from fekappa.kappa_decomposition import CM1_TO_MEV  # noqa: E402
from fekappa.kappa_kernel import heat_capacity_mode  # noqa: E402

BANDS = [("<105", 0, 105), ("105-145", 105, 145), ("145-160", 145, 160), (">=160", 160, np.inf)]
MISMATCH = (25, (10, 11))
T_K = 300.0


def splits(M, d):
    out = {"bands": {}, "branch": {}}
    for nm, lo, hi in BANDS:
        s = (M.omega_cm1 >= lo) & (M.omega_cm1 < hi)
        out["bands"][nm] = {"n": int(s.sum()), "delta": float(d[s].sum())}
    ac = M.branch <= 3
    out["branch"] = {"acoustic_1_3": {"n": int(ac.sum()), "delta": float(d[ac].sum())},
                     "optical_4_15": {"n": int((~ac).sum()), "delta": float(d[~ac].sum())}}
    out["sum_of_parts_bands"] = sum(v["delta"] for v in out["bands"].values())
    out["total"] = float(d.sum())
    return out


def main():
    rdir, tau_p, kl_arch, out = sys.argv[1:5]
    arch = {float(l.split()[0]): float(l.split()[1]) for l in open(kl_arch) if l.strip()}[T_K]
    kl = kl_row(rdir, "kl")["xx"]
    df = table(rdir)                                        # group-averaged widths; kappa_mode rescaled to them
    tau = np.loadtxt(tau_p)
    nq = df.q_index.nunique()
    blocks, used = blocks_from_tau(tau, nq)
    assert used == len(tau) == 1257, (used, len(tau))
    q0 = int(df.q_index.min())
    parts = []
    for (qi, g), b in zip(df.groupby("q_index"), blocks):
        g = (g.iloc[3:] if int(qi) == q0 else g).copy()
        assert len(g) == len(b)
        g["w_file"], g["tau_file_ps"] = b[:, 0], b[:, 2]
        # sorted-frequency pairing inside the group (stable sorts)
        oo = np.argsort(g.omega_cm1.to_numpy(), kind="stable")
        of = np.argsort(b[:, 0], kind="stable")
        ts = np.empty(len(b))
        ts[oo] = b[of, 2]
        g["tau_file_sorted_ps"] = ts
        g["pattern_differs"] = deg_pattern(g.omega_cm1.to_numpy()) != deg_pattern(b[:, 0])
        parts.append(g)
    import pandas as pd
    M = pd.concat(parts).reset_index(drop=True)
    excluded_gamma = df[(df.q_index == q0)].head(3)
    assert float(excluded_gamma.kappa_mode_W_mK.abs().sum()) == 0.0
    M["tau_ours_ps"] = tau_ps_from_gamma_hwhm_cm1(M.gamma_avg_cm1.to_numpy())
    M["C_ours"] = heat_capacity_mode(M.omega_cm1.to_numpy() * CM1_TO_MEV, np.full(len(M), T_K))
    M["C_file"] = heat_capacity_mode(M.w_file.to_numpy() * CM1_TO_MEV, np.full(len(M), T_K))
    mm = (M.q_index == MISMATCH[0]) & M.branch.isin(MISMATCH[1])
    p25 = M.q_index == MISMATCH[0]
    assert set(M.loc[M.pattern_differs, "q_index"]) == {MISMATCH[0]}, set(M.loc[M.pattern_differs, "q_index"])

    k0 = M.kappa_mode_W_mK.to_numpy()
    r_tau = np.where(mm, 1.0, M.tau_file_ps / M.tau_ours_ps)
    r_tau_sorted = np.where(mm, 1.0, M.tau_file_sorted_ps / M.tau_ours_ps)
    r_C = np.where(mm, 1.0, M.C_file / M.C_ours)
    k1 = k0 * r_tau
    k2 = k0 * r_tau * r_C
    k1s = k0 * r_tau_sorted
    d1, d2, d1s = k1 - k0, k2 - k0, k1s - k0

    def summ(v, keep):
        return float(v[keep].sum())

    allm = np.ones(len(M), bool)
    R = {"kl_xx": kl, "archived": arch, "n_modes": int(len(M)),
         "step0": summ(k0, allm), "step0_minus_kl": summ(k0, allm) - kl}
    if abs(R["step0_minus_kl"]) > 1e-4:
        R["STOP"] = "step 0 does not reproduce .kl within 1e-4"
        json.dump(R, open(out, "w"), indent=1)
        print(json.dumps(R, indent=1))
        sys.exit(1)
    for lab, keep in (("all", allm), ("without_p25_b10_11", ~mm.to_numpy()), ("without_point25", ~p25.to_numpy())):
        R[lab] = {"n": int(keep.sum()), "step0": summ(k0, keep), "step1": summ(k1, keep), "step2": summ(k2, keep),
                  "step1_sorted": summ(k1s, keep)}
    R["split_step1"] = splits(M, d1)
    R["split_step2"] = splits(M, d2)
    R["split_step1_sorted"] = splits(M, d1s)
    R["split_step1_sorted_minus_positional"] = splits(M, k1s - k1)
    R["sorted_pairing_modes_changed"] = int((~np.isclose(M.tau_file_sorted_ps, M.tau_file_ps, rtol=0, atol=0)).sum())
    R["sorted_pairing_points_changed"] = sorted(int(q) for q in
                                                M.loc[M.tau_file_sorted_ps != M.tau_file_ps, "q_index"].unique())
    R["mismatch_modes"] = M.loc[mm, ["q_index", "branch", "omega_cm1", "w_file", "tau_ours_ps", "tau_file_ps",
                                     "kappa_mode_W_mK"]].round(8).values.tolist()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(R, open(out, "w"), indent=1)
    M[["q_index", "branch", "omega_cm1", "w_file", "gamma_avg_cm1", "tau_ours_ps", "tau_file_ps", "tau_file_sorted_ps",
       "C_ours", "C_file", "kappa_mode_W_mK"]].assign(k1=k1, k2=k2, k1_sorted=k1s).to_csv(
        Path(out).with_name("p6_modes.csv"), index=False)
    print(json.dumps({k: v for k, v in R.items() if not k.startswith("split")}, indent=1))


if __name__ == "__main__":
    main()
