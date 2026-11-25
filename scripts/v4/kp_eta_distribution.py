#!/usr/bin/env python3
"""Distribution of kappa_P over eta = HWHM / Omega at one temperature.

Writes two outputs:
- the shares of kappa_P carried by modes with eta above 0.02, 0.05 and 0.1 (one row per run);
- the cumulative curve C(eta) = share of kappa_P carried by modes with eta' <= eta, for every run, on a log grid
  of eta and at each mode's own eta.

Same mode table as a3_transport_analysis.py: fekappa.kappa_decomposition.mode_kappa_table, modes with
Omega > 1e-3 cm^-1.

Option --gamma=avg: widths averaged over degenerate groups exactly as anphon does before forming kappa
(Conductivity::average_self_energy_at_degenerate_point: modes at the same irreducible q, in branch order, join the
group while |omega - omega_first| < 1e-7 Ry); the per-mode kappa is rescaled by gamma_raw / gamma_avg (tau = 1/2Gamma).
Default --gamma=raw: the per-branch widths printed in .result (previous behaviour, unchanged).

Usage: kp_eta_distribution.py [--gamma=raw|avg] T OUT_PREFIX LABEL=RTA_DIR [LABEL=RTA_DIR ...]
Writes OUT_PREFIX_shares.csv and OUT_PREFIX_cumulative.csv.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from fekappa.kappa_decomposition import mode_kappa_table

THRESH = (0.02, 0.05, 0.1)
from fekappa.modes import TOL_DEG_CM1, cell_volume_ang3, group_average_gamma, mesh_of  # noqa: E402,F401


def main():
    opts = [a for a in sys.argv[1:] if a.startswith("--gamma=")]
    argv = [a for a in sys.argv[1:] if not a.startswith("--gamma=")]
    gmode = opts[-1].split("=", 1)[1] if opts else "raw"
    assert gmode in ("raw", "avg"), gmode
    T, prefix = float(argv[0]), argv[1]
    shares, curves = [], []
    grid = np.logspace(-4, 0, 201)
    for spec in argv[2:]:
        lab, d = spec.split("=", 1)
        d = Path(d)
        mesh = mesh_of(d)
        df = mode_kappa_table(next(d.glob("*.result")), nk=mesh ** 3, V_cell_ang3=cell_volume_ang3(d))
        df = df[df.T_K == T].copy()
        if gmode == "avg":
            g = group_average_gamma(df)
            with np.errstate(divide="ignore", invalid="ignore"):
                df["kappa_mode_W_mK"] = np.where(df.kappa_mode_W_mK > 0, df.kappa_mode_W_mK * df.gamma_hwhm_cm1 / g, 0.0)
            df["gamma_hwhm_cm1"] = g
        df = df[df.omega_cm1 > 1e-3].copy()
        df["eta"] = df.gamma_hwhm_cm1 / df.omega_cm1
        kP = df.kappa_mode_W_mK.sum()
        top = df.loc[df.eta.idxmax()]
        row = {"run": lab, "T_K": T, "mesh": mesh, "gamma": gmode, "kappa_P_modesum": kP, "eta_max": df.eta.max(),
               "eta_max_q_index": int(top.q_index), "eta_max_branch": int(top.branch),
               "eta_max_omega_cm1": float(top.omega_cm1)}
        for t in THRESH:
            row[f"share_kP_eta_gt_{t}"] = df[df.eta > t].kappa_mode_W_mK.sum() / kP
        shares.append(row)
        s = df.sort_values("eta")
        cum = s.kappa_mode_W_mK.cumsum().to_numpy() / kP
        cg = np.interp(grid, s.eta.to_numpy(), cum, left=0.0, right=1.0)
        curves.append(pd.DataFrame({"run": lab, "kind": "grid", "eta": grid, "cum_share_kP_le_eta": cg}))
        curves.append(pd.DataFrame({"run": lab, "kind": "mode", "eta": s.eta.to_numpy(), "cum_share_kP_le_eta": cum}))
        print({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in row.items()})
    pd.DataFrame(shares).to_csv(f"{prefix}_shares.csv", index=False, float_format="%.6g")
    pd.concat(curves).to_csv(f"{prefix}_cumulative.csv", index=False, float_format="%.6g")


if __name__ == "__main__":
    main()
