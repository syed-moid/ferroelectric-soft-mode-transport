#!/usr/bin/env python3
"""P2: kappa at 300 K of the steps a-c-d-e against the archived Fig. 7 value, and
the mode comparison of one run (e-kappa) with the per-mode file of the 2015 study (TAU_FILE: Omega, Gamma, tau).

Per run: kappa_P = xx entry of the .kl row at 300 K (and trace/3); kappa_C = xx entry of .kl_coherent at 300 K; the
mode sum of the group-averaged table (kp_eta_distribution.group_average_gamma) as a closure check.

Mode comparison (--modes=LABEL): positional matching of P0.6: point by order of the
run's irreducible 12^3 list, branch by rank inside the point, the 3 acoustic modes at Gamma omitted. Our side:
group-averaged HWHM. Degeneracy patterns of the run against the file are recomputed (P0.6 recorded points 25 and
75 as different against model A); differences over: all rows; all rows except the points whose pattern differs;
the carriers of kappa_P in this run (frozen 50 % rule, whole degenerate groups at the boundary, as t4_analysis.py);
the lowest optical triplet at Gamma and the lowest triplet at R explicitly.

Usage: p2_analysis.py ARCH_KL_DAT TAU_FILE OUT_JSON --modes=LABEL LABEL=RTA_DIR [LABEL=RTA_DIR ...]
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fekappa.modes import cell_volume_ang3, group_average_gamma, mesh_of  # noqa: E402
from p0_archive_files import blocks_from_tau, deg_pattern  # noqa: E402
from fekappa.modes import deg_groups  # noqa: E402

from fekappa.kappa_decomposition import mode_kappa_table  # noqa: E402


def kl_row(d, suffix):
    f = next(Path(d).glob(f"*.{suffix}"))
    for l in f.read_text().splitlines():
        s = l.split()
        if s and not l.startswith("#") and abs(float(s[0]) - 300) < 1e-6:
            v = [float(x) for x in s[1:]]
            diag = (v[0], v[4], v[8]) if len(v) == 9 else tuple(v[:3])     # .kl: 9 components; .kl_coherent: xx yy zz
            return {"xx": diag[0], "trace_over_3": sum(diag) / 3, "file": f.name}
    raise ValueError(f)


def table(d):
    df = mode_kappa_table(next(Path(d).glob("*.result")), nk=mesh_of(Path(d)) ** 3, V_cell_ang3=cell_volume_ang3(Path(d)))
    df = df[df.T_K == 300].copy()
    g = group_average_gamma(df)
    with np.errstate(divide="ignore", invalid="ignore"):
        df["kappa_mode_W_mK"] = np.where(df.kappa_mode_W_mK > 0, df.kappa_mode_W_mK * df.gamma_hwhm_cm1 / g, 0.0)
    df["gamma_avg_cm1"] = g
    df = df.sort_values(["q_index", "branch"]).reset_index(drop=True)
    df["grp"] = np.concatenate([deg_groups(h.omega_cm1.to_numpy()) for _, h in df.groupby("q_index")])
    return df


def stats(d_abs, d_rel, rows):
    if len(d_abs) == 0:
        return None
    k = int(np.argmax(np.abs(d_abs)))
    r = rows.iloc[k]
    return {"n": int(len(d_abs)), "median_abs": float(np.median(np.abs(d_abs))),
            "q95_abs": float(np.quantile(np.abs(d_abs), 0.95)), "max_abs": float(np.abs(d_abs[k])),
            "median_abs_rel": float(np.median(np.abs(d_rel))), "q95_abs_rel": float(np.quantile(np.abs(d_rel), 0.95)),
            "max_abs_rel": float(np.abs(d_rel).max()),
            "max_abs_at": {"q_index": int(r.q_index), "q": [float(r.q_frac_x), float(r.q_frac_y), float(r.q_frac_z)],
                           "branch": int(r.branch), "omega_ours": float(r.omega_cm1), "omega_file": float(r.w_file),
                           "gamma_ours": float(r.gamma_avg_cm1), "gamma_file": float(r.g_file)}}


def main():
    opts = dict(a[2:].split("=", 1) for a in sys.argv[1:] if a.startswith("--"))
    argv = [a for a in sys.argv[1:] if not a.startswith("--")]
    arch_kl, tau_p, out = argv[:3]
    runs = dict(a.split("=", 1) for a in argv[3:])
    arch = {float(l.split()[0]): float(l.split()[1]) for l in open(arch_kl) if l.strip()}
    A = arch[300.0]
    R = {"archived_kappa_300K": A, "runs": {}}
    tabs = {}
    for lab, d in runs.items():
        kp, kc = kl_row(d, "kl"), kl_row(d, "kl_coherent")
        df = tabs[lab] = table(d)
        R["runs"][lab] = {"dir": d, "kappa_P_xx": kp["xx"], "kappa_P_trace3": kp["trace_over_3"],
                          "kappa_C_xx": kc["xx"], "kappa_C_trace3": kc["trace_over_3"],
                          "kappa_P_modesum_groupavg": float(df.kappa_mode_W_mK.sum()),
                          "minus_archived_abs": kp["xx"] - A, "minus_archived_rel": (kp["xx"] - A) / A,
                          "n_irred_q": int(df.q_index.nunique())}
    lab = opts.get("modes")
    if lab:
        df = tabs[lab]
        tau = np.loadtxt(tau_p)
        nq = df.q_index.nunique()
        blocks, used = blocks_from_tau(tau, nq)
        assert used == len(tau), (used, len(tau))
        q0 = int(df.q_index.min())
        parts, patdiff = [], []
        for (qi, g), b in zip(df.groupby("q_index"), blocks):
            g = g.iloc[3:] if int(qi) == q0 else g
            g = g.copy()
            g["w_file"], g["g_file"] = b[:, 0], b[:, 1]
            po, pf = deg_pattern(g.omega_cm1.to_numpy()), deg_pattern(b[:, 0])
            if po != pf:
                patdiff.append({"q_index": int(qi), "q": [float(x) for x in g.iloc[0][["q_frac_x", "q_frac_y", "q_frac_z"]]],
                                "pattern_ours": po, "pattern_file": pf,
                                "rows": g[["branch", "omega_cm1", "w_file", "gamma_avg_cm1", "g_file"]].round(6).values.tolist()})
            parts.append(g)
        M = pd.concat(parts).reset_index(drop=True)
        M["dW"], M["dG"] = M.omega_cm1 - M.w_file, M.gamma_avg_cm1 - M.g_file
        M["dW_rel"], M["dG_rel"] = M.dW / M.w_file, M.dG / M.g_file
        kP = df.kappa_mode_W_mK.sum()
        gk = df.groupby(["q_index", "grp"]).kappa_mode_W_mK.sum().sort_values(ascending=False)
        cum = gk.cumsum() / kP
        sel = set(cum.index[: int(np.searchsorted(cum.to_numpy(), 0.5)) + 1])
        M["carrier"] = [(q, gr) in sel for q, gr in zip(M.q_index, M.grp)]
        bad = {p["q_index"] for p in patdiff}
        subsets = {"all": M, "excluding_pattern_differences": M[~M.q_index.isin(bad)], "carriers_50pct": M[M.carrier]}
        C = {}
        for name, S in subsets.items():
            S = S.reset_index(drop=True)
            C[name] = {"omega": stats(S.dW.to_numpy(), S.dW_rel.to_numpy(), S),
                       "gamma": stats(S.dG.to_numpy(), S.dG_rel.to_numpy(), S)}
        qR = M.loc[(M.q_frac_x.abs() - 0.5).abs().lt(1e-6) & (M.q_frac_y.abs() - 0.5).abs().lt(1e-6)
                   & (M.q_frac_z.abs() - 0.5).abs().lt(1e-6), "q_index"].unique()
        ex = {"Gamma_lowest_optical_triplet": M[M.q_index == q0].head(3), "R_lowest_triplet": M[M.q_index == qR[0]].head(3)}
        R["modes"] = {"run": lab, "rows_matched": int(len(M)), "rows_file": int(len(tau)),
                      "pattern_differences": patdiff, "n_pattern_differences": len(patdiff),
                      "carriers": {"groups": len(sel), "modes_irred": int(M.carrier.sum()),
                                   "kappa_share": float(cum.iloc[len(sel) - 1])},
                      "comparison": C,
                      "explicit": {k: v[["q_index", "branch", "omega_cm1", "w_file", "gamma_avg_cm1", "g_file"]]
                                   .round(6).values.tolist() for k, v in ex.items()}}
    json.dump(R, open(out, "w"), indent=1)
    for k, v in R["runs"].items():
        print(f"{k:6s} kP_xx {v['kappa_P_xx']:.6f} (tr/3 {v['kappa_P_trace3']:.6f}, modesum {v['kappa_P_modesum_groupavg']:.6f}) "
              f"kC_xx {v['kappa_C_xx']:.6f}  - archived {v['minus_archived_abs']:+.6f} ({100 * v['minus_archived_rel']:+.3f} %)")
    if lab:
        m = R["modes"]
        print("modes", lab, "matched", m["rows_matched"], "/", m["rows_file"], "pattern diffs",
              [(p["q_index"], p["pattern_ours"], p["pattern_file"]) for p in m["pattern_differences"]])
        for name, c in m["comparison"].items():
            print(f" {name:30s} dW med {c['omega']['median_abs']:.3f} q95 {c['omega']['q95_abs']:.3f} max {c['omega']['max_abs']:.3f}"
                  f" | dG med {c['gamma']['median_abs']:.4f} q95 {c['gamma']['q95_abs']:.4f} max {c['gamma']['max_abs']:.4f}"
                  f" rel med {c['gamma']['median_abs_rel']:.3f}")
        print(" carriers", m["carriers"])
        for k, v in m["explicit"].items():
            print(" ", k, v)


if __name__ == "__main__":
    main()
