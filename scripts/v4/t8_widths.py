#!/usr/bin/env python3
"""Full-sum comparison: linewidths and κ_P of two RTA runs of the same
model and mesh, TRISYM = 1 (production) against TRISYM = 0. Per-branch `.result` widths and degenerate-group-averaged
widths (what κ uses): relative difference (TRISYM 0 − TRISYM 1) / TRISYM 1 over all modes with a non-zero width;
median, 95 % quantile and maximum of the absolute value, signed mean, location of the maximum; the R triplet (lowest modes at q = (½, ½, ½)); κ_P
(mode sum and `.kl` xx). Same measures as the T3 report tables (Markdown) (8³). Tables only.

Usage: t8_widths.py RTA_DIR_TRISYM1 RTA_DIR_TRISYM0 OUT_JSON  (prints a Markdown table)
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fekappa.modes import cell_volume_ang3, group_average_gamma, mesh_of  # noqa: E402
from t4_analysis import T_K  # noqa: E402

from fekappa.kappa_decomposition import mode_kappa_table  # noqa: E402


def load(d):
    d = Path(d)
    res = next(d.glob("*.result"))
    df = mode_kappa_table(res, nk=mesh_of(d) ** 3, V_cell_ang3=cell_volume_ang3(d))
    df = df[df.T_K == T_K].sort_values(["q_index", "branch"]).reset_index(drop=True)
    g = group_average_gamma(df)
    with np.errstate(divide="ignore", invalid="ignore"):          # κ with group-averaged widths, as in T4
        df["kappa_mode_W_mK"] = np.where(df.kappa_mode_W_mK > 0, df.kappa_mode_W_mK * df.gamma_hwhm_cm1 / g, 0.0)
    df["gamma_avg_cm1"] = g
    kl = [ln.split() for ln in open(next(d.glob("*.kl"))) if ln.strip() and not ln.startswith("#")]
    kxx = [float(r[1]) for r in kl if abs(float(r[0]) - T_K) < 1e-6][0]
    return df, kxx


def stats(a, b, df):
    m = (a > 0) & (b > 0)
    r = (b[m] - a[m]) / a[m]
    i = np.argmax(np.abs(r))
    loc = df[m].iloc[i]
    return {"n": int(m.sum()), "median_abs": float(np.median(np.abs(r))), "q95_abs": float(np.quantile(np.abs(r), 0.95)),
            "max_abs": float(np.abs(r[i])), "signed_mean": float(r.mean()),
            "max_at": {"q_index": int(loc.q_index), "branch": int(loc.branch), "omega_cm1": float(loc.omega_cm1)}}


def main():
    (d1, k1), (d0, k0) = load(sys.argv[1]), load(sys.argv[2])
    assert (d1.q_index.to_numpy() == d0.q_index.to_numpy()).all() and (d1.branch.to_numpy() == d0.branch.to_numpy()).all()
    assert np.allclose(d1.omega_cm1, d0.omega_cm1, rtol=0, atol=1e-9), "frequencies differ"
    out = {"T_K": T_K,
           "per_branch": stats(d1.gamma_hwhm_cm1.to_numpy(), d0.gamma_hwhm_cm1.to_numpy(), d1),
           "group_averaged": stats(d1.gamma_avg_cm1.to_numpy(), d0.gamma_avg_cm1.to_numpy(), d1),
           "kappa_P_modesum": [float(d1.kappa_mode_W_mK.sum()), float(d0.kappa_mode_W_mK.sum())],
           "kappa_kl_xx": [k1, k0]}
    atR = np.isclose(d1[["q_frac_x", "q_frac_y", "q_frac_z"]].abs().to_numpy(), 0.5).all(axis=1)
    lo = d1.omega_cm1[atR].min()                                   # R triplet: lowest modes at q = (½, ½, ½)
    sel = atR & np.isclose(d1.omega_cm1, lo, atol=1e-3)
    out["lowest_triplet"] = {"omega_cm1": float(lo), "n": int(sel.sum()),
                             "gamma_avg_cm1": [float(d1.gamma_avg_cm1[sel].mean()), float(d0.gamma_avg_cm1[sel].mean())]}
    json.dump(out, open(sys.argv[3], "w"), indent=1)
    g, p = out["group_averaged"], out["per_branch"]
    print("| quantity | TRISYM 1 | TRISYM 0 | relative difference |\n|---|---|---|---|")
    print(f"| κ_P mode sum (W/mK) | {out['kappa_P_modesum'][0]:.4f} | {out['kappa_P_modesum'][1]:.4f} | "
          f"{100 * (out['kappa_P_modesum'][1] / out['kappa_P_modesum'][0] - 1):+.3f} % |")
    print(f"| κ_P `.kl` xx (W/mK) | {k1:.4f} | {k0:.4f} | {100 * (k0 / k1 - 1):+.3f} % |")
    t = out["lowest_triplet"]
    print(f"| R triplet ({t['omega_cm1']:.2f} cm⁻¹, {t['n']} modes): group-averaged Γ (cm⁻¹) | "
          f"{t['gamma_avg_cm1'][0]:.4f} | {t['gamma_avg_cm1'][1]:.4f} | {100 * (t['gamma_avg_cm1'][1] / t['gamma_avg_cm1'][0] - 1):+.2f} % |")
    print("\n| widths | modes | median \\|Δ\\| | 95 % quantile | max | signed mean | max at (q index, branch, Ω cm⁻¹) |")
    print("|---|---|---|---|---|---|---|")
    for nm, s in (("group-averaged (what κ uses)", g), ("per branch (`.result`)", p)):
        a = s["max_at"]
        print(f"| {nm} | {s['n']} | {100 * s['median_abs']:.2f} % | {100 * s['q95_abs']:.2f} % | {100 * s['max_abs']:.1f} % | "
              f"{s['signed_mean']:+.2e} | {a['q_index']}, {a['branch']}, {a['omega_cm1']:.2f} |")


if __name__ == "__main__":
    main()
