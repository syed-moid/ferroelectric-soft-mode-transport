#!/usr/bin/env python3
"""Write the family-result tables of the T4 report tables (Markdown) from the t4_analysis summaries,
and the data of one figure per model. Prints nothing but the paths written.

Usage: t4_report.py SUMMARY_A SUMMARY_B MODES_A MODES_B SENS_A SENS_B CHECKS_JSON OUT_MD FIGDIR
"""
import json
import sys
from pathlib import Path

import pandas as pd

FAMS = ["neither", "pol", "rot", "mixed"]
DEFS = ["primary_c0.5", "secondary_c0.5", "primary_c0.4", "primary_c0.6"]


def pct(x):
    return f"{100 * x + 0.0:.2f}".replace("-0.00", "0.00")


def main():
    sA, sB = (json.load(open(p)) for p in sys.argv[1:3])
    mA, mB = (pd.read_csv(p) for p in sys.argv[3:5])
    senA, senB = (json.load(open(p)) for p in sys.argv[5:7])
    chk = json.load(open(sys.argv[7]))
    out, figdir = Path(sys.argv[8]), Path(sys.argv[9])
    figdir.mkdir(parents=True, exist_ok=True)
    S = {"A": sA, "B": sB}
    L = []
    w = L.append

    w("## 2. Results (opened here first)\n")
    w("All shares are κ_P-weighted means over the carriers (each carrier mode weighted by its κ_P contribution), of the")
    w("carrier's group-averaged three-phonon linewidth (decay + absorption). Soft = pol, rot or mixed. 300 K, 16³, q₁ 12³,")
    w("TRISYM = 1 (production widths).\n")
    w("### 2a. Carriers (frozen 50 % rule)\n")
    w("| | model A | model B |\n|---|---|---|")
    w(f"| κ_P mode sum (W/mK) | {sA['kappa_P_modesum']:.4f} | {sB['kappa_P_modesum']:.4f} |")
    w(f"| carrier modes (irreducible) / degenerate groups | {sA['carriers']} / {sA['carrier_groups']} | "
      f"{sB['carriers']} / {sB['carrier_groups']} |")
    w(f"| share of κ_P carried by the carriers | {pct(sA['carrier_kappa_share'])} % | {pct(sB['carrier_kappa_share'])} % |")
    for d in DEFS:
        cA, cB = (S[m]["definitions"][d]["carrier_family_counts"] for m in "AB")
        w(f"| carriers by family, {d} (neither / pol / rot / mixed) | "
          + " / ".join(str(cA[f]) for f in FAMS) + " | " + " / ".join(str(cB[f]) for f in FAMS) + " |")
    w("")
    w("### 2b. Linewidth of the carriers by number of soft partners\n")
    w("| definition | model | neither soft (%) | one soft (%) | both soft (%) | absorption part (%) |")
    w("|---|---|---|---|---|---|")
    for d in DEFS:
        for m in "AB":
            x = S[m]["definitions"][d]
            n = x["kappa_weighted_nsoft_frac"]
            w(f"| {d} | {m} | {pct(n['0'])} | {pct(n['1'])} | {pct(n['2'])} | {pct(x['kappa_weighted_absorption_share'])} |")
    w("\nMedians over the carriers (primary c0.5), neither / one / both soft: "
      + "; ".join(f"{m} " + " / ".join(pct(S[m]['definitions']['primary_c0.5']['median_nsoft_frac'][k])
                                       for k in "012") + " %" for m in "AB") + ".\n")
    w("### 2c. Linewidth of the carriers by partner family pair (primary c0.5; %)\n")
    pairs = list(sA["definitions"]["primary_c0.5"]["kappa_weighted_pair_frac"])
    w("| partner pair | model A | model B |\n|---|---|---|")
    for p in pairs:
        a, b = (S[m]["definitions"]["primary_c0.5"]["kappa_weighted_pair_frac"][p] for m in "AB")
        if a or b:
            w(f"| {p} | {pct(a)} | {pct(b)} |")
    w("\nPairs absent from the table are zero in both models. pol/pol: A −1.6 × 10⁻⁴ %, B 3.1 × 10⁻⁴ % (shown as 0.00).")
    w("Negative bins: 3 per model, all pol/pol absorption parts of carriers; largest magnitude 8.9 × 10⁻⁵ (A) and")
    w("1.2 × 10⁻⁹ (B) of the carrier's width. The absorption term −(f1 − f2)(W1 − W2) is not sign-definite pair by pair")
    w("when the two partners are close in frequency (the tetrahedron weights W1, W2 overlap); docs/v4/ANALYSIS_STEPS.md, step T3 §1")
    w("corrected. Other definitions: `analysis/t4/t4_{A,B}_summary.json`.\n")
    w("### 2d. κ_P carried by each family, next to κ_P carried by modes with η > 0.1 (all modes; %)\n")
    w("| definition | model | neither | pol | rot | mixed | η > 0.1 (group-averaged) |")
    w("|---|---|---|---|---|---|---|")
    for d in DEFS:
        for m in "AB":
            k = S[m]["definitions"][d]["kappa_P_share_by_family"]
            w(f"| {d} | {m} | " + " | ".join(pct(k[f]) for f in FAMS) + f" | {S[m]['share_kP_eta_gt_0.1'] * 100:.5f} |")
    w("")
    w("### 2e. Sensitivity to TRISYM, 8³\n")
    w("Same labels and settings, TRISYM = 1 against TRISYM = 0; each run's own 50 % carriers. Difference in percentage")
    w("points (TRISYM 1 − TRISYM 0).\n")
    w("| model | κ_P (W/mK) TRISYM 1 / 0 | carriers 1 / 0 (common) | neither / one / both soft, primary c0.5, TRISYM 1 (%) | same, TRISYM 0 (%) | largest abs. difference over all shares and definitions (pp) | at |")
    w("|---|---|---|---|---|---|---|")
    for m, s in (("A", senA), ("B", senB)):
        r = s["definitions"]["primary_c0.5"]
        w(f"| {m} | {s['kappa_P'][0]:.4f} / {s['kappa_P'][1]:.4f} | {s['carriers'][0]} / {s['carriers'][1]} "
          f"({s['carriers_common']}) | " + " / ".join(pct(r[f'nsoft{k}']['trisym1']) for k in "012") + " | "
          + " / ".join(pct(r[f'nsoft{k}']['trisym0']) for k in "012") + f" | {s['max_abs_diff_pp']:.3f} | {s['max_at']} |")
    esc = senA["escalate"] or senB["escalate"]
    w(f"\nEscalation rule (> 2 pp): **{'triggered' if esc else 'not triggered'}**"
      + ("" if esc else "; the 16³ TRISYM = 0 run goes on the list of the analysis plan") + "\n")

    # figure data: one per model, the carriers (primary c0.5)
    for m, M in (("A", mA), ("B", mB)):
        f = M[(M.definition == "primary_c0.5") & M.carrier].copy()
        f["kappa_share"] = f.kappa_mode_W_mK / S[m]["kappa_P_modesum"]
        cols = ["q_index", "branch", "q", "omega_cm1", "gamma_avg_cm1", "eta", "kappa_mode_W_mK", "kappa_share",
                "group_rank", "carrier_family", "nsoft0_frac", "nsoft1_frac", "nsoft2_frac",
                "nsoft0_decay", "nsoft0_absorption", "nsoft1_decay", "nsoft1_absorption",
                "nsoft2_decay", "nsoft2_absorption"]
        f.sort_values("omega_cm1")[cols].to_csv(figdir / f"t4_fig_{m}.csv", index=False, float_format="%.8g")

    body = out.read_text() if out.exists() else ""
    marker = "## 2. Results (opened here first)"
    if marker in body:
        body = body[: body.index(marker)]
    out.write_text(body + "\n".join(L) + "\n" + chk.get("tail", ""))
    print(out, figdir / "t4_fig_A.csv", figdir / "t4_fig_B.csv")


if __name__ == "__main__":
    main()
