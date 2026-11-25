#!/usr/bin/env python3
"""Write the results part of the T6 report tables (Markdown): every row of T4 §2a–2d, 16³ and 20³ side by side,
with the fixed reading of docs/v4/ANALYSIS_STEPS.md, step T6 (a share is mesh-stable if 16³ and 20³ differ by ≤ 2 pp; otherwise
both values and the range). Tables only; no interpretation.

Usage: t6_report.py S16_A S16_B S20_A S20_B OUT_MD_PART
"""
import json
import sys

FAMS = ["neither", "pol", "rot", "mixed"]
DEFS = ["primary_c0.5", "secondary_c0.5", "primary_c0.4", "primary_c0.6"]
TOL_PP = 2.0


def pct(x):
    return f"{100 * x + 0.0:.2f}".replace("-0.00", "0.00")


def row(name, model, a, b, digits=2):
    d = 100 * (b - a)
    st = "yes" if abs(d) <= TOL_PP else f"no; range {min(a, b) * 100:.{digits}f}–{max(a, b) * 100:.{digits}f} %"
    return (f"| {name} | {model} | {100 * a + 0.0:.{digits}f} | {100 * b + 0.0:.{digits}f} | "
            f"{d + 0.0:+.{digits}f} | {st} |").replace("-0.00 ", "0.00 ").replace("+-0.00", "0.00")


def main():
    S = {("A", 16): json.load(open(sys.argv[1])), ("B", 16): json.load(open(sys.argv[2])),
         ("A", 20): json.load(open(sys.argv[3])), ("B", 20): json.load(open(sys.argv[4]))}
    out = sys.argv[5]
    L = []
    w = L.append
    H = ("| quantity | model | 16³ (%) | 20³ (%) | 20³ − 16³ (pp) | mesh-stable (≤ 2 pp) |\n"
         "|---|---|---|---|---|---|")
    w("## 2. Results: 16³ and 20³ side by side (rows of T4 §2a–2d)\n")
    w("Same definitions as T4 (κ_P-weighted means over each mesh's own 50 % carriers; group-averaged widths; 300 K,")
    w("q₁ 12³, TRISYM = 1). Fixed reading: mesh-stable if |20³ − 16³| ≤ 2 pp; otherwise both values and the range.\n")
    w("### 2a. Carriers\n")
    w("| quantity | model | 16³ | 20³ |\n|---|---|---|---|")
    for m in "AB":
        a, b = S[(m, 16)], S[(m, 20)]
        w(f"| κ_P mode sum (W/mK) | {m} | {a['kappa_P_modesum']:.4f} | {b['kappa_P_modesum']:.4f} |")
        w(f"| carrier modes (irreducible) / degenerate groups | {m} | {a['carriers']} / {a['carrier_groups']} | "
          f"{b['carriers']} / {b['carrier_groups']} |")
    for d in DEFS:
        for m in "AB":
            ca, cb = (S[(m, n)]["definitions"][d]["carrier_family_counts"] for n in (16, 20))
            w(f"| carriers by family, {d} (neither / pol / rot / mixed) | {m} | "
              + " / ".join(str(ca[f]) for f in FAMS) + " | " + " / ".join(str(cb[f]) for f in FAMS) + " |")
    w("")
    w(H)
    for m in "AB":
        w(row("share of κ_P carried by the carriers", m, S[(m, 16)]["carrier_kappa_share"],
              S[(m, 20)]["carrier_kappa_share"]))
    w("")
    w("### 2b. Linewidth of the carriers by number of soft partners\n")
    w(H)
    for d in DEFS:
        for m in "AB":
            x, y = (S[(m, n)]["definitions"][d] for n in (16, 20))
            for k, nm in (("0", "neither soft"), ("1", "one soft"), ("2", "both soft")):
                w(row(f"{d}: {nm}", m, x["kappa_weighted_nsoft_frac"][k], y["kappa_weighted_nsoft_frac"][k]))
            w(row(f"{d}: absorption part", m, x["kappa_weighted_absorption_share"],
                  y["kappa_weighted_absorption_share"]))
    for m in "AB":
        x, y = (S[(m, n)]["definitions"]["primary_c0.5"]["median_nsoft_frac"] for n in (16, 20))
        for k, nm in (("0", "neither soft"), ("1", "one soft"), ("2", "both soft")):
            w(row(f"primary_c0.5 median: {nm}", m, x[k], y[k]))
    w("")
    w("### 2c. Linewidth of the carriers by partner family pair (primary c0.5)\n")
    w(H)
    pairs = list(S[("A", 16)]["definitions"]["primary_c0.5"]["kappa_weighted_pair_frac"])
    for p in pairs:
        v = {k: S[k]["definitions"]["primary_c0.5"]["kappa_weighted_pair_frac"][p] for k in S}
        if any(v.values()):
            for m in "AB":
                w(row(p, m, v[(m, 16)], v[(m, 20)]))
    w("\nPairs absent from the table are zero on both meshes in both models.\n")
    w("### 2d. κ_P carried by each family, next to κ_P carried by modes with η > 0.1 (all modes)\n")
    w(H)
    for d in DEFS:
        for m in "AB":
            x, y = (S[(m, n)]["definitions"][d]["kappa_P_share_by_family"] for n in (16, 20))
            for f in FAMS:
                w(row(f"{d}: {f}", m, x[f], y[f]))
    for m in "AB":
        w(row("η > 0.1 (group-averaged)", m, S[(m, 16)]["share_kP_eta_gt_0.1"], S[(m, 20)]["share_kP_eta_gt_0.1"],
              digits=5))
    w("")
    n = sum(1 for x in L if x.startswith("| ") and x.endswith("| yes |"))
    nn = sum(1 for x in L if "| no; range" in x)
    w(f"Count over the share rows above: {n} mesh-stable, {nn} not mesh-stable.\n")
    open(out, "w").write("\n".join(L))
    print(out)


if __name__ == "__main__":
    main()
