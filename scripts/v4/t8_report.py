#!/usr/bin/env python3
"""Write the share part of the T8 report tables (Markdown): every row of T4 §2a–2d for one model (A by default; B for T8b), 16³, 300 K,
with TRISYM = 1 (production) and TRISYM = 0 (full sum) side by side, each with its own 50 % carriers. Derived from
t7_report.py (same rows). Tables only; no interpretation.

Usage: t8_report.py S_TRISYM1 S_TRISYM0 OUT_MD_PART [MODEL]  (MODEL A by default; B for T8b, the analysis plan)
"""
import json
import sys

FAMS = ["neither", "pol", "rot", "mixed"]
DEFS = ["primary_c0.5", "secondary_c0.5", "primary_c0.4", "primary_c0.6"]


def pct(x):
    return f"{100 * x + 0.0:.2f}".replace("-0.00", "0.00")


def row(name, model, a, b, digits=2):
    d = 100 * (b - a)
    z = "0." + "0" * digits
    c = [f"{100 * a + 0.0:.{digits}f}", f"{100 * b + 0.0:.{digits}f}", f"{d + 0.0:+.{digits}f}"]
    c = [z if x in ("-" + z, "+" + z) else x for x in c]          # no signed zeros
    return f"| {name} | {MB(model)} | " + " | ".join(c) + " |"


def MB(m):
    return m


def main():
    M = sys.argv[4] if len(sys.argv) > 4 else "A"
    S = {(M, 16): json.load(open(sys.argv[1])), (M, 20): json.load(open(sys.argv[2]))}
    out = sys.argv[3]
    L = []
    w = L.append
    H = ("| quantity | model | TRISYM 1 (%) | TRISYM 0 (%) | TRISYM 0 − TRISYM 1 (pp) |\n"
         "|---|---|---|---|---|")
    w("## 3. T4 shares under both settings (rows of T4 §2a–2d)\n")
    w(f"Same definitions as T4 (κ_P-weighted means over each run's own 50 % carriers; group-averaged widths; model {M},")
    w(f"16³, q₁ 12³, 300 K; frozen labels t2_frozen/t3_code32_{M}_m16.txt for both runs).\n")
    w("### 3a. Carriers\n")
    w("| quantity | model | TRISYM 1 | TRISYM 0 |\n|---|---|---|---|")
    for m in M:
        a, b = S[(m, 16)], S[(m, 20)]
        w(f"| κ_P mode sum (W/mK) | {MB(m)} | {a['kappa_P_modesum']:.4f} | {b['kappa_P_modesum']:.4f} |")
        w(f"| carrier modes (irreducible) / degenerate groups | {MB(m)} | {a['carriers']} / {a['carrier_groups']} | "
          f"{b['carriers']} / {b['carrier_groups']} |")
    for d in DEFS:
        for m in M:
            ca, cb = (S[(m, n)]["definitions"][d]["carrier_family_counts"] for n in (16, 20))
            w(f"| carriers by family, {d} (neither / pol / rot / mixed) | {MB(m)} | "
              + " / ".join(str(ca[f]) for f in FAMS) + " | " + " / ".join(str(cb[f]) for f in FAMS) + " |")
    w("")
    w(H)
    for m in M:
        w(row("share of κ_P carried by the carriers", m, S[(m, 16)]["carrier_kappa_share"],
              S[(m, 20)]["carrier_kappa_share"]))
    w("")
    w("### 3b. Linewidth of the carriers by number of soft partners\n")
    w(H)
    for d in DEFS:
        for m in M:
            x, y = (S[(m, n)]["definitions"][d] for n in (16, 20))
            for k, nm in (("0", "neither soft"), ("1", "one soft"), ("2", "both soft")):
                w(row(f"{d}: {nm}", m, x["kappa_weighted_nsoft_frac"][k], y["kappa_weighted_nsoft_frac"][k]))
            w(row(f"{d}: absorption part", m, x["kappa_weighted_absorption_share"],
                  y["kappa_weighted_absorption_share"]))
    for m in M:
        x, y = (S[(m, n)]["definitions"]["primary_c0.5"]["median_nsoft_frac"] for n in (16, 20))
        for k, nm in (("0", "neither soft"), ("1", "one soft"), ("2", "both soft")):
            w(row(f"primary_c0.5 median: {nm}", m, x[k], y[k]))
    w("")
    w("### 3c. Linewidth of the carriers by partner family pair (primary c0.5)\n")
    w(H)
    pairs = list(S[(M, 16)]["definitions"]["primary_c0.5"]["kappa_weighted_pair_frac"])
    for p in pairs:
        v = {k: S[k]["definitions"]["primary_c0.5"]["kappa_weighted_pair_frac"][p] for k in S}
        if any(v.values()):
            for m in M:
                w(row(p, m, v[(m, 16)], v[(m, 20)]))
    w("\nPairs absent from the table are zero in both runs.\n")
    w("### 3d. κ_P carried by each family, next to κ_P carried by modes with η > 0.1 (all modes)\n")
    w(H)
    for d in DEFS:
        for m in M:
            x, y = (S[(m, n)]["definitions"][d]["kappa_P_share_by_family"] for n in (16, 20))
            for f in FAMS:
                w(row(f"{d}: {f}", m, x[f], y[f]))
    for m in M:
        w(row("η > 0.1 (group-averaged)", m, S[(m, 16)]["share_kP_eta_gt_0.1"], S[(m, 20)]["share_kP_eta_gt_0.1"],
              digits=5))
    w("")
    open(out, "w").write("\n".join(L))
    print(out)


if __name__ == "__main__":
    main()
