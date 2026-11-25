#!/usr/bin/env python3
"""Write the P3 report tables (Markdown) from the T4/T5/T5b analysis outputs of
model A, model B and the archived 2015 fit. Tables only; no interpretation. The table rules of the four observables
are those of docs/v4/ANALYSIS_STEPS.md, step P3 (copied from the analysis plan), applied as written.

Usage: p3_report.py CONFIG_JSON OUT_MD
CONFIG_JSON: {"date": ..., "checks": [[run, labels sha, time, RSS, Σ bins, total vs .result, outputs vs base], ...],
              "models": {"A": {...}, "B": {...}, "P": {...}}} with, per model: name, t4_summary, t5, t5b, sampling,
              labels_csv, rta_dir, prefix, eta_csv, eta_run, soft_csv, soft_run.
"""
import json
import sys
from pathlib import Path

import pandas as pd

FAMS4 = ["neither", "pol", "rot", "mixed"]
FAMS3 = ["neither", "pol", "rot"]
DEFS = ["primary_c0.5", "secondary_c0.5", "primary_c0.4", "primary_c0.6"]
ROT_WIN = (0.1059, 0.1659)     # the analysis plan: 13.59 % ± 3 pp (model A)
NOSOFT_WIN = (0.7313, 0.7913)  # the analysis plan: 76.13 % ± 3 pp (model A)


def pct(x):
    return "—" if x is None else f"{100 * x + 0.0:.2f}".replace("-0.00", "0.00")


def pk(x):  # per-mode share, % per 1000 modes
    return "—" if x is None else f"{1e5 * x:.3f}"


def rat(x):
    return "undefined" if x is None else f"{x:.2f}"


def ratio_class(r):
    if r is None:
        return "undefined"
    return "≥ 2.0" if r >= 2.0 else ("1.3 to < 2.0" if r >= 1.3 else "< 1.3")


def kl(d, prefix, ext):
    row = [ln for ln in open(Path(d) / f"{prefix}.{ext}") if not ln.startswith("#") and ln.strip()][0].split()
    return float(row[1])


def band_txt(b):
    if b.startswith("<"):
        return f"Ω < {b[1:]}"
    if b.startswith(">="):
        return f"Ω ≥ {b[2:]}"
    a, c = b.split("-")
    return f"{a} ≤ Ω < {c}"


def main():
    cfg = json.load(open(sys.argv[1]))
    out = Path(sys.argv[2])
    MS = list(cfg["models"])
    C = cfg["models"]
    T4 = {m: json.load(open(C[m]["t4_summary"])) for m in MS}
    T5 = {m: json.load(open(C[m]["t5"])) for m in MS}
    T5b = {m: json.load(open(C[m]["t5b"])) for m in MS}
    SMP = {m: json.load(open(C[m]["sampling"])) for m in MS}
    LAB = {m: pd.read_csv(C[m]["labels_csv"]) for m in MS}
    NM = {m: C[m]["name"] for m in MS}
    L = []
    w = L.append
    hdr = "| | " + " | ".join(NM[m] for m in MS) + " |\n|---|" + "---|" * len(MS)

    w(f"# P3: the archived 2015 fit under the present analysis ({cfg['date']})\n")
    w("Rules recorded before any run of the archived fit (docs/v4/ANALYSIS_STEPS.md, step P3).")
    w("**P3 is a fit-sensitivity test**: the anharmonic fit changes and the force data do not; many")
    w("force constants change together, so it is not a test of the R frequency alone. Each observable is reported")
    w("separately; there is no overall pass or fail. Tables only; the author writes the interpretation.\n")
    w("Archived fit (force constants of the 2015 study), distributed `BORN`, NONANALYTIC = 3, q₁ 12³, interpolation 2³, 300 K, transport 16³,")
    w("TRISYM = 1. Only the anharmonic force constants differ from model A (harmonic part identical, P0).\n")

    w("## 1. Runs and checks\n")
    w("OMP 12, one job at a time. Base runs: c5b015c (anphon 266b507a…); patched runs: a0a5d30c… as for T4/T5/T5b;")
    w("outputs checked against their SHA256SUMS. Labels: frozen rules (docs/v4/ANALYSIS_STEPS.md, step T2), same scripts")
    w("as models A and B, invariance 0 violations.\n")
    w("| run | label table sha256 | time | max RSS | Σ bins = total, per carrier group (max rel.) | total = `.result` (max rel.) | outputs vs base run |")
    w("|---|---|---|---|---|---|---|")
    for r in cfg["checks"]:
        w("| " + " | ".join(str(x) for x in r) + " |")
    w("")

    # ---- the four observables
    w("## 2. The four observables (300 K, 16³, primary definition c_min 0.5)\n")
    w("Rules as recorded before the runs (docs/v4/ANALYSIS_STEPS.md, step P3). Windows: rot share 10.59–16.59 %, no-soft-partner share")
    w("73.13–79.13 % (model A ± 3 percentage points; a practical tolerance, not an error bar). Ratio classes: ≥ 2.0;")
    w("1.3 to < 2.0; < 1.3.\n")
    w("| observable | " + " | ".join(NM[m] for m in MS) + " |\n|---|" + "---|" * len(MS))
    rot = {m: T4[m]["definitions"]["primary_c0.5"]["kappa_P_share_by_family"]["rot"] for m in MS}
    pol = {m: T4[m]["definitions"]["primary_c0.5"]["kappa_P_share_by_family"]["pol"] for m in MS}
    ns0 = {m: T4[m]["definitions"]["primary_c0.5"]["kappa_weighted_nsoft_frac"]["0"] for m in MS}

    def win(x, lo, hi, ref):
        return ("inside" if lo - 1e-12 <= x <= hi + 1e-12 else "outside") + f" (Δ vs A {100 * (x - ref):+.2f} pp)".replace("-", "−")
    w("| κ_P carried by the rotational family (%) | " + " | ".join(pct(rot[m]) for m in MS) + " |")
    w("| — class (window 10.59–16.59 %) | " + " | ".join(win(rot[m], *ROT_WIN, rot["A"]) for m in MS) + " |")
    w("| — κ_P carried by the polar family (%), for the separate rot vs pol statement | "
      + " | ".join(pct(pol[m]) for m in MS) + " |")
    w("| carrier-weighted linewidth share with no soft partner (%) | " + " | ".join(pct(ns0[m]) for m in MS) + " |")
    w("| — class (window 73.13–79.13 %) | " + " | ".join(win(ns0[m], *NOSOFT_WIN, ns0["A"]) for m in MS) + " |")
    w("| — above 50 % | " + " | ".join("yes" if ns0[m] > 0.5 else "no" for m in MS) + " |")
    for b, src in (("<105", T5), ("<80", T5b)):
        e = {m: src[m]["bands"][b] for m in MS}
        r = {m: e[m]["ratio_rot_over_nonrot"]["total"] for m in MS}
        w(f"| per-mode ratio rot / non-rot, {band_txt(b)} cm⁻¹ | " + " | ".join(rat(r[m]) for m in MS) + " |")
        w("| — class | " + " | ".join(ratio_class(r[m]) for m in MS) + " |")
        w("| — N rot / N non-rot | " + " | ".join(f"{e[m]['rot']['N_modes']} / {e[m]['nonrot']['N_modes']}"
                                                   for m in MS) + " |")
        w("| — S rot / S non-rot (%, unnormalized) | " + " | ".join(
            f"{pct(e[m]['rot']['S']['total'])} / {pct(e[m]['nonrot']['S']['total'])}" for m in MS) + " |")
        w("| — sampling test | " + " | ".join(
            "sensitive to sampling" if SMP[m]["bands"][b]["sensitive_to_sampling"] else "not sensitive"
            for m in MS) + " |")
    w("\nRatio = s_rot / s_non-rot with s = S / N (non-rot = neither ∪ pol), share of the total width (decay +")
    w("absorption), as in the T5 and T5b reports. Sampling test: §8.\n")

    # ---- bands side by side
    w("## 3. Band comparisons, one definition (primary c0.5), the three fits side by side\n")
    w("| band (cm⁻¹) | fit | N rot | N non-rot | S rot (%) | S non-rot (%) | s rot | s non-rot | ratio | class | sampling |")
    w("|---|---|---|---|---|---|---|---|---|---|---|")
    for b, src in (("<105", T5), ("<80", T5b), ("80-95", T5b), ("95-105", T5b)):
        for m in MS:
            e = src[m]["bands"][b]
            r = e["ratio_rot_over_nonrot"]["total"]
            w(f"| {band_txt(b)} | {NM[m]} | {e['rot']['N_modes']} | {e['nonrot']['N_modes']} | "
              f"{pct(e['rot']['S']['total'])} | {pct(e['nonrot']['S']['total'])} | {pk(e['rot']['per_mode']['total'])} | "
              f"{pk(e['nonrot']['per_mode']['total'])} | {rat(r)} | {ratio_class(r)} | "
              + ("sensitive to sampling" if SMP[m]["bands"][b]["sensitive_to_sampling"] else "not sensitive") + " |")
    w("\ns in % per 1000 modes. The ±3-point windows apply to §2 rows 1–2 only; the band ratios carry the classes.\n")

    # ---- kappa, eta, soft modes, rot family
    w("## 4. κ_P, κ_C, η (300 K, 16³, q₁ 12³)\n")
    eta = {}
    for m in MS:
        d = pd.read_csv(C[m]["eta_csv"])
        eta[m] = d[d.run == C[m]["eta_run"]].iloc[0]
    w(hdr)
    w("| κ_P (`.kl`, xx; W m⁻¹ K⁻¹) | " + " | ".join(f"{kl(C[m]['rta_dir'], C[m]['prefix'], 'kl'):.4f}" for m in MS) + " |")
    w("| κ_C (`.kl_coherent`, xx) | " + " | ".join(f"{kl(C[m]['rta_dir'], C[m]['prefix'], 'kl_coherent'):.4f}"
                                                 for m in MS) + " |")
    w("| κ_P mode sum (group-averaged widths) | " + " | ".join(f"{eta[m].kappa_P_modesum:.5f}" for m in MS) + " |")
    for t in ("0.02", "0.05", "0.1"):
        w(f"| share of κ_P with η > {t} (%) | " + " | ".join(f"{100 * eta[m][f'share_kP_eta_gt_{t}']:.5g}" for m in MS) + " |")
    w("| largest group-averaged η (at q index, branch; Ω cm⁻¹) | " + " | ".join(
        f"{eta[m].eta_max:.4f} ({int(eta[m].eta_max_q_index)}, {int(eta[m].eta_max_branch)}; {eta[m].eta_max_omega_cm1:.2f})"
        for m in MS) + " |")
    w("\nη = group-averaged HWHM / Ω (`scripts/kp_eta_distribution.py --gamma=avg`).\n")

    w("## 5. Soft-mode frequencies on q₁ 12³ (300 K; cm⁻¹; interpolation 2³)\n")
    w("| | " + " | ".join(NM[m] for m in MS) + " |\n|---|" + "---|" * len(MS))
    sm = {}
    for m in MS:
        d = pd.read_csv(C[m]["soft_csv"])
        sm[m] = d[(d.run == C[m]["soft_run"]) & (d.T_K == 300)].iloc[0]
    for k in ("Gamma", "X", "M", "R"):
        w(f"| Ω_{k.replace('Gamma', 'Γ')} | " + " | ".join(f"{sm[m][k]:.2f}" for m in MS) + " |")
    w("| lowest converged T (K) | " + " | ".join(str(int(sm[m].lowest_contiguous_converged_T)) for m in MS) + " |")
    w("\nExtraction rule of `scripts/scph_softmodes.py` (Γ 4th branch at k = 0; X, M, R lowest branch).\n")

    w("## 6. Rotational family: modes per band and frequency range (full 16³ mesh, primary c0.5)\n")
    bands_all = ["<80", "80-95", "95-105", "105-145", "145-160", ">=160"]
    w("| fit | " + " | ".join(band_txt(b) for b in bands_all) + " | total | Ω min | Ω max |")
    w("|---|" + "---|" * (len(bands_all) + 3))
    for m in MS:
        lab = LAB[m]
        rr = lab[lab["primary_c0.5"] == "rot"].Omega_cm1
        n = []
        for b in bands_all:
            src = T5b[m] if b in ("<80", "80-95", "95-105") else T5[m]
            n.append(str(src["classes"][f"rot|{b}"]["N_modes"]))
        w(f"| {NM[m]} | " + " | ".join(n) + f" | {len(rr)} | {rr.min():.2f} | {rr.max():.2f} |")
    w("\nFamily counts on the full 16³ mesh (61 440 modes), every definition (neither / pol / rot / mixed):\n")
    w("| definition | " + " | ".join(NM[m] for m in MS) + " |\n|---|" + "---|" * len(MS))
    for d in ("primary_c0.4", "primary_c0.5", "primary_c0.6", "secondary_c0.4", "secondary_c0.5", "secondary_c0.6"):
        w(f"| {d} | " + " | ".join(" / ".join(str(int((LAB[m][d] == f).sum())) for f in FAMS4) for m in MS) + " |")
    w("\nLabels: `t2_frozen/labels_{A,B}_m16.csv`, `t2_frozen/p15/labels_P_m16.csv`.\n")

    # ---- T4 layout
    w("## 7. Full tables in the layout of the T4, T5 and T5b reports\n")
    w("### 7a. Carriers (frozen 50 % rule) — T4 §2a\n")
    w(hdr)
    w("| κ_P mode sum (W/mK) | " + " | ".join(f"{T4[m]['kappa_P_modesum']:.4f}" for m in MS) + " |")
    w("| carrier modes (irreducible) / degenerate groups | " + " | ".join(
        f"{T4[m]['carriers']} / {T4[m]['carrier_groups']}" for m in MS) + " |")
    w("| share of κ_P carried by the carriers | " + " | ".join(f"{pct(T4[m]['carrier_kappa_share'])} %" for m in MS) + " |")
    for d in DEFS:
        w(f"| carriers by family, {d} (neither / pol / rot / mixed) | " + " | ".join(
            " / ".join(str(T4[m]["definitions"][d]["carrier_family_counts"][f]) for f in FAMS4) for m in MS) + " |")
    w("\n### 7b. Linewidth of the carriers by number of soft partners — T4 §2b\n")
    w("| definition | fit | neither soft (%) | one soft (%) | both soft (%) | absorption part (%) |")
    w("|---|---|---|---|---|---|")
    for d in DEFS:
        for m in MS:
            x = T4[m]["definitions"][d]
            n = x["kappa_weighted_nsoft_frac"]
            w(f"| {d} | {NM[m]} | {pct(n['0'])} | {pct(n['1'])} | {pct(n['2'])} | {pct(x['kappa_weighted_absorption_share'])} |")
    w("\nMedians over the carriers (primary c0.5), neither / one / both soft: " + "; ".join(
        f"{NM[m]} " + " / ".join(pct(T4[m]["definitions"]["primary_c0.5"]["median_nsoft_frac"][k]) for k in "012") + " %"
        for m in MS) + ".\n")
    w("### 7c. Linewidth of the carriers by partner family pair (primary c0.5; %) — T4 §2c\n")
    pairs = list(T4["A"]["definitions"]["primary_c0.5"]["kappa_weighted_pair_frac"])
    w("| partner pair | " + " | ".join(NM[m] for m in MS) + " |\n|---|" + "---|" * len(MS))
    for p in pairs:
        v = [T4[m]["definitions"]["primary_c0.5"]["kappa_weighted_pair_frac"][p] for m in MS]
        if any(v):
            w(f"| {p} | " + " | ".join(pct(x) for x in v) + " |")
    w("\nPairs absent from the table are zero in all three fits. Values shown as 0.00 may be small negatives (absorption")
    w("term not sign-definite pair by pair; T4 report §2c): " + "; ".join(
        f"{NM[m]} pol/pol {T4[m]['definitions']['primary_c0.5']['kappa_weighted_pair_frac']['pol/pol']:.2e}" for m in MS) + ".\n")
    w("### 7d. κ_P carried by each family, next to κ_P with η > 0.1 (all modes; %) — T4 §2d\n")
    w("| definition | fit | neither | pol | rot | mixed | η > 0.1 (group-averaged) |")
    w("|---|---|---|---|---|---|---|")
    for d in DEFS:
        for m in MS:
            k = T4[m]["definitions"][d]["kappa_P_share_by_family"]
            w(f"| {d} | {NM[m]} | " + " | ".join(pct(k[f]) for f in FAMS4) + f" | {T4[m]['share_kP_eta_gt_0.1'] * 100:.5f} |")
    w("")

    for tag, D, title in (("T5", T5, "7e. Band control, edges 105 / 145 / 160 cm⁻¹ — T5 §3–4"),
                          ("T5b", T5b, "7f. Lowest band divided at 80 / 95 cm⁻¹ — T5b §3–5")):
        w(f"### {title}\n")
        bands = list(D["A"]["bands"])
        for m in MS:
            assert list(D[m]["bands"]) == bands
        w("Share of the carriers' linewidth by partner class (family × band); S in %, s in % per 1000 modes.\n")
        w("| fit | family | band (cm⁻¹) | N_c | S_c total | of which decay | of which absorption | S_c of decay part | S_c of absorption part | s_c total | s_c decay | s_c absorption |")
        w("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for m in MS:
            for f in FAMS3:
                for b in bands:
                    x = D[m]["classes"][f"{f}|{b}"]
                    S, p = x["S"], x["per_mode"]
                    w(f"| {NM[m]} | {f} | {band_txt(b)} | {x['N_modes']} | {pct(S['total'])} | {pct(S['decay_of_total'])} | "
                      f"{pct(S['absorption_of_total'])} | {pct(S['decay_of_decay'])} | {pct(S['absorption_of_absorption'])} | "
                      f"{pk(p['total'])} | {pk(p['decay_of_decay'])} | {pk(p['absorption_of_absorption'])} |")
        w("")
        for key, nm in (("total", "decay + absorption"), ("decay_of_decay", "decay part"),
                        ("absorption_of_absorption", "absorption part")):
            w(f"Per band, rot against non-rot, {nm}:\n")
            w("| band (cm⁻¹) | fit | N rot | s rot | N non-rot | s non-rot | s rot / s non-rot | N neither | s neither | s rot / s neither |")
            w("|---|---|---|---|---|---|---|---|---|---|")
            for b in bands:
                for m in MS:
                    e = D[m]["bands"][b]
                    w(f"| {band_txt(b)} | {NM[m]} | {e['rot']['N_modes']} | {pk(e['rot']['per_mode'][key])} | "
                      f"{e['nonrot']['N_modes']} | {pk(e['nonrot']['per_mode'][key])} | "
                      f"{rat(e['ratio_rot_over_nonrot'][key])} | {e['neither']['N_modes']} | "
                      f"{pk(e['neither']['per_mode'][key])} | {rat(e['ratio_rot_over_neither'][key])} |")
            w("")
        if all(x.get("omega_cm1") is not None or x["N_modes"] == 0 for m in MS for x in D[m]["classes"].values()):
            fs = lambda o: "— | — | —" if not o else f"{o['mean']:.2f} | {o['min']:.2f} | {o['max']:.2f}"  # noqa: E731
            w("Frequency of each class (full 16³ mesh, unweighted; cm⁻¹):\n")
            w("| fit | family | band (cm⁻¹) | N | Ω mean | Ω min | Ω max |")
            w("|---|---|---|---|---|---|---|")
            for m in MS:
                for f in FAMS3:
                    for b in bands:
                        x = D[m]["classes"][f"{f}|{b}"]
                        w(f"| {NM[m]} | {f} | {band_txt(b)} | {x['N_modes']} | {fs(x.get('omega_cm1'))} |")
            w("")

    # ---- sampling
    w("## 8. Sampling test (count-only; `scripts/p3_sampling.py`)\n")
    w("Unit moved: one star-group (degenerate group × star of the 16³ mesh). The nearest rotational star-group, and")
    w("separately the nearest non-rotational one, is moved across each edge of the band (counts change, shares")
    w("unchanged) and the ratio recomputed. A move that takes the ratio across 2.0 or 1.3 marks the band \"sensitive to")
    w("sampling\"; a run with shifted edges would then be proposed, not run.\n")
    w("| band (cm⁻¹) | fit | ratio | moved | star-group Ω (cm⁻¹) | size | direction | edge | distance | N after (rot / non-rot) | ratio after | crosses 2.0 or 1.3 |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for b in ("<105", "<80", "80-95", "95-105"):
        for m in MS:
            rec = SMP[m]["bands"][b]
            for mv in rec["moves"]:
                w(f"| {band_txt(b)} | {NM[m]} | {rat(rec['ratio'])} | {mv['class_moved']} | {mv['star_group_Omega']:.3f} | "
                  f"{mv['size']} | {mv['direction']} | {mv['edge']:g} | {mv['distance_cm1']:.3f} | "
                  f"{mv['N_after']['rot']} / {mv['N_after']['nonrot']} | {rat(mv['ratio_after'])} | "
                  f"{'yes' if mv['crosses_2.0_or_1.3'] else 'no'} |")
    w("\nStar-groups within 2 cm⁻¹ of each edge: `analysis/p15/p3_sampling_{A,B,P}.json` (\"edges\").\n")
    out.write_text("\n".join(L))
    print(out)


if __name__ == "__main__":
    main()
