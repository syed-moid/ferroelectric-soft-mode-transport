#!/usr/bin/env python3
"""Write the T5 report tables (Markdown) from the t5_control_analysis outputs and a checks JSON.
Tables only; no interpretation.

Usage: t5_report.py T5_A_JSON T5_B_JSON CHECKS_JSON OUT_MD [T5b]
With T5b: header for the T5b report tables (Markdown); bands are read from the JSON ("band_edges_cm1"; T5 default
105/145/160) and a section with the frequency of each class is added when the JSON carries it.
CHECKS_JSON: {"rows": [[run, labels sha, time s, max RSS GB, sum bins, total vs .result, outputs vs reference], ...]}
"""
import json
import sys

BANDS, BANDTXT = [], {}


def set_bands(edges):
    BANDS[:] = [f"<{edges[0]:g}"] + [f"{a:g}-{b:g}" for a, b in zip(edges, edges[1:])] + [f">={edges[-1]:g}"]
    txt = [f"Ω < {edges[0]:g}"] + [f"{a:g} ≤ Ω < {b:g}" for a, b in zip(edges, edges[1:])] + [f"Ω ≥ {edges[-1]:g}"]
    BANDTXT.clear()
    BANDTXT.update(zip(BANDS, txt))
FAMS = ["neither", "pol", "rot"]
PARTS = [("total", "decay + absorption"), ("decay_of_decay", "decay"), ("absorption_of_absorption", "absorption")]


def pct(x):
    return "—" if x is None else f"{100 * x + 0.0:.2f}".replace("-0.00", "0.00")


def pk(x):  # per-mode share, in % per 1000 modes
    return "—" if x is None else f"{1e5 * x:.3f}"


def rat(x):
    return "—" if x is None else f"{x:.2f}"


def main():
    D = {"A": json.load(open(sys.argv[1])), "B": json.load(open(sys.argv[2]))}
    chk = json.load(open(sys.argv[3]))
    out = sys.argv[4]
    t5b = len(sys.argv) > 5 and sys.argv[5] == "T5b"
    set_bands(D["A"].get("band_edges_cm1", [105.0, 145.0, 160.0]))
    assert D["B"].get("band_edges_cm1", [105.0, 145.0, 160.0]) == D["A"].get("band_edges_cm1", [105.0, 145.0, 160.0])
    L = []
    w = L.append
    if t5b:
        w("# T5b report: control inside the lowest band (" + chk.get("date", "") + ")\n")
        w("the analysis plan; band edges fixed before the labels in docs/v4/ANALYSIS_STEPS.md, step T5b (terciles of model A's rot")
        w("frequencies below 105 cm⁻¹, rounded to 5 cm⁻¹; bands above 105 merged). Same quantities as T5")
        w("(the T5 report tables (Markdown)) plus the frequency of each class. Descriptive; no pass criterion. Tables only;")
        w("the author writes the interpretation.\n")
        w("## 1. Runs and checks\n")
        w("Node, patched anphon a0a5d30c…, OMP 12, TRISYM = 1, 16³, 300 K, production inputs (rta.in identical to the")
        w("T4 runs); `t3/run_t5b.sh`, log `queue/incoming/node_b/t5b/run_t5b.log`; outputs pulled to")
        w("`queue/incoming/node_b/t5b/`, sha256 of every file checked. Label tables `t5b_control/t5b_labels_{A,B}_m16.txt`")
        w("(label = 4 × family + band), invariance check 0 violations (`analysis/t5b_label_invariance.txt`).\n")
    else:
        w("# T5 report: frequency-matched control (21-Nov-2026)\n")
        w("the analysis plan; criteria and band edges fixed before the labels in docs/v4/ANALYSIS_STEPS.md, step T5. Descriptive; no")
        w("pass criterion. Tables only; the author writes the interpretation.\n")
        w("## 1. Runs and checks\n")
        w("Node, patched anphon a0a5d30c…, OMP 12, TRISYM = 1, 16³, production inputs (rta.in identical to the T4 runs);")
        w("`t3/run_t5.sh`, log `queue/incoming/node_b/t5/run_t5.log`; outputs pulled to `queue/incoming/node_b/t5/`,")
        w("sha256 of every file checked. Label tables `t5_control/t5_labels_{A,B}_m16.txt` (label = 4 × family + band),")
        w("invariance check 0 violations (`analysis/t5_label_invariance.txt`).\n")
    w("| run | label table sha256 | time | max RSS | Σ bins = total, per carrier group (max rel.) | total = `.result` (max rel.) | outputs vs reference |")
    w("|---|---|---|---|---|---|---|")
    for r in chk["rows"]:
        w("| " + " | ".join(str(x) for x in r) + " |")
    w("")
    w("Carriers: the T4 carriers (primary c0.5, 16³, frozen 50 % rule): A " + str(D["A"]["carriers"]) + ", B "
      + str(D["B"]["carriers"]) + " irreducible modes; κ_P weights equal to T4 (checked, rel. 10⁻⁸).\n")

    w("## 2. Definitions\n")
    w("- S_c: κ_P-weighted mean over the carriers of the share of the carrier's group-averaged linewidth from processes")
    w("  with at least one partner in class c (a process with both partners in c counted once). Shares of the total width")
    w("  (decay + absorption); of the decay part (decay processes / decay width); of the absorption part (absorption")
    w("  processes / absorption width).")
    w("- N_c: number of modes in c on the full 16³ mesh (star-weighted; 61 440 modes).")
    w("- s_c = S_c / N_c, given in **% per 1000 modes** (S_c in % × 1000 / N_c).")
    w("- Per band: rot against non-rot (neither ∪ pol; at least one partner in the union) and against neither alone;")
    w("  ratio = s_rot / s_nonrot (s_rot / s_neither).")
    w("- Bands (cm⁻¹): " + "; ".join(f"{i} {BANDTXT[b]}" for i, b in enumerate(BANDS)) + ".\n")

    w("## 3. Share of the carriers' linewidth by partner class (family × band)\n")
    for m in "AB":
        c = D[m]["classes"]
        w(f"### Model {m}\n")
        w("| family | band (cm⁻¹) | N_c | S_c total (%) | of which decay (%) | of which absorption (%) | S_c of decay part (%) | S_c of absorption part (%) | s_c total | s_c decay | s_c absorption |")
        w("|---|---|---|---|---|---|---|---|---|---|---|")
        for f in FAMS:
            for b in BANDS:
                x = c[f"{f}|{b}"]
                S, p = x["S"], x["per_mode"]
                w(f"| {f} | {BANDTXT[b]} | {x['N_modes']} | {pct(S['total'])} | {pct(S['decay_of_total'])} | "
                  f"{pct(S['absorption_of_total'])} | {pct(S['decay_of_decay'])} | {pct(S['absorption_of_absorption'])} | "
                  f"{pk(p['total'])} | {pk(p['decay_of_decay'])} | {pk(p['absorption_of_absorption'])} |")
        w("\ns_c in % per 1000 modes. Classes overlap (one process can have partners in two classes), so the S_c do not")
        w("sum to 100 %.\n")

    w("## 4. Per band: rot partners against non-rot partners, per mode\n")
    for key, name in PARTS:
        w(f"### {name[0].upper() + name[1:]}" + ("" if key == "total" else f" (share of the {name} part)") + "\n")
        w("| band (cm⁻¹) | model | N rot | s rot | N non-rot | s non-rot | s rot / s non-rot | N neither | s neither | s rot / s neither |")
        w("|---|---|---|---|---|---|---|---|---|---|")
        for b in BANDS:
            for m in "AB":
                e = D[m]["bands"][b]
                w(f"| {BANDTXT[b]} | {m} | {e['rot']['N_modes']} | {pk(e['rot']['per_mode'][key])} | "
                  f"{e['nonrot']['N_modes']} | {pk(e['nonrot']['per_mode'][key])} | "
                  f"{rat(e['ratio_rot_over_nonrot'][key])} | {e['neither']['N_modes']} | "
                  f"{pk(e['neither']['per_mode'][key])} | {rat(e['ratio_rot_over_neither'][key])} |")
        w("")
    ctl = "t5b" if t5b else "t5"
    w(f"s in % per 1000 modes; — where the class is empty. Data: `analysis/{ctl}/{ctl}_{{A,B}}.json`.\n")
    if any(x.get("omega_cm1") for m in "AB" for x in D[m]["classes"].values()):
        fs = lambda o: "— | — | —" if not o else f"{o['mean']:.2f} | {o['min']:.2f} | {o['max']:.2f}"  # noqa: E731
        w("## 5. Frequency of each class (full 16³ mesh, unweighted over the class's modes; cm⁻¹)\n")
        w("| model | family | band (cm⁻¹) | N | Ω mean | Ω min | Ω max |")
        w("|---|---|---|---|---|---|---|")
        for m in "AB":
            for f in FAMS:
                for b in BANDS:
                    x = D[m]["classes"][f"{f}|{b}"]
                    w(f"| {m} | {f} | {BANDTXT[b]} | {x['N_modes']} | {fs(x.get('omega_cm1'))} |")
        w("\nPer-band unions (as in §4):\n")
        w("| model | partners | band (cm⁻¹) | N | Ω mean | Ω min | Ω max |")
        w("|---|---|---|---|---|---|---|")
        for m in "AB":
            for b in BANDS:
                e = D[m]["bands"][b]
                for nm in ("rot", "nonrot", "neither"):
                    w(f"| {m} | {'non-rot (neither ∪ pol)' if nm == 'nonrot' else nm} | {BANDTXT[b]} | "
                      f"{e[nm]['N_modes']} | {fs(e[nm].get('omega_cm1'))} |")
        w("")
    open(out, "w").write("\n".join(L))
    print(out)


if __name__ == "__main__":
    main()
