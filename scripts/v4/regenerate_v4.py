#!/usr/bin/env python3
"""Regenerate the processed outputs of models A and B (parts i and ii) from the stored raw run directories.

RUNS_ROOT holds one directory per run, named as in data/provenance/v4/run_map.csv (column `run`), e.g.
RUNS_ROOT/A_partition_q12_m16_T300 is the raw output of job t4_a_m16_code32. The raw outputs are in the Zenodo
archives named in the run map (column `archive`); unpack them and link or rename the run directories accordingly.

Outputs (OUT = data/processed/v4 by default):
  OUT/a/  model A and its refits (part i)       OUT/b/  model B, and tables that hold A and B side by side (part ii)
Within each: labels/ (frozen family labels and label tables of the patched runs), t4/ (carrier analysis; m20/, T600/,
trisym0/, m8/), t5/, t5b/ (frequency-band analysis), t1/ (transport summaries), kp_eta_* (damping-ratio shares),
reports/ (Markdown tables of the report writers), figures/.

Usage (repository root): regenerate_v4.py RUNS_ROOT [OUT] [--source-bands DIR] [--a3-transport CSV] [--only STEP,...]
  --source-bands: data/SrTiO3/figures/bands of the unpacked model B dataset (for the soft-mode figure draft);
  --a3-transport: benchmark table A3_transport.csv (model A rows of the A-next-to-B comparison).
Steps: labels, t4, t5, kp_eta, t8, scph, t1, b4, census_tables, reports, counts, figures (all by default, in this
order).
"""
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable


def run(script, *args, env=None, stdout=None):
    e = dict(os.environ)
    e.update(env or {})
    cmd = [PY, str(HERE / script), *map(str, args)]
    print("+", " ".join(Path(c).name if i < 2 else str(c) for i, c in enumerate(cmd)), flush=True)
    with (open(stdout, "w") if stdout else open(os.devnull, "w")) as fh:
        subprocess.run(cmd, check=True, env=e, stdout=fh)


def one(d, pattern):
    hits = sorted(Path(d).glob(pattern))
    assert len(hits) == 1, (d, pattern, hits)
    return hits[0]


def main():
    argv = sys.argv[1:]
    opts = {}
    for k in ("--source-bands", "--a3-transport", "--only"):
        if k in argv:
            i = argv.index(k)
            opts[k] = argv[i + 1]
            del argv[i:i + 2]
    R = Path(argv[0])
    OUT = Path(argv[1]) if len(argv) > 1 else Path("data/processed/v4")
    steps = opts.get("--only", "labels,t4,t5,kp_eta,t8,scph,t1,b4,census_tables,reports,counts,figures").split(",")
    D = {m: OUT / m.lower() for m in "AB"}
    for d in D.values():
        for sub in ("labels/T600", "labels/m20", "labels/m8", "census", "t4/m20", "t4/T600", "t4/trisym0", "t4/m8/t1",
                    "t4/m8/t0", "t5", "t5b", "t1", "reports", "figures"):
            (d / sub).mkdir(parents=True, exist_ok=True)

    # ---- frozen family labels: census of the renormalized eigenvectors against the bare harmonic ones -------------
    cen = {  # tag: (model, renormalized evec, bare evec, label subdir, label tag)
        "A300": ("A", one(R / "A_evec_q12_m16_T300", "*.mesh.evec"), one(R / "A_harm_m16", "*.mesh.evec"), "", "A_m16"),
        "B300": ("B", one(R / "B_evec_q12_m16_T300", "*.mesh.evec"), one(R / "B_harm_m16", "*.mesh.evec"), "", "B_m16"),
        "A600": ("A", one(R / "A_evec_q12_m16_T600", "*.mesh.evec"), one(R / "A_harm_m16", "*.mesh.evec"), "T600", "A_m16_T600"),
        "B600": ("B", one(R / "B_evec_q12_m16_T600", "evec_*.mesh.evec"), one(R / "B_harm_m16", "*.mesh.evec"), "T600", "B_m16_T600"),
        "A20": ("A", one(R / "A_evec_q12_m20_T300", "*.mesh.evec"), one(R / "A_harm_m20", "*.mesh.evec"), "m20", "A_m20"),
        "B20": ("B", one(R / "B_evec_q12_m20_T300", "evec_*.mesh.evec"), one(R / "B_evec_q12_m20_T300", "harm_*.mesh.evec"), "m20", "B_m20"),
        "A8": ("A", one(R / "A_evec_q12_m8_T300", "*.mesh.evec"), one(R / "A_harm_m8", "*.mesh.evec"), "m8", "A_m8"),
        "B8": ("B", one(R / "B_evec_q12_m8_T300", "evec_*.mesh.evec"), one(R / "B_evec_q12_m8_T300", "harm_*.mesh.evec"), "m8", "B_m8"),
    }
    lab = {}
    for tag, (m, ren, bare, sub, ltag) in cen.items():
        lab[ltag] = D[m] / "labels" / sub / f"labels_{ltag}.csv"
        if "labels" in steps:
            run("t2_family_census.py", D[m] / "census" / "t2_census", f"{tag}={ren}:{bare}")
            run("t2_freeze_labels.py", D[m] / "census" / f"t2_census_{tag}.csv", ren, D[m] / "labels" / sub, ltag)
    t5lab = {m: D[m] / "labels" / f"t5_labels_{m}_m16.txt" for m in "AB"}
    t5blab = {m: D[m] / "labels" / f"t5b_labels_{m}_m16.txt" for m in "AB"}
    if "labels" in steps:
        for m in "AB":
            run("t5_control_labels.py", lab[f"{m}_m16"], t5lab[m])
            run("t5_control_labels.py", lab[f"{m}_m16"], t5blab[m], "80,95,105")

    # ---- carrier analysis (T4) on the patched runs -------------------------------------------------------------
    t4 = [  # (model, run, labels, outdir, env)
        *[(m, f"{m}_partition_q12_m16_T300", f"{m}_m16", "t4", {}) for m in "AB"],
        *[(m, f"{m}_partition_q12_m20_T300", f"{m}_m20", "t4/m20", {}) for m in "AB"],
        *[(m, f"{m}_partition_q12_m16_T600", f"{m}_m16_T600", "t4/T600", {"T4_T_K": "600"}) for m in "AB"],
        *[(m, f"{m}_partition_q12_m16_T300_fullsum", f"{m}_m16", "t4/trisym0", {}) for m in "AB"],
        *[(m, f"{m}_partition_q12_m8_T300_trisym{t}", f"{m}_m8", f"t4/m8/t{t}", {}) for m in "AB" for t in (1, 0)],
    ]
    if "t4" in steps:
        for m, rn, lt, sub, env in t4:
            rd = R / rn
            run("t4_analysis.py", m, rd / "t3_partition.txt", rd, lab[lt], D[m] / sub, env=env)
        for m in "AB":   # 8^3 file names as in the original layout (t4_{A,B}8_t{0,1}_*)
            for t in (1, 0):
                for kind in ("summary.json", "modes.csv"):
                    src = D[m] / f"t4/m8/t{t}/t4_{m}_{kind}"
                    if src.exists():
                        src.replace(D[m] / f"t4/m8/t{t}/t4_{m}8_t{t}_{kind}")
            run("t4_trisym_sensitivity.py", *(D[m] / f"t4/m8/t{t}/t4_{m}8_t{t}_summary.json" for t in (1, 0)),
                *(D[m] / f"t4/m8/t{t}/t4_{m}8_t{t}_modes.csv" for t in (1, 0)), D[m] / "t4/m8/trisym_sensitivity_{}.json".format(m))

    # ---- frequency-band analysis (T5, T5b) --------------------------------------------------------------------
    if "t5" in steps:
        for m in "AB":
            rd = R / f"{m}_partition_q12_m16_T300_bands"
            run("t5_control_analysis.py", m, rd / "t3_partition.txt", rd, t5lab[m], D[m] / f"t4/t4_{m}_modes.csv",
                D[m] / f"t5/t5_{m}.json")
            rd = R / f"{m}_partition_q12_m16_T300_lowbands"
            run("t5_control_analysis.py", m, rd / "t3_partition.txt", rd, t5blab[m], D[m] / f"t4/t4_{m}_modes.csv",
                D[m] / f"t5b/t5b_{m}.json", lab[f"{m}_m16"])

    # ---- damping-ratio shares (16^3, 300 K) ------------------------------------------------------------------
    if "kp_eta" in steps:
        ra = [f"a345={R / 'A_rta_q12_m16_T300'}", f"a198={R / 'A198_rta_q12_m16_T300'}", f"a600={R / 'A600_rta_q12_m16_T300'}"]
        for g, pre in (("avg", "kp_eta_300K_m16_gavg"), ("raw", "kp_eta_300K_m16")):
            run("kp_eta_distribution.py", f"--gamma={g}", 300, D["A"] / pre, *ra)
            run("kp_eta_distribution.py", f"--gamma={g}", 300, D["B"] / pre, f"H1={R / 'B_rta_q12_m16_T300'}")

    # ---- full sum against reduced sum (widths and kappa) -----------------------------------------------------
    if "t8" in steps:
        run("t8_widths.py", R / "A_rta_q12_m16_T300", R / "A_partition_q12_m16_T300_fullsum",
            D["A"] / "t4/trisym0/t8_widths.json", stdout=D["A"] / "reports/t8_widths_A.md")
        run("t8_widths.py", R / "B_rta_q12_m16_T300", R / "B_partition_q12_m16_T300_fullsum",
            D["B"] / "t4/trisym0/t8b_widths_B.json", stdout=D["B"] / "reports/t8_widths_B.md")
        for m in "AB":
            run("t8_widths.py", R / f"{m}_partition_q12_m8_T300_trisym1", R / f"{m}_partition_q12_m8_T300_trisym0",
                D[m] / f"t4/m8/t8_widths_{m}_m8.json")

    # ---- SCPH soft modes of model B (H1, H2) ------------------------------------------------------------------
    if "scph" in steps:
        run("scph_softmodes.py", D["B"] / "b3_scph_softmodes.csv",
            *[f"{f}_q{q}=" + str(R / ("H2" if f == "H2" else "B")) + f"_scph_q{q}"
              for f in ("H1", "H2") for q in (2, 4, 8, 12)])

    # ---- transport summaries (T1) -----------------------------------------------------------------------------
    if "t1" in steps:
        run("a3_transport_analysis.py", "H1_q12_m16", R / "B_rta_q12_m16_T300", R / "B_evec_q12_m16_T300",
            stdout=D["B"] / "t1/H1_q12_m16.json")
        run("a3_transport_analysis.py", "H1_q8_m12", R / "B_rta_q8_m12_T300", stdout=D["B"] / "t1/H1_q8_m12.json")
        for ph in ("cell", "atom"):
            for lab_, rta, ev, m in (("a345", "A_rta_q12_m16_T300", "A_evec_q12_m16_T300", "A"),
                                     ("a198", "A198_rta_q12_m16_T300", "A198_evec_q12_m16_T300", "A"),
                                     ("a600", "A600_rta_q12_m16_T300", "A600_evec_q12_m16_T300", "A"),
                                     ("H1", "B_rta_q12_m16_T300", "B_evec_q12_m16_T300", "B")):
                run("a3_transport_analysis.py", f"--phase={ph}", f"{lab_}_q12_m16", R / rta, R / ev,
                    stdout=D[m] / f"t1/charphase_{ph}_{lab_}.json")

    # ---- A next to B (300 K) -------------------------------------------------------------------------------------
    if "b4" in steps and "--a3-transport" in opts:
        harm = OUT / "_harm"
        for m in "AB":
            (harm / m).mkdir(parents=True, exist_ok=True)
            link = harm / m / f"harm_{m}_m16.mesh.evec"
            if not link.exists():
                link.symlink_to(one(R / f"{m}_harm_m16", "*.mesh.evec").resolve())
        both = D["B"] / "kp_eta_300K_m16_shares_AB.csv"   # per-branch shares of A and B in one table
        la, lb = ((D[m] / "kp_eta_300K_m16_shares.csv").read_text().splitlines() for m in "AB")
        both.write_text("\n".join(la + lb[1:]) + "\n")
        run("b4_model_comparison.py", D["B"] / "b4_model_comparison.csv", D["B"] / "t1/H1_q12_m16.json",
            D["B"] / "t1/H1_q8_m12.json", both, D["B"] / "b3_scph_softmodes.csv",
            opts["--a3-transport"], harm)

    # ---- census tables: family populations and frequency ranges (16^3, 300 and 600 K; four definitions) --------------------
    if "census_tables" in steps:
        for m in "AB":
            for pol in ("weighted", "plain"):
                for c in ("0.4", "0.5", "0.6"):
                    run("t2_census_tables.py", f"--polar={pol}", D[m] / "census/t2_census", D[m] / f"census/tab_{pol}_c{c}",
                        c, "0.5", f"{m}300", f"{m}600", stdout=D[m] / f"reports/census_{pol}_c{c}_{m}.md")

    # ---- report tables (Markdown; the SI tables are built from these) ---------------------------------------------
    if "reports" in steps:
        rp = D["B"] / "reports"
        chk = rp / "_checks.json"
        chk.write_text(json.dumps({"rows": []}))   # run checks are recorded per run in the manifest, not here
        run("t4_report.py", D["A"] / "t4/t4_A_summary.json", D["B"] / "t4/t4_B_summary.json",
            D["A"] / "t4/t4_A_modes.csv", D["B"] / "t4/t4_B_modes.csv", D["A"] / "t4/m8/trisym_sensitivity_A.json",
            D["B"] / "t4/m8/trisym_sensitivity_B.json", chk, rp / "t4_report.md", rp)
        for m in "AB":
            (rp / f"t4_fig_{m}.csv").replace(D[m] / f"t4/t4_fig_{m}.csv")
        run("t5_report.py", D["A"] / "t5/t5_A.json", D["B"] / "t5/t5_B.json", chk, rp / "t5_report.md")
        run("t5_report.py", D["A"] / "t5b/t5b_A.json", D["B"] / "t5b/t5b_B.json", chk, rp / "t5b_report.md", "T5b")
        run("t6_report.py", *(D[m] / "t4/t4_{}_summary.json".format(m) for m in "AB"),
            *(D[m] / "t4/m20/t4_{}_summary.json".format(m) for m in "AB"), rp / "t6_report.md")
        run("t7_report.py", *(D[m] / "t4/t4_{}_summary.json".format(m) for m in "AB"),
            *(D[m] / "t4/T600/t4_{}_summary.json".format(m) for m in "AB"), rp / "t7_report.md")
        for m in "AB":
            run("t8_report.py", D[m] / f"t4/t4_{m}_summary.json", D[m] / f"t4/trisym0/t4_{m}_summary.json",
                D[m] / f"reports/t8_report_{m}.md", m)

    # ---- mode counts by band (labels and frequencies only) ---------------------------------------------------------
    if "counts" in steps:
        both = D["B"] / "labels_both"   # links to the label tables of both models (the count scripts read both)
        (both / "T600").mkdir(parents=True, exist_ok=True)
        for m in "AB":
            for rel in (f"labels_{m}_m16.csv", f"T600/labels_{m}_m16_T600.csv"):
                link = both / rel
                if not link.exists():
                    link.symlink_to((D[m] / "labels" / rel).resolve())
            l5 = both / f"t5b_labels_{m}_m16.txt"
            if not l5.exists():
                l5.symlink_to(t5blab[m].resolve())
            j5 = both / f"t5b_{m}.json"
            if not j5.exists():
                j5.symlink_to((D[m] / f"t5b/t5b_{m}.json").resolve())
        run("rot_below_105.py", both, stdout=D["B"] / "reports/rot_below_105.md")
        run("v2_nonrot_below80.py", both, both, both, stdout=D["B"] / "reports/nonrot_below80.json")
        for m in "AB":   # count-only sampling test of the band ratios
            run("p3_sampling.py", m, lab[f"{m}_m16"], D[m] / f"t5/t5_{m}.json", D[m] / f"t5b/t5b_{m}.json",
                D[m] / f"p3_sampling_{m}.json")

    # ---- figures -------------------------------------------------------------------------------------------------
    if "figures" in steps:
        fd = OUT / "figures"
        run("figures/make_t4_figures.py", D["A"], D["B"], fd)
        run("figures/make_t4_figures.py", D["A"], D["B"], fd, "t5b")
        if "--source-bands" in opts:
            run("figures/make_figure_drafts.py", D["A"] / "kp_eta_300K_m16_gavg_cumulative.csv",
                D["B"] / "kp_eta_300K_m16_gavg_cumulative.csv", D["B"] / "b3_scph_softmodes.csv",
                opts["--source-bands"], fd)
    print("done:", OUT)


if __name__ == "__main__":
    main()
