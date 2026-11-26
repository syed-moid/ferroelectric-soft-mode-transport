#!/usr/bin/env python3
"""Number table of the paper (v4): every number with its source file and run, in three parts.

Format of data/processed/manuscript_numbers_v8.csv (key, value, unit, description, source, job, status). Values are
read from the processed outputs where a machine-readable source exists; the few that exist only in a fit record are
copied with status REPORT. Parts (written to separate files):
  a    model A, its refits and literature values          -> PROC/manuscript_numbers_v4_a.csv
  b    model B, and values of A and B side by side          -> PROC/manuscript_numbers_v4_b.csv
  held everything that needs the files of the 2015 study    -> HELD/manuscript_numbers_v4_held.csv (only with --held)
Usage (repository root): build_manuscript_numbers_v4.py PROC A3_TRANSPORT_CSV [--runs RUNS_ROOT] [--held HELD]
                         [--held-runs RUNS_ROOT_HELD]
  PROC: data/processed/v4 (a/ and b/ written by regenerate_v4.py); RUNS_ROOT: raw run directories by run-map name;
  HELD: the processed outputs of the steps that need the files of the 2015 study (not in this release).
"""
import json
import os
import sys
from datetime import date

import pandas as pd

_argv = sys.argv[1:]
_opt = {}
for _k in ("--runs", "--held", "--held-runs"):
    if _k in _argv:
        _i = _argv.index(_k)
        _opt[_k] = _argv[_i + 1]
        del _argv[_i:_i + 2]
PROC, A3CSV = _argv[0], _argv[1]
RUNS, HELD = _opt.get("--runs"), _opt.get("--held")
RUNS_HELD = _opt.get("--held-runs", RUNS)
PA, PB = f"{PROC}/a", f"{PROC}/b"
P = {"A": PA, "B": PB}
REPORT_A = "fit record of model A (example opt.log; not regenerated here)"
REPORT_B = "fit record of model B (ALM logs of h1_fit, h2_fit; not regenerated here)"
rows = []


def add(key, value, unit, desc, source, job="", status="FINAL"):
    rows.append({"key": key, "value": value, "unit": unit, "description": desc, "source": source, "job": job,
                 "status": status})


def r(x, n=4):
    return float(f"{x:.{n}g}")


# --- models -----------------------------------------------------------------------------------------------
add("a_A", 7.363, "bohr", "lattice parameter, model A (distributed example fit)", "sto_example template")
add("a_B", 7.35167842615, "bohr", "lattice parameter, model B (H1)", "model B SCPH template; T1 rta.in",
    "rta_mb_pm_q12_T300_m16")
add("born_B", "DFPT k6 header", "", "Born charges and dielectric tensor of model B (G-B2 decision)",
    REPORT_B, status="REPORT")
add("A_free_params", "3081 (698/2215/43/125)", "", "free anharmonic parameters, orders 3-6, model A",
    REPORT_A, status="REPORT")
add("A_nonzero_params", "530/1340/37/83", "", "nonzero parameters after LASSO, orders 3-6, model A",
    REPORT_A, status="REPORT")
add("A_residual", 2.30451, "%", "ALM training residual, model A", REPORT_A,
    status="REPORT")
add("B_ntrain", 128, "", "training configurations (DFSET_train), model B", REPORT_B,
    status="REPORT")
add("B_free_params", "698/2215/43/125 (FC2 fixed: H1; 45 free: H2)", "", "free parameters, orders 3-6",
    REPORT_B, status="REPORT")
for fit, alpha, cv, trE, trR, teE, teR in (("H1", 3.36918e-7, 2.140, 1.538, 10.55, 2.326, 15.78),
                                          ("H2", 3.48779e-7, 2.041, 1.476, 10.13, 2.180, 14.79)):
    add(f"{fit}_L1_alpha", alpha, "", f"LASSO penalty (4-fold CV minimum), {fit}", REPORT_B,
        status="REPORT")
    add(f"{fit}_cv_err", cv, "%", f"mean CV validation force error, {fit}", REPORT_B,
        status="REPORT")
    add(f"{fit}_train_EF_RMSE", f"{trE} / {trR}", "% / meV/Å", f"training force error / RMSE, {fit}",
        REPORT_B, status="REPORT")
    add(f"{fit}_test_EF_RMSE", f"{teE} / {teR}", "% / meV/Å", f"held-out test force error / RMSE, {fit}",
        REPORT_B, status="REPORT")
add("A_on_Btest_EF_RMSE", "4.19 / 28.4", "% / meV/Å", "model A evaluated on model B's test set (descriptive)",
    REPORT_B, status="REPORT")

# --- harmonic and SCPH frequencies --------------------------------------------------------------------------
cmp = pd.read_csv(f"{PB}/b4_model_comparison.csv")
for pt in ("Gamma", "R", "M"):
    row = cmp[cmp.quantity.str.startswith(f"harmonic frequencies at {pt}")].iloc[0]
    add(f"harm_{pt}_A", row.model_A, "cm-1", f"bare harmonic frequencies at {pt}, model A (negative = imaginary)",
        row.model_A_source)
    add(f"harm_{pt}_B", row.model_B_H1, "cm-1", f"bare harmonic frequencies at {pt}, model B",
        row.model_B_source)
ra = pd.read_csv("data/processed/v3/r_point_check.csv", comment="#")
ra = ra[(ra.PM == "on") & (ra.T_K == 300)].set_index("q1")
bs = pd.read_csv(f"{PB}/b3_scph_softmodes.csv")
bs = bs[bs.T_K == 300].set_index("run")
tt = pd.read_csv("data/literature/STO/tadano_2015/scph_convergence_300K.csv")
tt = tt[tt.method == "SCPH"].copy()
tt["q1"] = tt.branch.str.extract(r"_(\d+)x").astype(int)
for q in (2, 4, 8, 12):
    for pt in ("Gamma", "M", "R"):
        add(f"scph_{pt}_A_q{q}", ra.loc[q, pt], "cm-1", f"SCPH {pt} (lowest / TO1), 300 K, PM, q1 {q}^3, model A",
            "data/processed/v3/r_point_check.csv")
        for fit in ("H1", "H2"):
            add(f"scph_{pt}_{fit}_q{q}", r(bs.loc[f"{fit}_q{q}", pt], 6), "cm-1",
                f"SCPH {pt}, 300 K, PM, q1 {q}^3, {fit}", f"{PB}/b3_scph_softmodes.csv",
                f"scph_mb{'h2' if fit == 'H2' else ''}_pm_q{q}")
        v = tt[(tt.q_label == pt) & (tt.q1 == q)].omega_meV.iloc[0]
        add(f"scph_{pt}_2015_q{q}", round(v * 8.065544, 1), "cm-1", f"2015 Table I, {pt}, q1 {q}^3 ({v} meV)",
            "data/literature/STO/tadano_2015/scph_convergence_300K.csv")

# --- transport --------------------------------------------------------------------------------------------
a3 = pd.read_csv(A3CSV).set_index("label")
b16 = json.load(open(f"{PB}/t1/H1_q12_m16.json"))
b12 = json.load(open(f"{PB}/t1/H1_q8_m12.json"))
for lab, job in (("a345", "rta_pm_q12_T300_m16"), ("a198", "rta_rf_a198_pm_q12_T300_m16"),
                 ("a600", "rta_rf_a600_pm_q12_T300_m16")):
    x = a3.loc[f"{lab}_q12_m16"]
    add(f"kP_{lab}_m16", x.kappa_P_solver, "W/mK", f"kappa_P, 300 K, 16^3, q1 12^3, fit {lab}",
        "benchmark table A3_transport.csv (input)", job)
    add(f"kC_{lab}_m16", x.kappa_C_solver, "W/mK", f"kappa_C, 300 K, 16^3, fit {lab}",
        "benchmark table A3_transport.csv (input)", job)
    y = a3.loc[f"{lab}_q8_m12"]
    add(f"kP_{lab}_q8_m12", y.kappa_P_solver, "W/mK", f"kappa_P, 300 K, 12^3, q1 8^3, fit {lab}",
        "benchmark table A3_transport.csv (input)")
    add(f"kC_{lab}_q8_m12", y.kappa_C_solver, "W/mK", f"kappa_C, 300 K, 12^3, q1 8^3, fit {lab}",
        "benchmark table A3_transport.csv (input)")
add("kP_B_m16", b16["kappa_P_solver"], "W/mK", "kappa_P, 300 K, 16^3, q1 12^3, model B",
    f"{PB}/t1/H1_q12_m16.json", "rta_mb_pm_q12_T300_m16")
add("kC_B_m16", b16["kappa_C_solver"], "W/mK", "kappa_C, 300 K, 16^3, model B", f"{PB}/t1/H1_q12_m16.json",
    "rta_mb_pm_q12_T300_m16")
add("kP_B_q8_m12", b12["kappa_P_solver"], "W/mK", "kappa_P, 300 K, 12^3, q1 8^3, model B",
    f"{PB}/t1/H1_q8_m12.json", "rta_mb_pm_q8_T300_m12")
add("kC_B_q8_m12", b12["kappa_C_solver"], "W/mK", "kappa_C, 300 K, 12^3, q1 8^3, model B",
    f"{PB}/t1/H1_q8_m12.json", "rta_mb_pm_q8_T300_m12")
add("kP_2015", 8.7, "W/mK", "2015 study, 300 K, 12^3 BTE (digitized Fig. 7)",
    "data/literature/STO/tadano_2015/kappa_total_T.csv")
km = pd.read_csv("data/processed/sto_production/kappa_mesh_T.csv", comment="#")
km = km[km.T_K == 300]
for _, x in km.iterrows():
    add(f"kP_A_mesh{int(x.mesh)}", x.kappa_P, "W/mK", f"kappa_P, 300 K, {int(x.mesh)}^3, model A",
        "data/processed/sto_production/kappa_mesh_T.csv", x.job)
    add(f"kC_A_mesh{int(x.mesh)}", x.kappa_C, "W/mK", f"kappa_C, 300 K, {int(x.mesh)}^3, model A",
        "data/processed/sto_production/kappa_mesh_T.csv", x.job)

# --- eta (group-averaged widths, the analysis plan; per-branch values kept as secondary) ---------------------------
g = pd.concat([pd.read_csv(f"{p}/kp_eta_300K_m16_gavg_shares.csv") for p in (PA, PB)]).set_index("run")
w = pd.concat([pd.read_csv(f"{p}/kp_eta_300K_m16_shares.csv") for p in (PA, PB)]).set_index("run")
jobs = {"a345": "rta_pm_q12_T300_m16", "a198": "rta_rf_a198_pm_q12_T300_m16",
        "a600": "rta_rf_a600_pm_q12_T300_m16", "H1": "rta_mb_pm_q12_T300_m16"}
for run in g.index:
    nm = "B" if run == "H1" else run
    for t in (0.02, 0.05, 0.1):
        add(f"share_eta_gt_{t}_{nm}", r(100 * g.loc[run, f"share_kP_eta_gt_{t}"]), "%",
            f"share of kappa_P from modes with eta > {t}, 300 K, 16^3, group-averaged widths",
            f"{P['B' if run == 'H1' else 'A']}/kp_eta_300K_m16_gavg_shares.csv", jobs[run])
        add(f"share_eta_gt_{t}_{nm}_perbranch", r(100 * w.loc[run, f"share_kP_eta_gt_{t}"]), "%",
            "same with per-branch widths (released convention)", f"{P['B' if run == 'H1' else 'A']}/kp_eta_300K_m16_shares.csv",
            jobs[run], "SECONDARY")
    add(f"eta_max_{nm}", r(g.loc[run, "eta_max"]), "",
        f"largest group-averaged eta on 16^3, 300 K (R triplet, {g.loc[run, 'eta_max_omega_cm1']} cm-1)",
        f"{P['B' if run == 'H1' else 'A']}/kp_eta_300K_m16_gavg_shares.csv", jobs[run])
    add(f"eta_max_{nm}_perbranch", r(w.loc[run, "eta_max"]), "", "same, per-branch widths (released convention)",
        f"{P['B' if run == 'H1' else 'A']}/kp_eta_300K_m16_shares.csv", jobs[run], "SECONDARY")
    add(f"R_omega_{nm}_m16", g.loc[run, "eta_max_omega_cm1"], "cm-1", "R triplet frequency in the 16^3 run",
        f"{P['B' if run == 'H1' else 'A']}/kp_eta_300K_m16_gavg_shares.csv", jobs[run])

# --- conclusion 2 under the existing rules ------------------------------------------------------------------
for run, nm in (("a345", "A"), ("H1", "B"), ("a198", "a198"), ("a600", "a600")):
    for ph in ("cell", "atom"):
        _pt = f"{P['B' if run == 'H1' else 'A']}/t1/charphase_{ph}_{run}.json"
        d = json.load(open(_pt))
        if ph == "cell":
            add(f"hh_freq_{nm}", r(100 * d["share_kC_hardhard_freq"]), "%",
                "hard/hard share of kappa_C, frequency rule, 300 K, 16^3", _pt,
                jobs[run])
        add(f"hh_char_{ph}_{nm}", r(100 * d["share_kC_hardhard_char"]), "%",
            "hard/hard share of kappa_C, character rule, " +
            ("published (cell) convention" if ph == "cell" else "star-invariant (atom-position) convention"),
            _pt, jobs[run], "FINAL" if ph == "atom" else "SECONDARY")

# --- frozen family census (characters and frequencies only) --------------------------------------------------
for m in ("A", "B"):
    lab = pd.read_csv(f"{P[m]}/labels/labels_{m}_m16.csv")
    for c in (0.4, 0.5, 0.6):
        vc = lab[f"primary_c{c}"].value_counts()
        add(f"census_{m}_primary_c{c}", f"pol {vc.get('pol', 0)} / rot {vc.get('rot', 0)} / mixed "
            f"{vc.get('mixed', 0)}", "modes", f"star-weighted family counts, 16^3, primary definition, c_min {c}",
            f"{P[m]}/labels/labels_{m}_m16.csv")

# --- T4 / T6 / T7: who carries the heat and who scatters the carriers -------------------------
def t4_rows(path, tag, desc, job, status="FINAL"):
    if not os.path.exists(path):
        return
    s = json.load(open(path))
    d = s["definitions"]["primary_c0.5"]
    pre = f"{desc}; primary c_min 0.5; carriers = frozen 50 % rule"
    add(f"kP_modesum_{tag}", r(s["kappa_P_modesum"]), "W/mK", f"kappa_P (mode sum = solver), {desc}", path, job, status)
    add(f"carriers_{tag}", f"{s['carriers']} / {s['carrier_groups']}", "modes / groups", f"carriers, {desc}", path, job,
        status)
    for k, nm in (("0", "no"), ("1", "one"), ("2", "both")):
        add(f"nsoft{k}_{tag}", r(100 * d["kappa_weighted_nsoft_frac"][k]), "%",
            f"carriers' linewidth with {nm} soft partner(s), kappa_P-weighted, {pre}", path, job, status)
    add(f"pair_neither_rot_{tag}", r(100 * d["kappa_weighted_pair_frac"]["neither/rot"]), "%",
        f"carriers' linewidth from partner pair neither/rot, {pre}", path, job, status)
    for f in ("neither", "pol", "rot"):
        add(f"kP_share_{f}_{tag}", r(100 * d["kappa_P_share_by_family"][f]), "%",
            f"share of kappa_P carried by family {f} (all modes), {pre}", path, job, status)
    add(f"kP_share_eta01_{tag}", r(100 * s["share_kP_eta_gt_0.1"], 6), "%",
        f"share of kappa_P carried by modes with group-averaged eta > 0.1, {desc}", path, job, status)


for m in ("A", "B"):
    t4_rows(f"{P[m]}/t4/t4_{m}_summary.json", f"{m}_m16", f"model {m}, 300 K, 16^3, TRISYM 1",
            "T4")
    t4_rows(f"{P[m]}/t4/m20/t4_{m}_summary.json", f"{m}_m20", f"model {m}, 300 K, 20^3, TRISYM 1",
            "T6")
    t4_rows(f"{P[m]}/t4/T600/t4_{m}_summary.json", f"{m}_m16_T600", f"model {m}, 600 K, 16^3, TRISYM 1"
            + ("; outside the training distribution" if m == "B" else ""), "T7",
            "FINAL" if m == "A" else "OUTSIDE_TRAINING")
t4_rows(f"{PA}/t4/trisym0/t4_A_summary.json", "A_m16_trisym0", "model A, 300 K, 16^3, TRISYM 0 (full sum)",
        "T8")
_w = f"{PA}/t4/trisym0/t8_widths.json"
if os.path.exists(_w):
    w = json.load(open(_w))
    for key, nm in (("group_averaged", "group-averaged"), ("per_branch", "per-branch")):
        for st, lab in (("median_abs", "median"), ("q95_abs", "95% quantile"), ("max_abs", "maximum")):
            add(f"trisym0_m16_A_width_{key}_{st}", r(100 * w[key][st], 3), "%",
                f"{lab} of |relative difference| of {nm} widths, TRISYM 0 vs 1, model A, 300 K, 16^3", _w,
                "T8")
    k = w["kappa_P_modesum"]
    add("trisym0_m16_A_kP_rel", r(100 * (k[1] / k[0] - 1), 4), "%", "relative change of kappa_P, TRISYM 0 vs 1, model A, "
        "300 K, 16^3", _w, "T8")

# --- T5 / T5b frequency-matched controls (per-mode share, % per 1000 modes) ------------------------------------
for ctl, band_names, job in (("t5", ("<105", "105-145", "145-160"), "T5"),
                             ("t5b", ("<80", "80-95", "95-105", ">=105"), "T5b")):
    for m in ("A", "B"):
        path = f"{P[m]}/{ctl}/{ctl}_{m}.json"
        if not os.path.exists(path):
            continue
        bd = json.load(open(path))["bands"]
        for b in band_names:
            e = bd[b]
            if e["rot"]["per_mode"]["total"] is None:
                continue
            add(f"{ctl}_{m}_{b}_s_rot", r(1e5 * e["rot"]["per_mode"]["total"]), "% per 1000 modes",
                f"carriers' linewidth per rot partner mode, band {b} cm-1, model {m}, 300 K, 16^3", path, job)
            add(f"{ctl}_{m}_{b}_s_nonrot", r(1e5 * e["nonrot"]["per_mode"]["total"]), "% per 1000 modes",
                f"carriers' linewidth per non-rot partner mode, band {b} cm-1, model {m}", path, job)
            for part, nm in (("total", "total"), ("decay_of_decay", "decay"), ("absorption_of_absorption", "absorption")):
                add(f"{ctl}_{m}_{b}_ratio_{nm}", r(e["ratio_rot_over_nonrot"][part], 3), "",
                    f"s_rot / s_nonrot ({nm}), band {b} cm-1, model {m}", path, job)
            for nm in ("rot", "nonrot"):
                if e[nm].get("omega_cm1"):
                    add(f"{ctl}_{m}_{b}_omega_mean_{nm}", r(e[nm]["omega_cm1"]["mean"], 5), "cm-1",
                        f"mean frequency of the {nm} partner modes, band {b} cm-1, model {m}, full 16^3 mesh", path, job)

# --- the analysis plan: numbers still missing from Draft 1 (existing outputs only) ------------------------------------
_kc = f"{RUNS}/B_rta_q12_m20_T300/rta_mb_pm_q12_T300_m20.kl_coherent" if RUNS else ""
if _kc and os.path.exists(_kc):
    _v = [ln.split() for ln in open(_kc) if ln.strip() and not ln.startswith("#")][0]
    add("kC_B_m20", _v[1], "W/mK", "kappa_C (xx), model B, 300 K, 20^3 (base run bB; patched tB identical)", _kc,
        "T6")
for m, t, path in (("A", 300, f"{PA}/labels/labels_A_m16.csv"), ("B", 300, f"{PB}/labels/labels_B_m16.csv"),
                   ("A", 600, f"{PA}/labels/T600/labels_A_m16_T600.csv"), ("B", 600, f"{PB}/labels/T600/labels_B_m16_T600.csv")):
    L = pd.read_csv(path)
    rot = L[L["primary_c0.5"] == "rot"]
    st = "OUTSIDE_TRAINING" if (m, t) == ("B", 600) else "FINAL"
    job = "scripts/v4/rot_below_105.py (labels and frequencies only)"
    add(f"rot_below105_{m}_T{t}", str(int((rot.Omega_cm1 < 105).sum())), "modes",
        f"rot-family modes (primary c_min 0.5) below 105 cm-1, model {m}, {t} K, full 16^3 mesh (61440 modes)",
        path, job, st)
    add(f"rot_range_{m}_T{t}", f"{rot.Omega_cm1.min():.2f}–{rot.Omega_cm1.max():.2f}", "cm-1",
        f"frequency range of the rot family (primary c_min 0.5), model {m}, {t} K, 16^3", path, job, st)
_exp = "data/processed/v3/exp_300K_crossread.csv"
if os.path.exists(_exp):
    E = pd.read_csv(_exp, comment="#")
    for _, e in E.iterrows():
        add(f"exp_kappa_{e.dataset}_{e.figure.split()[0].replace('&', '')}_{e.T_K:g}K", f"{e.kappa_W_mK} ± {e.unc_W_mK}",
            "W/mK", f"measured total kappa, {e.dataset}, at {e.T_K:g} K, read from {e.figure} ({e.note}); "
            + ("used in the v3 paper (Fig. 3; Muta 300 K: 11.0–11.1 in the text)" if e.dataset in ("Muta2005", "Popuri2014")
               else "cross-read in the same file; not cited in the v3 paper"), _exp,
            "earlier release (v3)", "READ_FROM_FIGURE")

# --- the analysis plan: archived 2015 fit (P1, P1b, P2; rows appear when their analysis exists) ---------------------------
_CASE = {"b": "distributed fit, NONANALYTIC 2, BORN", "c": "archived 2015 fit, NONANALYTIC 3, BORN",
         "d": "archived 2015 fit, NONANALYTIC 2, BORN", "e": "archived 2015 fit, NONANALYTIC 2, Born data of the 2015 study"}
for _f, _job in (("p1_softmodes.csv", "P1"), ("p1b_softmodes.csv", "P1b")):
    _f = f"{HELD}/{_f}" if HELD else ""
    if _f and os.path.exists(_f):
        PP = pd.read_csv(_f)
        for run, g in PP.groupby("run", sort=False):
            c, q = run.split("_q")
            g300 = g[g.T_K == 300].iloc[0]
            for pt in ("Gamma", "M", "R"):
                add(f"p15_{run}_{pt}", r(g300[pt], 6), "cm-1",
                    f"SCPH {pt} ({'lowest optical' if pt == 'Gamma' else 'lowest'}) at 300 K, case {c} ({_CASE[c]}), q1 {q}^3",
                    _f, _job)
            add(f"p15_{run}_lowest_conv_T", f"{g300.lowest_contiguous_converged_T:g}", "K",
                f"lowest contiguous converged SCPH temperature, case {c}, q1 {q}^3", _f, _job)
_p2 = f"{HELD}/p2_analysis.json" if HELD else ""
if _p2 and os.path.exists(_p2):
    P2 = json.load(open(_p2))
    add("p15_kappa_archived_300K", r(P2["archived_kappa_300K"], 7), "W/mK",
        "kappa at 300 K of the 2015 study (numerical data of its Fig. 7, provided by its author)",
        "2015 kappa(T) curve (input of the held part)", "P2", "COMMUNICATED")
    for k, v in P2["runs"].items():
        if k == "a":
            continue
        add(f"p15_kP_{k}_q8_m12", r(v["kappa_P_xx"], 5), "W/mK",
            f"kappa_P (xx) at 300 K, step {k}-kappa ({_CASE[k]}), q1 8^3, transport 12^3", _p2,
            "P2")
        add(f"p15_kC_{k}_q8_m12", r(v["kappa_C_xx"], 5), "W/mK",
            f"kappa_C (xx) at 300 K, step {k}-kappa, q1 8^3, transport 12^3", _p2,
            "P2")
        add(f"p15_kP_{k}_minus_archived_rel", r(100 * v["minus_archived_rel"], 4), "%",
            f"(kappa_P of step {k}-kappa - archived 8.756089) / archived", _p2,
            "P2")

# --- the analysis plan: P3 numbers (archived 2015 fit under the settings of model A; the P3 report tables (Markdown)) ---------------
_t4p = f"{HELD}/t4/t4_P_summary.json"
if HELD and RUNS_HELD and all(os.path.exists(f"{HELD}/{x}") for x in (
        "t4/t4_P_summary.json", "t5_P.json", "t5b_P.json", "p3_kp_eta_300K_m16_gavg_shares.csv", "p3_scph_q12_softmodes.csv")):
    _job = "P3"
    _d = "archived 2015 fit, settings of model A (NONANALYTIC 3, BORN, q1 12^3), 300 K, 16^3"
    T4P = json.load(open(_t4p))
    _pr = T4P["definitions"]["primary_c0.5"]
    _src4 = _t4p
    add("p15_p3_kP_rot_share", r(100 * _pr["kappa_P_share_by_family"]["rot"], 6), "%",
        f"share of kappa_P carried by the rotational family, primary c0.5, {_d}", _src4, _job)
    add("p15_p3_kP_pol_share", r(100 * _pr["kappa_P_share_by_family"]["pol"], 6), "%",
        f"share of kappa_P carried by the polar family, primary c0.5, {_d}", _src4, _job)
    add("p15_p3_lw_nosoft", r(100 * _pr["kappa_weighted_nsoft_frac"]["0"], 6), "%",
        f"carrier-weighted linewidth share with no soft partner, primary c0.5, {_d}", _src4, _job)
    add("p15_p3_lw_neither_rot", r(100 * _pr["kappa_weighted_pair_frac"]["neither/rot"], 6), "%",
        f"carrier-weighted linewidth share, partner pair neither/rot, primary c0.5, {_d}", _src4, _job)
    add("p15_p3_carriers", T4P["carriers"], "modes", f"carriers (irreducible modes, frozen 50 % rule), {_d}", _src4, _job)
    add("p15_p3_carrier_groups", T4P["carrier_groups"], "groups", f"carrier degenerate groups, {_d}", _src4, _job)
    for _f, _bands in ((f"{HELD}/t5_P.json", ("<105", "105-145", "145-160")),
                       (f"{HELD}/t5b_P.json", ("<80", "80-95", "95-105"))):
        _B = json.load(open(_f))["bands"]
        for b in _bands:
            k = b.replace("<", "lt").replace("-", "_")
            add(f"p15_p3_ratio_{k}", r(_B[b]["ratio_rot_over_nonrot"]["total"], 6), "",
                f"per-mode ratio rot / non-rot partners, band {b} cm-1 (total width), {_d}", _f, _job)
            add(f"p15_p3_N_rot_{k}", _B[b]["rot"]["N_modes"], "modes", f"rotational modes in band {b} cm-1 (16^3 mesh), {_d}",
                _f, _job)
            add(f"p15_p3_N_nonrot_{k}", _B[b]["nonrot"]["N_modes"], "modes",
                f"non-rotational modes in band {b} cm-1 (16^3 mesh), {_d}", _f, _job)
    _rd = f"{RUNS_HELD}/P_rta_q12_m16_T300_settingsA"
    for ext, nm in (("kl", "kP"), ("kl_coherent", "kC")):
        _row = [l for l in open(f"{_rd}/p15_rta_q12_T300_m16.{ext}") if not l.startswith("#") and l.strip()][0].split()
        add(f"p15_p3_{nm}_m16", float(_row[1]), "W/mK", f"{nm} (xx) at 300 K, {_d}",
            f"P_rta_q12_m16_T300_settingsA/p15_rta_q12_T300_m16.{ext}", "p15_rta_q12_T300_m16")
    _E = pd.read_csv(f"{HELD}/p3_kp_eta_300K_m16_gavg_shares.csv").set_index("run").loc["P"]
    for t in ("0.02", "0.05", "0.1"):
        add(f"p15_p3_kP_eta_gt_{t}", r(100 * _E[f"share_kP_eta_gt_{t}"], 6), "%",
            f"share of kappa_P from modes with group-averaged eta > {t}, {_d}",
            f"{HELD}/p3_kp_eta_300K_m16_gavg_shares.csv", _job)
    add("p15_p3_eta_max", r(_E.eta_max, 6), "", f"largest group-averaged eta (R triplet, {_E.eta_max_omega_cm1:.2f} cm-1), {_d}",
        f"{HELD}/p3_kp_eta_300K_m16_gavg_shares.csv", _job)
    _S = pd.read_csv(f"{HELD}/p3_scph_q12_softmodes.csv")
    _S = _S[(_S.run == "arch_q12") & (_S.T_K == 300)].iloc[0]
    for pt in ("Gamma", "X", "M", "R"):
        add(f"p15_p3_q12_{pt}", r(_S[pt], 6), "cm-1",
            f"SCPH {pt} ({'lowest optical' if pt == 'Gamma' else 'lowest'}) at 300 K, archived 2015 fit, NONANALYTIC 3, BORN, q1 12^3",
            f"{HELD}/p3_scph_q12_softmodes.csv", "p15_scph_q12")
    add("p15_p3_q12_lowest_conv_T", f"{_S.lowest_contiguous_converged_T:g}", "K",
        "lowest contiguous converged SCPH temperature, archived 2015 fit, NONANALYTIC 3, BORN, q1 12^3",
        f"{HELD}/p3_scph_q12_softmodes.csv", "p15_scph_q12")

# --- measured R frequency at 300 K (literature; open circles of Fig. 1(d) of Fauqué et al. 2022) -------------------
_fq = "data/literature/STO/fauque_2022/dho_fits_T.csv"
_F = pd.read_csv(_fq)
_F = _F[(_F.T_K == 300) & (_F.q_label == "R") & (_F.method == "experimental")]
assert len(_F) == 1, _F
add("lit_R_omega_300K_meV", float(_F.omega_meV.iloc[0]), "meV",
    "measured R-point soft-mode energy at 300 K, Q = (1.5, 0.5, 0.5) (Fauqué et al. 2022, Fig. 1(d), open circles; "
    "reading uncertainty +-0.62 meV, see the README of the data folder)", _fq)
add("lit_R_omega_300K_cm1", round(float(_F.omega_meV.iloc[0]) * 8.065544, 2), "cm-1",
    "the same in cm-1 (x 8.065544)", _fq)

# --- the 300 K lifetime test (P6): our velocities and heat capacities with the 2015 per-mode lifetimes (2015 files) ---
_p6 = f"{HELD}/p6/p6.json"
if HELD and os.path.exists(_p6):
    _job, P6 = "P6", json.load(open(_p6))
    _d = "300 K, configuration of the 2015 study, q1 8^3, transport 12^3, present group velocities"
    _arch = P6["archived"]
    for st, what in (("step0", "present lifetimes and heat capacities"),
                     ("step1", "2015 lifetimes, present heat capacities"),
                     ("step2", "2015 lifetimes, heat capacities at the 2015 frequencies")):
        v = P6["all"][st]
        add(f"p15_p6_kP_{st}", r(v, 6), "W/mK", f"kappa_P mode sum, {what}; {_d}", _p6, _job)
        add(f"p15_p6_kP_{st}_minus_archived", r(v - _arch, 6), "W/mK",
            f"kappa_P mode sum ({what}) minus the archived 300 K value {_arch}", _p6, _job)
        add(f"p15_p6_kP_{st}_minus_archived_rel", r(100 * (v / _arch - 1), 4), "%",
            f"(kappa_P mode sum ({what}) - archived) / archived", _p6, _job)
    add("p15_p6_step1_minus_step0", r(P6["all"]["step1"] - P6["all"]["step0"], 6), "W/mK",
        f"change of kappa_P on substituting the 2015 lifetimes; {_d}", _p6, _job)
    add("p15_p6_step2_minus_step1", r(P6["all"]["step2"] - P6["all"]["step1"], 6), "W/mK",
        f"change of the substituted sum on taking heat capacities at the 2015 frequencies; {_d}", _p6, _job)
    add("p15_p6_n_modes", P6["n_modes"], "modes", "irreducible modes in the sums (acoustic modes at Gamma excluded)",
        _p6, _job)
    _M = pd.read_csv(f"{HELD}/p6/p6_modes.csv")
    add("p15_p6_n_tau_2015_longer", int((_M.tau_file_ps > _M.tau_ours_ps).sum()), "modes",
        "modes whose 2015 lifetime is longer than the present one (file-order pairing)", f"{HELD}/p6/p6_modes.csv", _job)
    add("p15_p6_sorted_pairing_modes_changed", P6["sorted_pairing_modes_changed"], "modes",
        "modes whose 2015 partner changes under frequency-sorted pairing within each wave-vector group", _p6, _job)
    add("p15_p6_step1_sorted_minus_positional", r(P6["all"]["step1_sorted"] - P6["all"]["step1"], 6), "W/mK",
        "substituted sum, frequency-sorted minus file-order pairing", _p6, _job)
    for sub, what in (("without_point25", "without the wave-vector group of unmatched degeneracy (15 modes)"),
                      ("without_p25_b10_11", "without the two unmatched branches of that group")):
        add(f"p15_p6_step1_minus_step0_{sub}", r(P6[sub]["step1"] - P6[sub]["step0"], 6), "W/mK",
            f"change on substituting the 2015 lifetimes, {what}", _p6, _job)

df = pd.DataFrame(rows)
assert df.key.is_unique


def part(r):
    src, key = str(r.source), r.key
    if (HELD and src.startswith(HELD)) or key.startswith("p15_") or "held part" in src:
        return "held"
    if src.startswith(PB) or "model B" in src or key.startswith(("B_", "H1_", "H2_")) or "_B_" in f"_{key}_" \
            or "_H1_" in f"_{key}_" or "_H2_" in f"_{key}_" or key.endswith("_B_perbranch"):
        return "b"
    return "a"


df["part"] = [part(r) for r in df.itertuples()]
# sources relative to the repository (processed outputs), or named by role (held outputs, raw run directories)
for _p, _n in ((PROC, "data/processed/v4"), (HELD, "<held>"), (RUNS_HELD, "<runs>"), (RUNS, "<runs>")):
    if _p:
        df["source"] = df.source.astype(str).str.replace(str(_p).rstrip("/"), _n, regex=False)
df["source"] = df.source.str.replace(str(A3CSV), "data/processed/v4/a/benchmark/A3_transport.csv", regex=False)
for prt, dest in (("a", f"{PROC}/manuscript_numbers_v4_a.csv"), ("b", f"{PROC}/manuscript_numbers_v4_b.csv"),
                  ("held", f"{HELD}/manuscript_numbers_v4_held.csv" if HELD else None)):
    sub = df[df.part == prt].drop(columns="part")
    if dest is None or sub.empty:
        continue
    with open(dest, "w") as f:
        f.write(f"# script: scripts/v4/build_manuscript_numbers_v4.py; part {prt}\n# generated: {date.today()}\n")
        sub.to_csv(f, index=False)
    print(prt, len(sub), "rows ->", dest)
