"""Rebuild the v8 placeholder values (DRAFT_v8_SKELETON placeholder map) from production data only
(own SCPH q1 = 12^3, PM; main-text range 300-1000 K; 250 K is an SI boundary point).
Output: data/processed/manuscript_numbers_v8.csv (placeholder, value, units, description, source, status).
status: FINAL (computed here), CHOICE (a choice is pending; the value shows the default), PENDING (input missing).
The v7.3 index data/processed/manuscript_numbers.csv is left untouched."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PROD = ROOT / "data/processed/sto_production"
PJ = ROOT / "data/processed/revision"
rd = lambda f: pd.read_csv(f, comment="#")
man = json.loads((ROOT / "data/provenance/manifest_final.json").read_text())["jobs"]
km = rd(PROD / "kappa_mesh_T.csv")
eta = rd(PROD / "eta_validity_T.csv")
sec = rd(PROD / "kappa_sectors_T.csv")
coh = rd(PROD / "coherence_blocks_T.csv")
d4 = rd(PJ / "d4_bare_framings.csv")
rpc = rd(PJ / "r_point_check.csv")
vogt = rd(ROOT / "data/literature/STO/vogt_1995/omega_T.csv")
MAIN = list(range(300, 1001, 50))
rows = []


def add(ph, v, units, desc, src, status="FINAL"):
    rows.append({"placeholder": ph, "value": v, "units": units, "description": desc, "source": src, "status": status})


def k(T, m, col):
    return float(km[(km.T_K == T) & (km.mesh == m)][col].iloc[0])


add("T_min", 300, "K", "lowest main-text temperature (q1 = 12^3 meets the D0 rule against 16^3 from 300 K)",
    "n3_q1_16_check.csv")
add("q1_prod", "12x12x12", "", "SCPH q1 grid, production", "docs/recorded_rules.md (R1)")
add("scph_tol", "1e-10", "", "TOL_SCPH (MAXITER 500, MIXALPHA 0.2)", "docs/recorded_rules.md (R1) (scph.log)")
add("a_bohr", 7.363, "bohr", "fixed cubic lattice parameter (3.8963 A)", "rta.in")
add("kP_300_16_prod", round(k(300, 16, "kappa_P"), 3), "W/mK", "kappa_P(300 K), 16^3", "kappa_mesh_T.csv")
big = km[(km.T_K == 300) & (km.mesh >= 15)]
kP_best, kC_best = big.kappa_P.mean(), big.kappa_C.mean()
add("kP_300_best", f"{kP_best:.2f} (observed spread {km[(km.T_K == 300) & (km.mesh >= 12)].kappa_P.min():.2f}-"
    f"{km[(km.T_K == 300) & (km.mesh >= 12)].kappa_P.max():.2f} over 12^3-{km[km.T_K == 300].mesh.max()}^3)", "W/mK",
    f"kappa_P(300 K) best estimate = mean of the >= 15^3 meshes ({','.join(str(m) for m in sorted(big.mesh))})", "kappa_mesh_T.csv")
add("ktot_300_best", round(kP_best + kC_best, 2), "W/mK", "kappa_P + kappa_C, means of the >= 15^3 meshes at 300 K",
    "kappa_mesh_T.csv")
add("exp_gap_pct", round(100 * ((kP_best + kC_best) / 11.0 - 1), 1), "%",
    "kappa_total(300 K) vs 11 +- 1 W/mK (Tachibana 2008, single crystal, approx.)",
    "data/literature/STO/tachibana_2008/kappa_total_T.csv", "CHOICE")
add("dkP_sampling_pct", "16 (without PM: 10.31 vs 8.92 W/mK); 18 (with PM: 10.93 vs 9.24 W/mK)", "%",
    "kappa_P(300 K, 12^3 mesh) raised by SCPH q1 2^3 relative to the dense-grid value at matched settings "
    "(base = q1 8^3; +15.6 / +18.3)",
    "runs rta_on_k2_m12, rta_on_k8_m12, rta_off_k8_m12; AUDIT_REPORT_W1 (10.31)")
add("dkP_PM_pct", round(100 * (9.2384 / 8.9195 - 1), 1), "%", "PM on vs off at q1 = 8^3, 12^3 mesh, 300 K",
    "runs rta_on_k8_m12, rta_off_k8_m12")
add("dkP_q1_8_12_pct", round(100 * (9.1112 / 9.2384 - 1), 1), "%", "kappa_P(300 K, 12^3 mesh), SCPH q1 8^3 -> 12^3",
    "docs/recorded_rules.md (R1)")
add("dkP_q1_12_16_pct", round(100 * (man["rta_pm_q16_T300_m12"]["metrics"]["kappa_P"] /
                                   man["rta_pm_q12_T300_m12"]["metrics"]["kappa_P"] - 1), 1), "%",
    "kappa_P(300 K, 12^3 mesh), SCPH q1 12^3 -> 16^3", "n3_q1_16_check.csv")
s300 = km[(km.T_K == 300) & (km.mesh >= 12)].kappa_P
add("kP_mesh_spread_pct", round(100 * (s300.max() - s300.min()) / k(300, 16, "kappa_P"), 1), "%",
    "observed kappa_P(300 K) spread over 12^3-20^3, relative to 16^3", "kappa_mesh_T.csv")
r3 = rpc[(rpc.PM == "on") & (rpc.T_K == 300)].set_index("q1").R
add("R_2cubed", round(float(r3[2]), 1), "cm-1", "Omega_R(300 K), q1 = 2^3, PM", "r_point_check.csv")
add("R_8cubed", round(float(r3[8]), 1), "cm-1", "Omega_R(300 K), q1 = 8^3, PM", "r_point_check.csv")
add("dkP_12_16_pct", round(100 * (k(300, 16, "kappa_P") / k(300, 12, "kappa_P") - 1), 1), "%",
    "kappa_P(300 K) 12^3 -> 16^3", "kappa_mesh_T.csv")
seq = km[km.T_K == 300].sort_values("mesh")
add("kC_seq_300", "; ".join(f"{m}^3: {v:.3f}" for m, v in zip(seq.mesh, seq.kappa_C) if m >= 12), "W/mK",
    "kappa_C(300 K) by mesh", "kappa_mesh_T.csv")
big5 = km[(km.T_K == 300) & (km.mesh >= 15)]
add("kC_300_range", f"{big5.kappa_C.min():.2f}-{big5.kappa_C.max():.2f} (mean {big5.kappa_C.mean():.3f}, "
    f"{len(big5)} meshes: {', '.join(f'{m}^3' for m in sorted(big5.mesh))})", "W/mK",
    "kappa_C(300 K) observed range over the >= 15^3 meshes", "kappa_mesh_T.csv")
tot5 = big5.kappa_P + big5.kappa_C
add("ktot_300_range", f"{tot5.min():.2f}-{tot5.max():.2f}", "W/mK", "kappa_total(300 K) observed spread, >= 15^3 meshes",
    "kappa_mesh_T.csv")
add("dkC_20_24_pct", round(100 * (k(300, 24, "kappa_C") / k(300, 20, "kappa_C") - 1), 1), "%",
    "300 K; pre-registered prediction failed (|D| >= 13.1 %)", "kappa_mesh_T.csv; docs/recorded_rules.md (R2)")
d2s = rd(PJ / "d2_summary.csv")
rr = lambda T, m: float(d2s[(d2s.T_K == T) & (d2s.mesh == m)].kw_median_r.iloc[0])
add("r_criterion", f"r <= 0.6: 13-23 % (either sign); r = {rr(500, 16):.2f}: "
    f"{100 * abs(k(500, 20, 'kappa_C') / k(500, 16, 'kappa_C') - 1):.1f} %; r >= {rr(650, 16):.2f}: <= "
    f"{max(100 * abs(k(T, 20, 'kappa_C') / k(T, 16, 'kappa_C') - 1) for T in (650, 800, 1000)):.1f} % (650, 800, 1000 K)", "",
    "|kappa_C change| between successive meshes vs kappa_C-weighted median r of the coarser mesh",
    "d2_summary.csv; kappa_mesh_T.csv")
add("kC_settled_T", f"650, 800 and 1000 K (tested on 20^3); 900-950 K by the criterion (r = {rr(900, 16):.2f}-"
    f"{rr(950, 16):.2f} at 16^3)", "K", "kappa_C settled",
    "d2_summary.csv; docs/recorded_rules.md (R8) (1000 K 20^3 outcome)")
add("r16_900_1000", "; ".join(f"{T} K: {rr(T, 16):.2f}" for T in (900, 950, 1000)), "",
    "r(16^3) at 900-1000 K (1000 K tested on 20^3; 900/950 K by the criterion)", "d2_summary.csv")
_c = lambda T, m: float(d2s[(d2s.T_K == T) & (d2s.mesh == m)].kappa_C.iloc[0])
add("dkC_16_20_1000_pct", round(100 * (_c(1000, 20) / _c(1000, 16) - 1), 2), "%",
    "1000 K, pair sums (pre-registered |D| <= 2.5 %: holds); solver values give "
    f"{100 * (k(1000, 20, 'kappa_C') / k(1000, 16, 'kappa_C') - 1):.2f} %", "d2_summary.csv; docs/recorded_rules.md (R8)")
top = d2s[d2s.T_K == 300].set_index("mesh").top20_share
add("top20_share_300", "; ".join(f"{m}^3: {top[m]:.2f}" for m in (15, 16, 17, 20, 24)), "",
    "share of kappa_C(300 K) in the top-20 (q, branch-pair) terms (SI descriptive)", "d2_summary.csv")
add("top20_range_300", f"{100 * min(top[m] for m in (15, 16, 17, 20, 24)):.0f}-{100 * max(top[m] for m in (15, 16, 17, 20, 24)):.0f}",
    "%", "share of kappa_C(300 K) in the 20 largest pair terms, range over 15^3-24^3", "d2_summary.csv")
add("dkC_12_16_pct", round(100 * (k(300, 16, "kappa_C") / k(300, 12, "kappa_C") - 1), 1), "%", "300 K",
    "kappa_mesh_T.csv")
add("dkC_16_20_pct", round(100 * (k(300, 20, "kappa_C") / k(300, 16, "kappa_C") - 1), 1), "%", "300 K",
    "kappa_mesh_T.csv")
e16 = eta[(eta.mesh == 16) & eta.T_K.isin(MAIN)].set_index("T_K")
add("eta_max_300", round(float(e16.eta_max[300]), 3), "", "max Gamma3/Omega, 16^3, 300 K (R-point AFD mode)",
    "eta_validity_T.csv")
add("eta_max_1000", round(float(e16.eta_max[1000]), 3), "", "16^3, 1000 K", "eta_validity_T.csv")
n3 = rd(PJ / "n3_q1_16_check.csv").set_index(["T_K", "quantity"]).shift_pct
add("q1_250_shift", f"Omega_R {n3[(250, 'omega_R')]:.0f} %; kappa_P {n3[(250, 'kappa_P')]:.1f} %", "",
    "250 K, SCPH q1 12^3 -> 16^3 (12^3 transport mesh): the D0 rule fails (kappa_P > 3 %, R > 5 %)",
    "n3_q1_16_check.csv")
_hi = eta[(eta.T_K > 850) & eta.mesh.isin([12, 16])]
_f = sorted({int(round(v, -1)) for v in _hi.eta_max_omega_cm1})
add("most_damped_hard_cm", "; ".join(f"{v}" for v in _f), "cm-1",
    "frequencies (rounded to 10 cm-1) of the most damped mode at 900-1000 K on 12^3 / 16^3 (hard optical)",
    "eta_validity_T.csv")
_e1000 = eta[(eta.T_K == 1000) & eta.mesh.isin([12, 16])].eta_max
add("eta_max_1000_range", f"{_e1000.min():.2f}-{_e1000.max():.2f}", "", "largest three-phonon Gamma/Omega at 1000 K "
    "on 12^3 and 16^3 (hard optical modes)", "eta_validity_T.csv")
_e300 = eta[(eta.T_K == 300) & eta.mesh.isin([12, 16, 20])].set_index("mesh").eta_max
add("eta_max_global", f"{_e300.min():.2f}-{_e300.max():.2f} at 300 K (R mode; 16^3, 12^3, 20^3)", "",
    "largest three-phonon Gamma/Omega over 300-1000 K, at 300 K on every production mesh with all T or 300 K "
    f"({', '.join(f'{m}^3: {v:.3f}' for m, v in _e300.items())}); 16^3 max over 300-1000 K: "
    f"{e16.eta_max.max():.3f} at {int(e16.eta_max.idxmax())} K", "eta_validity_T.csv")
add("f005_300", round(100 * float(e16.frac_eta_gt_0p05_kappa[300]), 2), "%", "kappa-weighted fraction eta > 0.05, 300 K, 16^3",
    "eta_validity_T.csv")
add("f005_1000", round(100 * float(e16.frac_eta_gt_0p05_kappa[1000]), 2), "%", "1000 K, 16^3", "eta_validity_T.csv")
above = [T for T in MAIN if e16.eta_max[T] > 0.1]
add("T_eta01", f"eta_max > 0.1 at every T 300-1000 K on 16^3 ({len(above)}/{len(MAIN)})" if len(above) == len(MAIN)
    else ",".join(map(str, above)), "K", "temperatures where some mode has eta > 0.1 (definition to confirm)",
    "eta_validity_T.csv", "CHOICE")
add("most_damped_sector", "R-point AFD mode (300-850 K); hard optic near 282 cm-1 above 850 K", "",
    "mode at eta_max; label 'A' (branch rule) / 'hard' (character rule)", "eta_validity_T.csv; kappa_mesh_T.csv",
    "CHOICE")
# Soft-sector linewidth sensitivity (production 16^3, 300 K): measured / calculated zone-centre width. Vogt 1995 gives
# a damped-harmonic-oscillator damping constant (hyper-Raman; p.4, Fig. 6), which is a FWHM in the usual DHO form;
# our Gamma is a HWHM. Both conventions are reported.
vg = rd(ROOT / "data/literature/STO/vogt_1995/gamma_T.csv"); vg = vg[vg.phase_mode == "cubic_F1u"]
vo = rd(ROOT / "data/literature/STO/vogt_1995/omega_T.csv"); vo = vo[vo.phase_mode == "cubic_F1u"]
mt = rd(PROD / "mode_table_T.csv.gz")
x16 = mt[(mt.T_K == 300) & (mt.mesh == 16)]
to1 = x16[(x16.q_frac_x.abs() < 1e-6) & (x16.q_frac_y.abs() < 1e-6) & (x16.q_frac_z.abs() < 1e-6) & (x16.branch == 4)].iloc[0]
g_exp = float(np.interp(300, vg.T_K, vg.gamma_meV)) * 8.06554
r_om = float(np.interp(300, vo.T_K, vo.Omega_meV)) * 8.06554 / to1.omega_cm1
soft = (x16.sector_char == "soft").to_numpy()
rg = {c: g_exp * f / to1.gamma_hwhm_cm1 for c, f in (("FWHM", 0.5), ("HWHM", 1.0))}
emr = {c: (np.nanmax(np.where(soft, r * x16.eta, x16.eta)), np.nanmax(np.where(soft, r * x16.eta / r_om, x16.eta)))
       for c, r in rg.items()}
add("gamma_ratio_hwhm_reading", round(rg["HWHM"], 1), "", "SI only: ratio if Vogt's DHO damping constant were read "
    f"as a HWHM; largest rescaled ratio then {emr['HWHM'][0]:.2f} (Gamma side) / {emr['HWHM'][1]:.2f} (with Omega)",
    "vogt_1995/gamma_T.csv; mode_table_T.csv.gz", "FINAL (SI only)")
add("gamma_ratio", round(rg["FWHM"], 1), "",
    "convention: DHO damping constant gamma = FWHM = 2 Gamma_HWHM; Vogt 1995 cubic F1u damping {g_exp:.1f} cm-1 at ~298 K over our Gamma-TO1 HWHM {to1.gamma_hwhm_cm1:.3f} cm-1 "
    "(three-phonon, 16^3, 300 K)", "vogt_1995/gamma_T.csv; mode_table_T.csv.gz")
add("eta_max_rescaled", f"{emr['FWHM'][0]:.2f}; {emr['FWHM'][1]:.2f}", "",
    f"FWHM / HWHM reading: {emr['FWHM'][0]:.2f} / {emr['HWHM'][0]:.2f} (Gamma side); {emr['FWHM'][1]:.2f} / "
    f"{emr['HWHM'][1]:.2f} (Gamma and Omega side, Omega_Vogt/Omega_calc = {r_om:.2f}). max eta over all modes after scaling soft-sector (character rule) eta by gamma_ratio (FWHM / HWHM reading), and "
    "additionally by Omega_calc/Omega_Vogt; 16^3, 300 K; FWHM reading (Gamma side; with Omega); all < 0.5",
    "vogt_1995; mode_table_T.csv.gz")
for lab in ("char", "freq"):
    s = sec[(sec.T_K == 300) & (sec.mesh == 16) & (sec.labels == lab)].iloc[0]
    st = "FINAL" if lab == "char" else "FINAL (SI robustness)"
    for ph, col in (("sh_A_300", "kappa_P_A"), ("sh_H_300", "kappa_P_hard"), ("sh_S_300", "kappa_P_soft")):
        add(ph + ("" if lab == "char" else "_freqrule"), round(100 * s[col] / s.kappa_P, 1), "%",
            f"share of kappa_P(300 K, 16^3), {lab} labels", "kappa_sectors_T.csv", st)
m16 = km[(km.mesh == 16) & km.T_K.isin(MAIN)]
sh = 100 * m16.kappa_C / m16.kappa_total
add("kC_share_range", f"{sh.min():.1f}-{sh.max():.1f}", "%", "kappa_C/kappa_total, 16^3, 300-1000 K", "kappa_mesh_T.csv")
for lab in ("char", "freq"):
    c = coh[(coh.mesh >= 15) & (coh.labels == lab) & coh.T_K.isin(MAIN)]
    hh = 100 * c.kC_hard_hard / c.kappa_C_pairsum
    add("hh_share_range" + ("" if lab == "char" else "_freqrule"), f"{hh.min():.0f}-{hh.max():.0f}", "%",
        f"hard/hard share of kappa_C, {lab} labels, meshes >= 15^3 with eigenvectors, 300-1000 K "
        f"(units: {len(c)} (T, mesh); 24^3 has no eigenvector job, so it is not in this range)", "coherence_blocks_T.csv",
        "FINAL")
e1 = rd(PJ / "e1_prod_block_factors.csv")
n = e1[(e1.mesh == 16) & (e1.T_K == 300) & (e1.labels == "char") & (e1.block == "hard/hard")].n_pairs
add("n_hh_pairs", int(n.iloc[0]), "", "hard/hard pairs (irreducible q, branch pair), 16^3, 300 K, char labels",
    "e1_prod_block_factors.csv")
# Four-phonon audit, soft modes at 300 K (R4; tolerance check R7)
a1 = rd(PJ / "a1_per_mesh.csv")
a1 = a1[a1.mesh.isin([6, 8])]  # the meshes that enter the §8 rule
a1d = rd(PJ / "a1_decision.csv").set_index("mode")
star = lambda df: df[["q_frac_x", "q_frac_y", "q_frac_z"]].abs().apply(lambda r: tuple(sorted(r)), axis=1)
SOFT = {"Gamma-TO1": ((0, 0, 0), 4), "M": ((0, 0.5, 0.5), 1), "R": ((0.5, 0.5, 0.5), 1)}
for mode, (st, b) in SOFT.items():
    g3 = []
    for mesh in (12, 16, 20):
        x = mt[(mt.T_K == 300) & (mt.mesh == mesh)]
        r = x[(star(x) == st) & (x.branch == b)].iloc[0]
        g3.append(r.gamma_hwhm_cm1); om = r.omega_cm1
    g4 = a1[a1["mode"] == mode].Gamma4
    tag = {"Gamma-TO1": "TO1", "M": "M", "R": "R"}[mode]
    lo, hi = (min(g3) + g4.min()) / om, (max(g3) + g4.max()) / om
    dig = 2 if mode == "R" else 3
    add(f"eta34_{tag}_300", f"{lo:.{dig}f}-{hi:.{dig}f}", "",
        f"eta_3+4 = (Gamma3 + Gamma4)/Omega at 300 K, {mode}: Gamma4 over the 6^3/8^3 audit jobs, Gamma3 over the "
        f"production meshes 12/16/20^3 ({', '.join(f'{v:.3f}' for v in g3)} cm-1); 16^3 Gamma3 alone: "
        f"{(g3[1] + g4.min()) / om:.3f}-{(g3[1] + g4.max()) / om:.3f}" + (" (in this model)" if mode == "R" else ""),
        "a1_per_mesh.csv; mode_table_T.csv.gz; R4")
    add(f"a1_outcome_{tag}", a1d.outcome[mode], "", f"pre-registered rule (R4) outcome, {mode}; "
        f"Gamma4_max/Gamma4_thr = {a1d.Gamma4_max[mode] / a1d.Gamma4_thr_prereg[mode]:.3f}", "a1_decision.csv")
add("g4g3_soft_300", f"{a1.G4_over_G3.min():.1f}-{a1.G4_over_G3.max():.1f}", "",
    "Gamma4/Gamma3 at the same settings, soft modes (Gamma-TO1, M, R), 6^3/8^3 audit jobs, 300 K", "a1_per_mesh.csv")
_a2 = rd(PJ / "a1_per_mesh_T1000.csv"); _cal = pd.concat([a1.Gamma3_calib_pct, _a2.Gamma3_calib_pct])
add("g3_calib_range", f"{_cal.min():.0f} to +{_cal.max():.0f}", "%",
    "Gamma3 on the 6^3/8^3 audit meshes relative to production 16^3: soft modes at 300 K, soft + X b7 / R b13 at "
    f"1000 K (300 K alone: {a1.Gamma3_calib_pct.min():.0f} to +{a1.Gamma3_calib_pct.max():.0f})",
    "a1_per_mesh.csv; a1_per_mesh_T1000.csv")
emx = eta[eta.T_K >= 300]
ei = emx.eta_max.idxmax()
add("eta3_max_any", round(float(emx.eta_max.max()), 3), "",
    f"largest three-phonon eta at >= 300 K on any production mesh ({emx.mesh[ei]}^3, {int(emx.T_K[ei])} K, "
    "R-point AFD mode)", "eta_validity_T.csv")
# meshes outside the damping map (no eigenvector data): largest eta from their per-mode widths (eta_max_outside_map.py)
_mx = rd(PJ / "eta_max_outside_map.csv").groupby("mesh").eta_max.max().to_dict()
add("eta3_max_8_24", f"{_mx[8]:.3f}; {_mx[24]:.3f}", "", "largest three-phonon eta from the per-mode widths of the 8^3 "
    "(300-1000 K) and 24^3 (300 K) runs, which are not part of the damping map", "eta_max_outside_map.csv")
_e20 = eta[(eta.mesh == 20) & (eta.T_K == 800)].iloc[0]
add("eta20_800_omega", f"{float(_e20.eta_max_omega_cm1):.0f}", "cm-1", "largest-eta mode on 20^3 at 800 K (q index "
    f"{int(_e20.eta_max_q_index)}, branch {int(_e20.eta_max_branch)}, {_e20.eta_max_sector_char} optic)", "eta_validity_T.csv")
_sn = rd(ROOT / "data/processed/si_v4/si_numbers_v4.csv").set_index("placeholder")
add("strict_diag_R", _sn.value["strict_diag_R"], "", _sn.description["strict_diag_R"], "si_numbers_v4.csv (" +
    _sn.source["strict_diag_R"] + ")")
for _k in ("bench_R2", "bench_R_ge8", "bench_G_max", "bench_G_8", "bench_M_max"):   # Sec. 3.1 = SI Table S27
    add(_k, _sn.value[_k], "%", _sn.description[_k], "si_numbers_v4.csv (" + _sn.source[_k] + ")")
add("g4g3_required", round(0.5 / float(emx.eta_max.max()) - 1, 2), "",
    "Gamma4/Gamma3 a mode would need for eta_3+4 = 0.5: 0.5/eta3_max_any - 1 (ratio form; superseded in the main "
    "text by the additive eta4 form; SI supporting table only)", "eta_validity_T.csv",
    "FINAL (SI only)")
q129 = mt[(mt.mesh == 16) & (mt.T_K == 1000) & (mt.q_index == 129) & (mt.branch == 7)].iloc[0]
add("g4g3_required_q129", round(0.5 / q129.eta - 1, 1), "",
    f"same for q129 branch 7 at 1000 K (16^3 only; eta3 = {q129.eta:.3f}), the 1000 K eta3 maximum (SI only)",
    "mode_table_T.csv.gz", "FINAL (SI only)")
e4 = rd(PJ / "eta4_margin_cases.csv")
add("g4g3_max_sampled", round(float(e4.G4_over_G3.max()), 2), "",
    "largest Gamma4/Gamma3 (same settings) over all sampled cases (Gamma-TO1, 1000 K, 6^3, eps 2); "
    "ratio statistic, SI only", "eta4_margin_cases.csv", "FINAL (SI only)")
tol = []
for e in (2, 5):
    ja, jb = man[f"m4ph_audit_q12_T300_m6_e{e}"]["metrics"]["modes"], man[f"m4ph_audittol6_q12_T300_m6_e{e}"]["metrics"]["modes"]
    tol += [abs(jb[k]["Gamma4_hwhm"][0] / ja[k]["Gamma4_hwhm"][0] - 1) for k in ja]
add("tol_sens_max_pct", round(100 * max(tol), 1), "%",
    "largest |change| of Gamma4 when the identity tolerance is widened from 1e-9 to 1e-6 Ry (6^3, eps 2/5, soft modes; "
    "all changes are decreases)", "data/provenance/manifest_final.json m4ph_audittol6_*; R7")
# Four-phonon at 1000 K (R5) and N4 heat carriers (R6); additive margin
a2 = rd(PJ / "a1_per_mesh_T1000.csv")
A2M = {"Gamma-TO1": ((0, 0, 0), 4), "M": ((0, 0.5, 0.5), 1), "R": ((0.5, 0.5, 0.5), 1), "X-b7": ((0, 0, 0.5), 7),
       "R-b13": ((0.5, 0.5, 0.5), 13)}
hi_all = 0
for mode, (st, b) in A2M.items():
    g3 = []
    for mesh in (12, 16):
        x = mt[(mt.T_K == 1000) & (mt.mesh == mesh)]
        r = x[(star(x) == st) & (x.branch == b)].iloc[0]
        g3.append(r.gamma_hwhm_cm1); om = r.omega_cm1
    g4 = a2[a2["mode"] == mode].Gamma4
    lo, hi = (min(g3) + g4.min()) / om, (max(g3) + g4.max()) / om
    hi_all = max(hi_all, hi)
    add(f"eta34_{mode.replace('Gamma-', '').replace('-', '')}_1000", f"{lo:.3f}-{hi:.3f}", "",
        f"eta_3+4 at 1000 K, {mode}: Gamma4 over the 6^3/8^3 audit jobs, Gamma3 over production 12/16^3 "
        f"({', '.join(f'{v:.3f}' for v in g3)} cm-1)" + (" (in this model)" if mode == "R" else ""),
        "a1_per_mesh_T1000.csv; mode_table_T.csv.gz; R5")
add("eta34_soft_1000", f"<= {hi_all:.2f}", "", "largest eta_3+4 at 1000 K over Gamma-TO1, M, R, X b7, R b13 "
    "(R, in this model); all five in the first outcome class of R5", "a1_per_mesh_T1000.csv; a1_decision_T1000.csv")
s2 = a2[a2["mode"].isin(["Gamma-TO1", "M", "R"])]
add("g4g3_soft_1000", f"{s2.G4_over_G3.min():.2f}-{(s2.Gamma4 / s2.Gamma3).max():.2f}", "",
    "Gamma4/Gamma3 at the same settings, soft modes, 6^3/8^3, 1000 K (Gamma-TO1 largest)", "a1_per_mesh_T1000.csv")
add("g4g3_soft_300_2dp", f"{a1.G4_over_G3.min():.2f}-{a1.G4_over_G3.max():.2f}", "",
    "as g4g3_soft_300, two decimals (soft-mode sentence)", "a1_per_mesh.csv")
n4 = rd(PJ / "n4_results.csv")
for T in (300, 800):
    x = n4[n4["T"] == T]
    add(f"g4g3_{T}", f"{x.r.min():.2f}-{x.r.max():.2f}", "", f"Gamma4/Gamma3 at the same settings, twelve N4 heat "
        f"carriers, {T} K, all eps", "n4_results.csv; R6")
    add(f"life_short_{T}", f"{100 * (1 - x.life.max()):.0f}-{100 * (1 - x.life.min()):.0f}", "%",
        f"lifetime shortening 1 - 1/(1 + Gamma4/Gamma3), N4, {T} K", "n4_results.csv")
add("n_pilot", "12 heat carriers (6 on 8^3, 6 on 6^3; plus 3 soft and 3 hard-optic modes)", "",
    "four-phonon sample size", "n4_results.csv; R6")
nk = n4.drop_duplicates(["sample", "q", "b"])
add("n4_kP_share", "1.2-1.3", "%", f"share of kappa_P carried by the N4 sample ({len(nk)} modes; R6)",
    "n4_selection.csv; R6", "FINAL (copied from R6)")
e4t = rd(PJ / "eta4_margin.csv")
top = e4.iloc[e4.eta4.idxmax()]
add("eta4_max", round(float(e4.eta4.max()), 2), "", f"largest Gamma4/Omega over {len(e4.drop_duplicates(['mode', 'T_K']))} "
    f"mode-temperature combinations, {len(e4)} mesh/eps cases ({top['mode']}, {int(top.T_K)} K, {int(top.mesh)}^3, "
    f"eps {int(top.eps)}: {top.eta4:.4f}; in this model)", "eta4_margin_cases.csv")
add("n_combos_4ph", len(e4.drop_duplicates(["mode", "T_K"])), "", f"mode-temperature combinations with diagram (c) "
    f"evaluated ({len(e4)} mesh/eps cases on 6^3/8^3)", "eta4_margin_cases.csv")
hc = e4[e4.group.str.startswith("heat")]
add("eta4_max_hc", round(float(hc.eta4.max()), 3), "", "largest Gamma4/Omega, heat carriers (800 K)", "eta4_margin_cases.csv")
for r in e4t.itertuples():
    add(f"eta4_needed_{r.T_K}", round(r.eta4_needed, 2), "", f"0.5 - largest production eta3 at {r.T_K} K "
        f"({r.eta3_max:.3f}, {r.eta3_max_mesh}^3); {r.needed_over_eta4_max:.1f}x the largest sampled eta4",
        "eta4_margin.csv")
v300 = vogt[vogt.phase_mode == "cubic_F1u"].iloc[(vogt[vogt.phase_mode == "cubic_F1u"].T_K - 300).abs().argmin()]
w0 = float(rpc[(rpc.PM == "on") & (rpc.q1 == 12) & (rpc.T_K == 300)].Gamma.iloc[0])
add("Omega0_ratio", round(w0 / (v300.Omega_meV * 8.06554), 2), "",
    f"Omega_TO1(300 K) ours {w0:.1f} cm-1 / Vogt 1995 {v300.Omega_meV} meV at {v300.T_K} K (read from figure)",
    "r_point_check.csv; data/literature/STO/vogt_1995/omega_T.csv")
p = d4[(d4.variant == "production_pm_q12") & d4.T_K.isin(MAIN)]
add("bare_factor_range", f"{p.ratio_full.min():.2f}-{p.ratio_full.max():.2f} (full set); "
    f"{p.ratio_subset_char_frac.min():.2f}-{p.ratio_subset_char_frac.max():.2f} (character-matched); "
    f"{p.ratio_absw.min():.2f}-{p.ratio_absw.max():.2f} (|w| substitution)", "", "SCPH/bare kappa_P, 300-1000 K",
    "d4_bare_framings.csv")
# Benchmark against the published value and experiment (Sec. 5.1; SI Sec. S7b)
KPUB = 8.7  # Tadano 2015 Fig. 7 (12^3 BTE, q1 8^3), tracked row "digitized_Fig7_12x12x12"; re-read 8.71 +- 0.11
add("kP_pub", "8.7", "W/mK", "published SCPH kappa_P(300 K), read from figure (Tadano 2015 Fig. 7, +-0.11; the earlier "
    "+-0.3 has no recorded basis); Guo & Guo 2023 p.6 quote Tadano's value as 8.75 W/mK", "SI Sec. S7b and Sec. 5.1 (exp_datasets.csv); "
    "SI Sec. S10 (prior_work_comparison.md)")
add("kP_pub_settings", 9.24, "W/mK", "ours at the published settings (q1 8^3, PM, 12^3 mesh), 300 K",
    "run rta_on_k8_m12 / rta_pm_q8_T300_m12")
add("kP_pub_settings_pct", f"+{100 * (9.2384 / KPUB - 1):.0f}", "%", "ours at published settings vs 8.7", "")
add("kP_vs_pub_pct", f"+{100 * (kP_best / KPUB - 1):.0f}", "%", "best estimate vs 8.7", "kappa_mesh_T.csv")
xg = rd(PJ / "exp_gap_vs_T.csv")
xr = rd(PJ / "exp_300K_crossread.csv")
m300 = xr[(xr.dataset == "Muta2005") & xr.T_K.between(295, 305)].kappa_W_mK
kmu = (min(float(xg[xg.T_K == 300]["k_Muta2005"].iloc[0]), m300.min()), max(m300.max(), 11.0))
kmu = (round(kmu[0], 2), round(m300.max(), 2))
gaps = [100 * (1 - t / ke) for ke in (11.0, *kmu) for t in (tot5.min(), tot5.max())]
gP = [100 * (1 - t / ke) for ke in (11.0, *kmu) for t in (big5.kappa_P.min(), big5.kappa_P.max())]
add("exp_gap_range", f"{min(gaps):.0f}-{max(gaps):.0f} %", "%",
    f"kappa_total(300 K) observed spread {tot5.min():.2f}-{tot5.max():.2f} below the single-crystal values (Muta 2005 "
    f"{kmu[0]:.2f}-{kmu[1]:.2f}, read from three figures; Tachibana 2008 11); kappa_P {min(gP):.0f}-{max(gP):.0f} % below; "
    "with Tachibana's +-1: 3.5-22 %", "exp_300K_crossread.csv; exp_gap_vs_T.csv; SI Sec. S7b and Sec. 5.1 (exp_datasets.csv)")
add("exp_gap_kP", f"{min(gP):.0f}-{max(gP):.0f} %", "%", "kappa_P(300 K), >= 15^3 meshes, below the same values",
    "exp_300K_crossread.csv")
add("exp_datasets_300", f"{kmu[0]:.1f}-{kmu[1]:.1f}", "W/mK",
    f"Muta 2005 at 300 K, {kmu[0]:.2f}-{kmu[1]:.2f}, read from figure in Tadano Fig. 7, Wang Fig. 4 and Guo & Guo "
    "Fig. 10 (+-0.13-0.32; reduced/La-doped single crystal per its title); Tachibana 2008: 11 +- 1 (tracked CSV); "
    "Popuri 2014 has no point below 321 K (9.8 at 321-325 K)",
    "exp_300K_crossread.csv; exp_datasets.csv")
# Guo & Guo 2023 (prior_work_comparison.md; values stated on pp.5-6): four-phonon reduction computed from their stated values
_gg = {"300": [(13.9, 8.7), (13.0, 8.4)], "900": [(9.9, 5.5), (9.8, 5.4)]}
_r = {T: [100 * (1 - b / a) for a, b in v] for T, v in _gg.items()}
add("g4_lit_guo", f"{min(_r['300']):.0f}-{max(_r['300']):.0f} % (300 K); {min(_r['900']):.0f}-{max(_r['900']):.0f} % "
    "(900 K)", "%", "kappa reduction by four-phonon scattering, computed from Guo & Guo's stated 3ph and 3+4ph values "
    "(ITA / RTA)", "prior_work_comparison.md (Guo & Guo p.5)")
add("g4_lit_band", f"15-20; {min(min(v) for v in _r.values()):.0f}-{max(max(v) for v in _r.values()):.0f}", "%", "literature band of the four-phonon kappa reduction",
    "SI Sec. S10 (prior_work_comparison.md)")
mfp = rd(PJ / "mfp_vs_eta.csv")
for T in (300, 900):
    z = mfp[(mfp.T_K == T) & mfp.mfp_cut.str.startswith("2.8")].iloc[0]
    add(f"mfp28_{T}", f"{100 * z.frac_modes:.0f} % of modes; {100 * z.kappaP_share:.1f} % of kappa_P; eta <= "
        f"{z.eta_max:.2f} (median {z.eta_median:.3f})", "", f"modes with MFP = |v| tau < 2.8 A, 16^3, {T} K",
        "mfp_vs_eta.csv")
wc = rd(PJ / "wang_comparisons.csv")
add("alpha_300_800", "0.50-0.51", "", "kappa_P ~ T^-alpha, 300-800 K, meshes 12/16/20 (Muta read from figure: 0.69)",
    "wang_comparisons.csv; SI Sec. S7b and Sec. 5.1 (exp_datasets.csv)")
add("Tc_R_fit", 245, "K", "mean-field fit Omega_R^2 = A(T - Tc), 300-600 K, q1 = 12^3", "r_point_check.csv")
# ---------------------------------------------------------------- draft 5
t300 = km[(km.T_K == 300)].set_index("mesh")
dense = t300[t300.index >= 15]
add("kP_ref_16", f"{k(300, 16, 'kappa_P'):.2f}", "W/mK", "designated reference kappa_P(300 K): 16^3, the only mesh >= 15^3 at every T",
    "kappa_mesh_T.csv")
_kp = t300[t300.index >= 12].kappa_P
add("kP_300_tested", f"{_kp.min():.2f}-{_kp.max():.2f}", "W/mK", "kappa_P(300 K) over the tested meshes 12^3-24^3",
    "kappa_mesh_T.csv")
add("kC_300_tested", f"{dense.kappa_C.min():.3f}-{dense.kappa_C.max():.3f}", "W/mK",
    "kappa_C(300 K) over the tested meshes 15^3-24^3 (observed range, not an error bar)", "kappa_mesh_T.csv")
_tot = dense.kappa_P + dense.kappa_C
add("ktot_300_paired", f"{_tot.min():.2f}-{_tot.max():.2f}", "W/mK", "kappa_P + kappa_C paired per mesh, 15^3-24^3, 300 K",
    "kappa_mesh_T.csv")
# reference 16^3; three separate comparisons with the published 8.7
_pub = 8.7
_ref = k(300, 16, 'kappa_P')
add("kP_ref_pub_pct", f"+{100 * (_ref / _pub - 1):.0f}", "%", "reference kappa_P(300 K, 16^3) vs published 8.7 (read from figure)",
    "kappa_mesh_T.csv; kP_pub")
add("kP_tested_pub_pct", f"+{100 * (_kp.min() / _pub - 1):.0f} to +{100 * (_kp.max() / _pub - 1):.0f}", "%",
    "tested meshes 12^3-24^3 vs published 8.7", "kappa_mesh_T.csv; kP_pub")
_t16 = k(300, 16, 'kappa_P') + k(300, 16, 'kappa_C')
add("ktot_ref_16", f"{_t16:.2f}", "W/mK", "kappa_P + kappa_C at the reference 16^3 mesh, 300 K", "kappa_mesh_T.csv")
add("exp_gap_ref", f"{100 * (1 - _t16 / 11.05):.0f}", "%", "reference total below the midpoint of 11.0-11.1 "
    f"({100 * (1 - _t16 / 11.0):.1f}-{100 * (1 - _t16 / 11.1):.1f} %)", "kappa_mesh_T.csv; exp_datasets.csv")
_e = [float(x) for x in "11.0-11.1".split("-")]
add("exp_gap_paired", f"{100 * (1 - _tot.max() / _e[0]):.0f}-{100 * (1 - _tot.min() / _e[1]):.0f}", "%",
    "paired totals below the 300 K single-crystal values 11.0-11.1 (smallest gap: max total vs 11.0; largest: min total vs 11.1)",
    "kappa_mesh_T.csv; exp_datasets.csv")
_c = coh[(coh.mesh >= 15) & (coh.T_K >= 300)].copy()
_c["hh"] = 100 * _c.kC_hard_hard / _c.kappa_C_pairsum
for lab, nm in (("char", "hh_char_15up"), ("freq", "hh_freq_15up")):
    z = _c[_c.labels == lab].hh
    add(nm, f"{z.min():.1f}-{z.max():.1f}", "%", f"hard/hard share of kappa_C, {lab} rule, every mesh >= 15^3, 300-1000 K "
        f"({len(z)} temperature-mesh units)", "coherence_blocks_T.csv")
bf = rd(PJ / "e1_prod_block_factors.csv")
_b = bf[(bf.mesh == 16) & (bf.T_K == 300) & (bf.labels == "char")].set_index("block")
hh, hs = _b.loc["hard/hard"], _b.loc["hard/soft"]
add("vsum_ratio_hh_hs", f"{hh.sum_V / hs.sum_V:.2f}", "", "ratio of summed velocity factors V, hard/hard : hard/soft, 16^3, "
    "300 K, character rule", "e1_prod_block_factors.csv")
add("npairs_hh_hs", f"{int(hh.n_pairs)}; {int(hs.n_pairs)}", "", "ordered branch pairs s != s' at the 165 irreducible q of 16^3, "
    "hard/hard; hard/soft (300 K, character rule)", "e1_prod_block_factors.csv")
add("vpair_hh_hs", f"{hh.sum_V / hh.n_pairs:.1f}; {hs.sum_V / hs.n_pairs:.1f}", "W/mK/cm^-1", "summed V per ordered pair, "
    "hard/hard; hard/soft", "e1_prod_block_factors.csv")
_b1 = bf[(bf.mesh == 16) & (bf.T_K == 1000) & (bf.labels == "char")].set_index("block")
add("H_hh_300_1000", f"{hh.kw_mean_H:.2f}; {_b1.loc['hard/hard'].kw_mean_H:.2f}", "", "kappa-weighted mean H of hard/hard pairs, "
    "300; 1000 K (H -> 2 in the classical limit)", "e1_prod_block_factors.csv")
fz = rd(PJ / "e1_prod_factor_freezing.csv")
fz = fz[(fz.mesh == 16) & (fz.T_ref == 300) & (fz.labels == "char") & (fz.block == "hard/hard")].pivot(
    index="T_K", columns="frozen", values="kappa_block")
add("ff_hh_none", f"{fz['none'].min():.3f}-{fz['none'].max():.3f}", "W/mK", "hard/hard kappa_C, 16^3, character rule, 300-1000 K",
    "e1_prod_factor_freezing.csv")
add("ff_hh_1000", f"{fz.loc[1000, 'none']:.3f}; {fz.loc[1000, 'H']:.3f}; {fz.loc[1000, 'L']:.3f}; {fz.loc[1000, 'V']:.3f}",
    "W/mK", "hard/hard at 1000 K: all factors at T; H; L; V held at 300 K", "e1_prod_factor_freezing.csv")
s3p = rd(ROOT / "data/processed/si_v4/s3_pairs.csv")
_hi = s3p[s3p.T_K.isin([650, 800, 1000]) & (s3p.pair == "16³→20³")]
add("kC_16_20_highT_max", f"{_hi.dkC_pct.abs().max():.1f}", "%", "largest |kappa_C(20^3)/kappa_C(16^3) - 1| at 650, 800, 1000 K",
    "s3_pairs.csv")
add("r_900_950", "1.09; 1.10", "", "kappa_C-weighted median r on 16^3 at 900; 950 K (no finer mesh run)", "d2_summary.csv")
_ev = eta[(eta.T_K >= 300)]
add("f_eta01_kappa_max", f"{100 * _ev.frac_eta_gt_0p1_kappa.max():.2f}", "%", "largest kappa-weighted fraction of modes with "
    "eta_3 > 0.1, any T 300-1000 K and mesh", "eta_validity_T.csv")
sd = rd(PJ / "sector_degeneracy.csv")
sd = sd[sd.T_K >= 300]
add("deg_mixed", f"{100 * sd.n_subspaces_mixed.sum() / sd.n_subspaces_deg.sum():.0f}; {sd.kappaP_share_mixed_pct.max():.2f}; "
    f"{(sd.hh_share_inv_pct - sd.hh_share_char_pct).abs().max():.1f}", "%; %; pp",
    "exactly degenerate subspaces whose members get different per-eigenvector labels (share); largest kappa_P share in them; "
    "largest change of the hard/hard share with subspace-summed labels", "sector_degeneracy.csv")
rl = rd(ROOT / "data/processed/si_v4/s5_r_mode_label.csv").set_index("T_K").loc[300]
add("R_label_300", f"{rl.omega_cm1:.1f}; {rl.c_A:.2f}; {rl.c_soft:.2f}; {rl.c_hard:.2f}", "cm-1; -; -; -",
    "R-point lowest mode on 16^3 at 300 K: frequency and projections c_A, c_soft, c_hard", "s5_r_mode_label.csv")
add("ovl_frac_hh_hs", f"{100 * hh.frac_pairs_dw_lt_G:.1f}; {100 * hs.frac_pairs_dw_lt_G:.1f}", "%",
    "pairs with |dw| < G1 + G2: hard/hard; hard/soft (16^3, 300 K, character rule)", "e1_prod_block_factors.csv")
add("kpair_hh_hs", f"{1e5 * hh.kappa_block / hh.n_pairs:.2f}; {1e5 * hs.kappa_block / hs.n_pairs:.2f}", "1e-5 W/mK",
    "block contribution per ordered pair, hard/hard; hard/soft (16^3, 300 K, character rule)", "e1_prod_block_factors.csv")
pcn = rd(PJ / "e1_pair_counts.csv")
szs = rd(PJ / "e1_sector_sizes.csv")
for rule in ("char", "freq"):
    z = pcn[(pcn.mesh == 16) & (pcn.T_K == 300) & (pcn.labels == rule)].set_index("block").loc["hard/hard"]
    n = szs[(szs.mesh == 16) & (szs.T_K == 300) & (szs.labels == rule)].iloc[0]
    add(f"hh_pairs_{rule}", f"{100 * z.pair_frac:.0f}; {100 * z.share_of_kC:.0f}; {z.per_pair_rel:.1f}; {n.mean_n_hard:.1f}",
        "%; %; -; -", f"hard/hard, 16^3, 300 K, {rule} rule: share of ordered pairs; share of kappa_C; per-pair contribution "
        "relative to the all-pair mean; mean number of hard modes per irreducible q (of 15)", "e1_pair_counts.csv; e1_sector_sizes.csv")
for rule in ("char", "freq"):   # star-weighted (full-zone pairs)
    z = pcn[(pcn.mesh == 16) & (pcn.T_K == 300) & (pcn.labels == rule)].set_index("block").loc["hard/hard"]
    add(f"hh_pairs_w_{rule}", f"{100 * z.pair_frac_w:.0f}; {100 * z.share_of_kC:.0f}; {z.per_pair_rel_w:.1f}", "%; %; -",
        f"hard/hard, 16^3, 300 K, {rule} rule, star-weighted: share of full-zone ordered pairs; share of kappa_C; per-pair "
        "contribution relative to the all-pair mean", "e1_pair_counts.csv")
_w = pcn[(pcn.mesh == 16) & (pcn.T_K == 300) & (pcn.labels == "char")].set_index("block")
add("npairs_w_hh_hs", f"{int(_w.n_pairs_w['hard/hard'])}; {int(_w.n_pairs_w['hard/soft'])}; {int(_w.n_pairs_w.sum())}", "",
    "full-zone ordered pairs (star-weighted), hard/hard; hard/soft; all (= 210 N_k), 16^3, 300 K, character rule",
    "e1_pair_counts.csv")
add("kpair_w_hh_hs", f"{1e7 * _w.kappa_per_pair_w['hard/hard']:.2f}; {1e7 * _w.kappa_per_pair_w['hard/soft']:.2f}",
    "1e-7 W/mK", "contribution per full-zone ordered pair, hard/hard; hard/soft (16^3, 300 K, character rule)",
    "e1_pair_counts.csv")
add("vpair_w_hh_hs", f"{_w.V_per_pair_w['hard/hard']:.2f}; {_w.V_per_pair_w['hard/soft']:.2f}", "W/mK/cm^-1",
    "summed V per full-zone ordered pair, hard/hard; hard/soft (16^3, 300 K, character rule)", "e1_pair_counts.csv")
# factor freezing: pair sets at 300 and 1000 K, matched by (q_index, branch) on 16^3 (as e1_production.py does)
_m = pd.read_csv(PROD / "mode_table_T.csv.gz", usecols=["q_index", "branch", "T_K", "mesh", "sector_char"])
_m = _m[(_m.mesh == 16) & _m.T_K.isin([300, 1000])]
_h = {T: {q: (d.sort_values("branch").sector_char.to_numpy() == "hard") for q, d in g.groupby("q_index")}
      for T, g in _m.groupby("T_K")}
_n = {"300": 0, "1000": 0, "both": 0}
for q in _h[300]:
    a, b = np.outer(_h[300][q], _h[300][q]), np.outer(_h[1000][q], _h[1000][q])
    np.fill_diagonal(a, False), np.fill_diagonal(b, False)
    _n["300"] += int(a.sum()); _n["1000"] += int(b.sum()); _n["both"] += int((a & b).sum())
add("ff_pairsets", f"{_n['300']}; {_n['1000']}; {_n['both']}; {_n['1000'] - _n['both']}", "",
    "hard/hard ordered pairs on 16^3 (character rule) at 300 K; at 1000 K; at 1000 K and hard/hard at 300 K under the same "
    "(q, branch) indices; at 1000 K entering the sector only there", "mode_table_T.csv.gz")
_z = szs[(szs.mesh == 16) & (szs.T_K >= 300)]
add("n_hard_range", f"{_z.mean_n_hard.min():.1f}-{_z.mean_n_hard.max():.1f}", "", "mean number of hard modes per irreducible q, "
    "16^3, 300-1000 K, both rules", "e1_sector_sizes.csv")
_z = pcn[(pcn.mesh == 16) & (pcn.T_K >= 300) & (pcn.block == "hard/hard")]
add("hh_pairfrac_range", f"{100 * _z.pair_frac.min():.0f}-{100 * _z.pair_frac.max():.0f}", "%", "hard/hard share of ordered "
    "pairs, 16^3, 300-1000 K, both rules", "e1_pair_counts.csv")
mi = rd(ROOT / "data/processed/si_v4/s7_margin_inequality.csv").set_index(["T_K", "mode"])
for m, nm in (("Gamma-TO1", "TO1"), ("M", "M"), ("R", "R")):
    add(f"eta34_{nm}_300_v", f"{mi.loc[(300, m), 'eta34_max']:.3f}", "", f"eta_3 (16^3 production) + eta_4,max, {m}, 300 K "
        "(plain sum; Gamma_4,max over 6^3 and 8^3 runs)", "s7_margin_inequality.csv")
sx = rd(ROOT / "data/processed/si_v4/s7b_experiment.csv").set_index("T_K")
add("R_exp_deficit", "; ".join(f"{-sx.loc[T, 'R_deficit_pct']:.1f}" for T in (300, 600, 1000)), "%",
    "R-point frequency below experiment (read from Fauqué Fig. 1(d), tracked CSV) at 300; 600; 1000 K", "s7b_experiment.csv")
out = pd.DataFrame(rows)
f = ROOT / "data/processed/manuscript_numbers_v8.csv"
with open(f, "w") as fh:
    fh.write("# v8 placeholder values from production only (own SCPH q1 = 12^3, PM); 300-1000 K main text\n"
             "# status FINAL / CHOICE / PENDING\n")
    out.to_csv(fh, index=False)
from fekappa.provenance import stamp  # noqa: E402
stamp(f, __file__, ["data/processed/sto_production/", "data/processed/revision/", "data/processed/si_v4/",
                    "data/provenance/manifest_final.json", "data/literature/STO/"])
pd.set_option("display.width", 250, "display.max_colwidth", 110)
print(out[["placeholder", "value", "status"]].to_string(index=False))
