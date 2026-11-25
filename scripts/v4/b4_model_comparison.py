#!/usr/bin/env python3
"""B4 / T1 comparison table, model A next to model B (H1), 300 K.

Model A values come from tracked tables only (data/processed/...); model B values from the T1 / B3 outputs and
the T1 analysis JSON lines (a3_transport_analysis.py). Rows that were not run (20^3, 24^3; conclusion 3) are
written with value "not run". Nothing is classified by polar or rotational family here.

Usage: b4_model_comparison.py OUT_CSV H1_M16_JSON H1_M12_JSON KP_ETA_SHARES_CSV B3_SOFTMODES_CSV A3_TRANSPORT_CSV HARM_DIR
  HARM_DIR: bare 16^3 MODE=phonons runs, HARM_DIR/A/harm_A_m16.mesh.evec and HARM_DIR/B/harm_B_m16.mesh.evec.
  Run from the repository root (model A values are read from data/processed/).
"""
import json
import sys

import numpy as np
import pandas as pd

from fekappa.evec_io import read_evec

RY_CM = 109737.3157
P = "data/processed"
PTS = {"Gamma": (0, 0, 0), "M": (0.5, 0.5, 0), "R": (0.5, 0.5, 0.5)}


def tracked(path):
    return pd.read_csv(path, comment="#")


def harmonic_lowest(evec, pt, nbr):
    """Lowest nbr frequencies (signed, cm^-1) at a high-symmetry point of a bare 16^3 MODE=phonons run."""
    xk, eig, _ = read_evec(evec)
    q = np.array(PTS[pt])
    d = np.abs(((xk - q + 0.5) % 1.0) - 0.5).sum(axis=1)
    k = int(np.argmin(d))
    assert d[k] < 1e-8
    w = np.sign(eig[k]) * np.sqrt(np.abs(eig[k])) * RY_CM
    return np.sort(w)[:nbr]


def main():
    out, j16, j12, shares, b3csv, a3csv, BH = sys.argv[1:8]
    b16, b12 = json.loads(open(j16).read()), json.loads(open(j12).read())
    sh = pd.read_csv(shares).set_index("run")
    km = tracked(f"{P}/sto_production/kappa_mesh_T.csv")
    km = km[km.T_K == 300].set_index("mesh")
    ev = tracked(f"{P}/sto_production/eta_validity_T.csv")
    ev = ev[ev.T_K == 300].set_index("mesh")
    cb = tracked(f"{P}/sto_production/coherence_blocks_T.csv")
    cb = cb[(cb.T_K == 300) & (cb.mesh == 16)].set_index("labels")
    ah = tracked(f"{P}/v3/a2_harmonic_xcheck.csv")
    ah = ah[ah.alamode_fc2 == "STO222"]
    fs = tracked(f"{P}/v3/d_q1_freq_shift.csv")
    fs = fs[fs.T_K == 300].set_index("point")
    bm = tracked(f"{P}/v3/c1_benchmark_matrix.csv")
    bm = bm[(bm.config == "on_k8") & (bm.transport_mesh == 12)].iloc[0]
    bs = pd.read_csv(b3csv)
    bs = bs[bs.T_K == 300].set_index("run")
    a3 = pd.read_csv(a3csv).set_index("label")

    rows = []

    def add(q, mesh, a, a_src, b, b_src, note=""):
        rows.append({"quantity": q, "mesh": mesh, "model_A": a, "model_A_source": a_src, "model_B_H1": b,
                     "model_B_source": b_src, "note": note})

    add("lattice parameter (Bohr)", "", 7.363, "sto_example template (distributed fit)", 7.35167842615,
        "model B SCPH template / rta.in (T1)")
    add("Born charges and dielectric tensor", "", "distributed example BORN (sha256 e283...)", "",
        "DFPT k6 header (BORN sha256 1961c459...)", "G-B2 decision")
    for pt, nb in (("Gamma", 6), ("R", 3), ("M", 1)):
        a = harmonic_lowest(f"{BH}/A/harm_A_m16.mesh.evec", pt, nb)
        t = ah[ah.point == pt].sort_values("branch").omega_alamode_cm1.to_numpy()[:nb]
        same = bool(np.allclose(a, t, atol=6e-3))
        b = harmonic_lowest(f"{BH}/B/harm_B_m16.mesh.evec", pt, nb)
        add(f"harmonic frequencies at {pt}, lowest {nb} (cm^-1; negative = imaginary)", "exact point, no NA direction",
            ", ".join(f"{x:.2f}" for x in a), "harm_A_m16 (bare 16^3, MODE=phonons)",
            ", ".join(f"{x:.2f}" for x in b), "harm_B_m16 (bare 16^3, MODE=phonons)",
            "model A equals data/processed/v3/a2_harmonic_xcheck.csv (STO222)" if same else
            "a2_harmonic_xcheck.csv approaches Gamma along [100] (LO split): " + ", ".join(f"{x:.2f}" for x in t))
    for q1, run, col in ((8, "H1_q8", "omega_q1_8_cm1"), (12, "H1_q12", "omega_q1_12_cm1")):
        for pt in ("Gamma", "M", "R"):
            add(f"SCPH frequency at {pt}, 300 K, PM (cm^-1)", f"q1 {q1}^3", fs.loc[pt, col],
                f"{P}/v3/d_q1_freq_shift.csv", round(bs.loc[run, pt], 4),
                f"b3_scph_softmodes.csv ({run})")
    add("kappa_P (W/mK)", "12^3 (q1 12^3)", km.loc[12, "kappa_P"], km.loc[12, "job"], "not run", "",
        "T1 has no 12^3 run at q1 = 12^3")
    add("kappa_P (W/mK)", "16^3 (q1 12^3)", km.loc[16, "kappa_P"], km.loc[16, "job"], b16["kappa_P_solver"],
        "rta_mb_pm_q12_T300_m16")
    add("kappa_P (W/mK)", "20^3 (q1 12^3)", km.loc[20, "kappa_P"], km.loc[20, "job"], "not run", "",
        "the analysis plan: 20^3 not run now")
    add("kappa_P (W/mK)", "12^3 (q1 8^3)", bm.kappa_P, f"{P}/v3/c1_benchmark_matrix.csv (on_k8)",
        b12["kappa_P_solver"], "rta_mb_pm_q8_T300_m12")
    for mesh, src_b, kb in ((12, "rta_mb_pm_q8_T300_m12 (q1 8^3)", b12), (16, "rta_mb_pm_q12_T300_m16", b16)):
        add("kappa_C (W/mK)", f"{mesh}^3", km.loc[mesh, "kappa_C"] if mesh == 16 else bm.kappa_C,
            km.loc[16, "job"] if mesh == 16 else f"{P}/v3/c1_benchmark_matrix.csv (on_k8, q1 8^3)",
            kb["kappa_C_solver"], src_b)
    add("kappa_C (W/mK)", "15^3-24^3", "0.442-0.541", f"{P}/sto_production/kappa_mesh_T.csv", "not run", "")
    add("R frequency, lowest branch (cm^-1)", "16^3 (q1 12^3)", fs.loc["R", "omega_q1_12_cm1"],
        f"{P}/v3/d_q1_freq_shift.csv", b16["R_omega_cm1"], "rta_mb_pm_q12_T300_m16 .result")
    add("R mode HWHM, lowest branch (cm^-1)", "16^3", a3.loc["a345_q12_m16", "R_hwhm_cm1"],
        "A3_transport.csv (production run)", b16["R_hwhm_cm1"], "rta_mb_pm_q12_T300_m16")
    add("eta of R mode, lowest branch", "16^3", a3.loc["a345_q12_m16", "R_eta"], "A3_transport.csv",
        b16["R_eta"], "rta_mb_pm_q12_T300_m16")
    for mesh in (12, 16, 20):
        r = ev.loc[mesh]
        b = "not run"
        if mesh == 16:
            b = f"{sh.loc['H1', 'eta_max']:.6g} (R, q index 165, branch 2, 22.85 cm^-1)"
        add("largest eta and its mode", f"{mesh}^3", f"{r.eta_max:.6g} (q index {int(r.eta_max_q_index)}, "
            f"branch {int(r.eta_max_branch)}, {r.eta_max_omega_cm1} cm^-1)", f"{P}/sto_production/eta_validity_T.csv",
            b, "rta_mb_pm_q12_T300_m16" if mesh == 16 else "",
            "per-branch .result widths inside the degenerate R triplet (before anphon's averaging)")
    add("largest eta", "12^3 (q1 8^3)", "", "", b12["eta_max"], "rta_mb_pm_q8_T300_m12")
    for t in (0.02, 0.05, 0.1):
        a = {0.05: ev.loc[16, "frac_eta_gt_0p05_kappa"], 0.1: ev.loc[16, "frac_eta_gt_0p1_kappa"]}.get(t, "")
        a_src = f"{P}/sto_production/eta_validity_T.csv" if a != "" else ""
        if a == "":
            a, a_src = sh.loc["a345", f"share_kP_eta_gt_{t}"], "kp_eta_300K_m16_shares.csv (a345)"
        add(f"share of kappa_P from modes with eta > {t}", "16^3", a, a_src, sh.loc["H1", f"share_kP_eta_gt_{t}"],
            "kp_eta_300K_m16_shares.csv (H1)")
    for rule in ("freq", "char"):
        a_kc, a_tot = cb.loc[rule, "kC_hard_hard"], cb.loc[rule, "kappa_C_pairsum"]
        add(f"hard/hard kappa_C ({rule} rule, W/mK)", "16^3", a_kc, f"{P}/sto_production/coherence_blocks_T.csv",
            b16[f"kC_hardhard_{rule}"], "rta_mb_pm_q12_T300_m16 + evec_mb_pm_q12_T300_m16")
        add(f"hard/hard share of kappa_C ({rule} rule)", "16^3", a_kc / a_tot,
            f"{P}/sto_production/coherence_blocks_T.csv", b16[f"share_kC_hardhard_{rule}"], "same")
        add(f"hard/hard pairs ({rule} rule; star-weighted)", "16^3",
            f"{int(a3.loc['a345_q12_m16', f'n_pairs_hardhard_{rule}'])} "
            f"({int(a3.loc['a345_q12_m16', f'n_pairs_w_hardhard_{rule}'])})", "A3_transport.csv (production run)",
            f"{b16[f'n_pairs_hardhard_{rule}']} ({b16[f'n_pairs_w_hardhard_{rule}']})", "same")
        add(f"mean kappa_C per hard/hard pair ({rule} rule; star-weighted, W/mK)", "16^3",
            a3.loc["a345_q12_m16", f"kC_per_pair_w_hardhard_{rule}"], "A3_transport.csv",
            b16[f"kC_per_pair_w_hardhard_{rule}"], "same")
    add("hard/hard share of kappa_C (freq rule)", "12^3 (q1 8^3)",
        a3.loc["a345_q8_m12", "share_kC_hardhard_freq"], "A3_transport.csv (A2 rerun)",
        b12["share_kC_hardhard_freq"], "rta_mb_pm_q8_T300_m12", "no eigenvector run at this setting")
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"{len(rows)} rows -> {out}")


if __name__ == "__main__":
    main()
