#!/usr/bin/env python3
"""Regenerate the outputs that need the files of the 2015 study (not in this release) from the raw run directories.

The files of the 2015 study are arguments; nothing here names them. RUNS_ROOT holds the run directories by run-map
name (P_* runs of the archived fit, and the model A runs A_rta_q12_m16_T300, A_rta_q8_m12_T300, A_scph_q12).
PROC is the processed output of parts a and b (regenerate_v4.py), used for the model A and B columns of the
side-by-side tables.

Usage (repository root): regenerate_v4_held.py RUNS_ROOT PROC HELD_OUT --kappa-curve FILE --mode-file FILE
         [--arch-xml FILE --a-alm-in FILE --a-dfset FILE --b-dfset FILE] [--source-bands DIR]
  --kappa-curve: the 2015 kappa(T) curve (T, kappa); --mode-file: the 2015 per-mode file at 300 K (Omega, Gamma,
  tau); --arch-xml: the archived force constants (ALM XML); --a-alm-in / --a-dfset: model A's alm input and its 40
  training configurations; --b-dfset: model B's 128 test configurations.
Writes HELD_OUT: p1_softmodes.csv, p1b_softmodes.csv, p3_scph_q12_softmodes.csv, p3_kp_eta_300K_m16_gavg_*,
labels/, t4/, t5_P.json, t5b_P.json, p3_sampling_P.json, p2_analysis.json, p5_kappa_T.csv, p6/, p4/, figures/.
"""
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable


def run(script, *args, stdout=None):
    cmd = [PY, str(HERE / script), *map(str, args)]
    print("+", Path(script).name, *(Path(str(a)).name for a in args), flush=True)
    with (open(stdout, "w") if stdout else open(os.devnull, "w")) as fh:
        subprocess.run(cmd, check=True, stdout=fh)


def one(d, pattern):
    hits = sorted(Path(d).glob(pattern))
    assert len(hits) == 1, (d, pattern, hits)
    return hits[0]


def main():
    argv, opts = sys.argv[1:], {}
    for k in ("--kappa-curve", "--mode-file", "--arch-xml", "--a-alm-in", "--a-dfset", "--b-dfset", "--source-bands"):
        if k in argv:
            i = argv.index(k)
            opts[k] = argv[i + 1]
            del argv[i:i + 2]
    R, PROC, H = Path(argv[0]), Path(argv[1]), Path(argv[2])
    for sub in ("labels", "census", "t4", "p6", "p4", "figures", "p5_runs"):
        (H / sub).mkdir(parents=True, exist_ok=True)

    # SCPH soft modes of the archived fit (configuration steps b-e; grids) and under the settings of model A
    run("scph_softmodes.py", H / "p1_softmodes.csv",
        *[f"{c}_q{q}={R / f'P_scph_{c}_q{q}'}" for c, q in (("b", 8), ("c", 2), ("c", 8), ("d", 8), ("e", 8), ("e", 2))])
    run("scph_softmodes.py", H / "p1b_softmodes.csv", *[f"e_q{q}={R / f'P_scph_e_q{q}'}" for q in (4, 12)])
    run("scph_softmodes.py", H / "p3_scph_q12_softmodes.csv", f"arch_q12={R / 'P_scph_q12_settingsA'}",
        f"A_q12={R / 'A_scph_q12'}")

    # damping-ratio shares, model A next to the archived fit (settings of model A)
    run("kp_eta_distribution.py", "--gamma=avg", 300, H / "p3_kp_eta_300K_m16_gavg", f"A={R / 'A_rta_q12_m16_T300'}",
        f"P={R / 'P_rta_q12_m16_T300_settingsA'}")

    # family labels of the archived spectrum and the carrier / band analyses
    ren = one(R / "P_evec_q12_m16_T300_settingsA", "*.mesh.evec")
    run("t2_family_census.py", H / "census/t2_census", f"P={ren}:{one(R / 'P_harm_m16', '*.mesh.evec')}")
    run("t2_freeze_labels.py", H / "census/t2_census_P.csv", ren, H / "labels", "P_m16")
    run("t5_control_labels.py", H / "labels/labels_P_m16.csv", H / "labels/t5_labels_P_m16.txt")
    run("t5_control_labels.py", H / "labels/labels_P_m16.csv", H / "labels/t5b_labels_P_m16.txt", "80,95,105")
    rd = R / "P_partition_q12_m16_T300"
    run("t4_analysis.py", "P", rd / "t3_partition.txt", rd, H / "labels/labels_P_m16.csv", H / "t4")
    for tag, rn, lt in (("t5", "P_partition_q12_m16_T300_bands", "t5_labels_P_m16.txt"),
                        ("t5b", "P_partition_q12_m16_T300_lowbands", "t5b_labels_P_m16.txt")):
        rd = R / rn
        run("t5_control_analysis.py", "P", rd / "t3_partition.txt", rd, H / "labels" / lt, H / "t4/t4_P_modes.csv",
            H / f"{tag}_P.json", H / "labels/labels_P_m16.csv")
    run("p3_sampling.py", "P", H / "labels/labels_P_m16.csv", H / "t5_P.json", H / "t5b_P.json", H / "p3_sampling_P.json")

    # 300 K against the 2015 data; kappa(T); where the difference sits
    if "--kappa-curve" in opts and "--mode-file" in opts:
        run("p2_analysis.py", opts["--kappa-curve"], opts["--mode-file"], H / "p2_analysis.json", "--modes=e",
            f"a={R / 'A_rta_q8_m12_T300'}", *[f"{c}={R / f'P_rta_{c}_q8_m12_T300'}" for c in "cde"])
        inc = H / "p5_runs"
        for t in (200, 300, 400, 500, 600, 700, 800, 900, 1000):
            src = R / (f"P_rta_e_q8_m12_T{t}" if t != 300 else "P_rta_e_q8_m12_T300")
            link = inc / f"p15_e_q8_T{t}_m12"
            if not link.exists():
                link.symlink_to(src.resolve())
        run("p5_kappa_T.py", inc, opts["--kappa-curve"], H / "p5_kappa_T.csv")
        run("p6_kappa_split.py", R / "P_rta_e_q8_m12_T300", opts["--mode-file"], opts["--kappa-curve"],
            H / "p6/p6.json")
        run("figures/make_kappa_T_2015.py", H / "p5_kappa_T.csv", H / "figures/fig_kappa_T_2015.png")

    # forces of the archived fit on model A's training set and model B's test set
    if all(k in opts for k in ("--arch-xml", "--a-alm-in", "--a-dfset", "--b-dfset")):
        run("b1_force_eval.py", opts["--arch-xml"], opts["--a-alm-in"], opts["--a-dfset"], 40, "--xml",
            "--save", H / "p4/pred_arch_A40.npy", stdout=H / "p4/arch_A40.txt")
        run("b1_force_eval.py", opts["--arch-xml"], opts["--a-alm-in"], opts["--b-dfset"], 128, "--xml",
            "--save", H / "p4/pred_arch_Btest.npy", stdout=H / "p4/arch_Btest.txt")

    # Fig. 2 with the archived series
    if "--source-bands" in opts:
        run("figures/make_fig2_softmode_grid.py", PROC / "a/kp_eta_300K_m16_gavg_cumulative.csv",
            PROC / "b/kp_eta_300K_m16_gavg_cumulative.csv", PROC / "b/b3_scph_softmodes.csv", opts["--source-bands"],
            H / "figures", H / "p1_softmodes.csv", H / "p1b_softmodes.csv")
    print("done:", H)


if __name__ == "__main__":
    main()
