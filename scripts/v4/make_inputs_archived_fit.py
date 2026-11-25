#!/usr/bin/env python3
"""Inputs of the runs with the archived 2015 fit (job names p15_*; see data/provenance/v4/run_map.csv).

Every input is the input of an existing production run of model A with named lines replaced; each replacement must
hit exactly once, and the diff against the reference input is written next to it (input.diff). Output: one directory
per job under OUT with the input file, job.env (what the runner needs) and input.diff. Nothing is run here.
The file names of the 2015 study's force constants and Born data are arguments; the files themselves are not read.

Usage (repository root): make_inputs_archived_fit.py OUT ARCH_XML_NAME ARCH_BORN_NAME EXAMPLE_REFERENCE_DIR
                         HARM_A_PH_IN [LABELS_DIR]
  EXAMPLE_REFERENCE_DIR: example/SrTiO3/reference of the ALAMODE checkout (phband.in);
  HARM_A_PH_IN: ph.in of the bare 16^3 model A run harm_A_m16; LABELS_DIR: label tables of the patched runs
  (t3_code32_P_m16.txt, t5_labels_P_m16.txt, t5b_labels_P_m16.txt), when they exist.
"""
import difflib
import shutil
import sys
from pathlib import Path

REPO = Path(".")
A = sys.argv[1:] + [None] * 6
RAW = REPO / "data/raw/alamode_sto/revision"
EX = Path(A[3]) if A[3] else Path("example/SrTiO3/reference")
REF = {"scph8": RAW / "pm_q8/scph_pm_q8/scph.in", "scph12": RAW / "pm_q12/scph_pm_q12/scph.in",
       "rta12": RAW / "pm_q8/rta_pm_q8_T300_m12/rta.in", "rta16": RAW / "pm_q12/rta_pm_q12_T300_m16/rta.in",
       "evec16": RAW / "pm_q12/evec_pm_q12_T300_m16/ph.in", "band": EX / "phband.in",
       "harm16": Path(A[4]) if A[4] else Path("ph.in")}
ARCH = A[1] or "ARCHIVED_FIT.xml"
DIST = "STO_anharm.xml"
FIT = {"dist": DIST, "arch": ARCH}
BORNF = {"dist": "BORN", "p15": A[2] or "ARCHIVED_BORN"}


def sub(s, old, new):
    assert s.count(old) == 1, (old, s.count(old))
    return s.replace(old, new)


def common(s, prefix_old, job, fit, na, born):
    s = sub(s, f"PREFIX = {prefix_old}", f"PREFIX = {job}")
    s = s.replace(f"FCSXML = {DIST} \n", f"FCSXML = {DIST}\n")       # trailing blank of the example inputs
    s = sub(s, f"FCSXML = {DIST}\n", f"FCSXML = {FIT[fit]}\n")
    s = sub(s, "NONANALYTIC = 3; BORNINFO = BORN", f"NONANALYTIC = {na}; BORNINFO = {BORNF[born]}")
    return s


def scph(job, q1, fit, na, born):
    ref = REF["scph12"] if q1 == 12 else REF["scph8"]
    s = ref.read_text()
    s = common(s, "scph_pm_q12" if q1 == 12 else "scph_on_k8", job, fit, na, born)
    if q1 not in (8, 12):
        s = sub(s, "KMESH_SCPH = 8 8 8", f"KMESH_SCPH = {q1} {q1} {q1}")
    return ref, "scph.in", s, {"KIND": "scph", "FIT": fit, "BORNF": BORNF[born]}


def rta(job, scph_job, mesh, fit, na, born, tk=300):
    ref = REF["rta16"] if mesh == 16 else REF["rta12"]
    s = ref.read_text()
    s = common(s, "rta_pm_q12_T300_m16" if mesh == 16 else "rta_on_k8_m12", job, fit, na, born)
    if mesh == 12:
        s = sub(s, "FC2XML = renorm_300K.xml", "FC2XML = renorm.xml")
    if tk != 300:                                                  # P5: other temperatures
        s = sub(s, "TMIN = 300; TMAX = 300;", f"TMIN = {tk}; TMAX = {tk};")
    return ref, "rta.in", s, {"KIND": "rta", "FIT": fit, "BORNF": BORNF[born], "SCPH": scph_job, "TK": str(tk)}


def evec(job, scph_job):
    ref = REF["evec16"]
    s = sub(ref.read_text(), "PREFIX = evec_pm_q12_T300_m16", f"PREFIX = {job}")
    return ref, "ph.in", s, {"KIND": "evec", "FIT": "none", "BORNF": "BORN", "SCPH": scph_job, "TK": "300"}


def patched(job, labels):
    """P3 patched run (T4 / T5 / T5b encoding): the input of the base run p15_rta_q12_T300_m16 unchanged (same PREFIX,
    so that the outputs can be compared byte for byte), patched build, label table copied into the job directory."""
    ref, name, s, env = JOBS["p15_rta_q12_T300_m16"]
    return ref, name, s, dict(env, KIND="patched", LABELS=Path(labels).name, LABELS_SRC=str(labels))


def harm(job):
    """Bare harmonic 16^3 eigenvectors (census input, as T2's harm_A_m16) with the archived XML."""
    ref = REF["harm16"]
    s = sub(ref.read_text(), "PREFIX = harm_A_m16", f"PREFIX = {job}")
    s = sub(s, "FCSXML = STO222.xml\n", f"FCSXML = {ARCH}\n")
    return ref, "ph.in", s, {"KIND": "harm", "FIT": "arch", "BORNF": "BORN"}


def band(job):
    ref = REF["band"]
    s = common(ref.read_text(), "STO222_NA3", job, "arch", 3, "dist")
    return ref, "phband.in", s, {"KIND": "band", "FIT": "arch", "BORNF": "BORN"}


JOBS = {
    # setup check
    "p15_setup_band": band("p15_setup_band"),
    "p15_c_q2": scph("p15_c_q2", 2, "arch", 3, "dist"),
    # P1
    "p15_b_q8": scph("p15_b_q8", 8, "dist", 2, "dist"),
    "p15_c_q8": scph("p15_c_q8", 8, "arch", 3, "dist"),
    "p15_d_q8": scph("p15_d_q8", 8, "arch", 2, "dist"),
    "p15_e_q8": scph("p15_e_q8", 8, "arch", 2, "p15"),
    "p15_e_q2": scph("p15_e_q2", 2, "arch", 2, "p15"),
    # P2
    "p15_c_q8_T300_m12": rta("p15_c_q8_T300_m12", "p15_c_q8", 12, "arch", 3, "dist"),
    "p15_d_q8_T300_m12": rta("p15_d_q8_T300_m12", "p15_d_q8", 12, "arch", 2, "dist"),
    "p15_e_q8_T300_m12": rta("p15_e_q8_T300_m12", "p15_e_q8", 12, "arch", 2, "p15"),
    # P3 (SCPH, eigenvectors, base transport; the patched runs are staged after the labels exist)
    "p15_harm_m16": harm("p15_harm_m16"),
    "p15_scph_q12": scph("p15_scph_q12", 12, "arch", 3, "dist"),
    "p15_evec_q12_T300_m16": evec("p15_evec_q12_T300_m16", "p15_scph_q12"),
    "p15_rta_q12_T300_m16": rta("p15_rta_q12_T300_m16", "p15_scph_q12", 16, "arch", 3, "dist"),
    # P1b
    "p15_e_q4": scph("p15_e_q4", 4, "arch", 2, "p15"),
    "p15_e_q12": scph("p15_e_q12", 12, "arch", 2, "p15"),
}
# P5: kappa(T) of case e at the other temperatures of the 2015 curve
for _t in (200, 400, 500, 600, 700, 800, 900, 1000):
    JOBS[f"p15_e_q8_T{_t}_m12"] = rta(f"p15_e_q8_T{_t}_m12", "p15_e_q8", 12, "arch", 2, "p15", _t)
LAB = Path(A[5]) if A[5] else None
if LAB and (LAB / "t3_code32_P_m16.txt").exists():               # P3 patched runs, once the labels exist
    JOBS["p15_t4_code32_m16"] = patched("p15_t4_code32_m16", LAB / "t3_code32_P_m16.txt")
    JOBS["p15_t5_m16"] = patched("p15_t5_m16", LAB / "t5_labels_P_m16.txt")
    JOBS["p15_t5b_m16"] = patched("p15_t5b_m16", LAB / "t5b_labels_P_m16.txt")


def main():
    out = Path(sys.argv[1])
    for job, (ref, name, s, env) in JOBS.items():
        d = out / job
        d.mkdir(parents=True, exist_ok=True)
        (d / name).write_text(s)
        env = dict(env, INPUT=name, JOB=job)
        src = env.pop("LABELS_SRC", None)
        if src:
            shutil.copyfile(src, d / env["LABELS"])
        (d / "job.env").write_text("".join(f"{k}={v}\n" for k, v in env.items()))
        diff = difflib.unified_diff(ref.read_text().splitlines(), s.splitlines(), str(ref.relative_to(ref.parents[2])),
                                    f"{job}/{name}", lineterm="", n=0)
        (d / "input.diff").write_text("\n".join(diff) + "\n")
        n = sum(1 for l in (d / "input.diff").read_text().splitlines() if l[:1] in "+-" and l[:3] not in ("+++", "---"))
        print(f"{job:24s} {name:10s} changed lines {n // 2:2d}  {env}")


if __name__ == "__main__":
    main()
