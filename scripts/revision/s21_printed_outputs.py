"""Full-precision inputs of SI Table S21 (identity tolerance 1e-9 Ry vs 1e-6 Ry; 6³, 300 K, ε = 2 and 5 cm⁻¹):
for each of the four runs, the printed width files (*.Gamma.N: 2Γ₃ and 2Γ₄ in cm⁻¹ as anphon prints them) and the
AUDIT lines of rta.log (excluded-set counters), copied verbatim, with the mode list of the run.
Inputs: data/raw/alamode_sto/revision/pm_q12_4ph_audit{,_tol1e6}/ (zenodo:pm_q12_4ph_audit.tar.gz);
data/provenance/manifest_final.json (output sha256 of each file, checked here).
Output: data/processed/revision/s21_tolerance_printed.txt
"""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/raw/alamode_sto/revision"
man = json.loads((ROOT / "data/provenance/manifest_final.json").read_text())["jobs"]
RUNS = [("1e-9 Ry (rule)", "pm_q12_4ph_audit", f"m4ph_audit_q12_T300_m6_e{e}") for e in (2, 5)] + \
       [("1e-6 Ry (sensitivity)", "pm_q12_4ph_audit_tol1e6", f"m4ph_audittol6_q12_T300_m6_e{e}") for e in (2, 5)]

out = ["# SI Table S21 inputs: printed outputs of the four tolerance runs, verbatim (see the script docstring).",
       "# Each block: job id, identity tolerance, anphon build id, then the files as printed. Γ = Im Σ_code; the",
       "# files give 2Γ (FWHM) and the table gives Γ (HWHM).", ""]
for tol, var, jid in RUNS:
    j, d = man[jid], RAW / var / jid
    out += [f"=== job {jid}; identity tolerance {tol}; build {j['build']}", "--- modes.txt", (d / "modes.txt").read_text().rstrip()]
    out += ["--- rta.log, AUDIT lines"] + [ln.rstrip() for ln in (d / "rta.log").read_text().splitlines() if "AUDIT" in ln]
    for f in sorted(d.glob(f"{jid}.Gamma.*")):
        rec = j["outputs"].get(f.name)
        got = hashlib.sha256(f.read_bytes()).hexdigest()
        if rec and rec != got:
            raise SystemExit(f"{f}: sha256 differs from the manifest")
        out += [f"--- {f.name} (sha256 {got[:12]}{', as in the manifest' if rec else ''})", f.read_text().rstrip()]
    out.append("")
p = ROOT / "data/processed/revision/s21_tolerance_printed.txt"
p.write_text("\n".join(out) + "\n")
print(f"wrote {p.relative_to(ROOT)}: {len(out)} lines")
