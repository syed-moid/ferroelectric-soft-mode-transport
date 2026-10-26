"""The independent count of the four-phonon excluded set (scripts/revision/a1_identity_count.py) reproduces the
counters printed by the audit build for the 4^3 job m4ph_audit_q12_T300_m4_e2 (data/provenance/manifest_final.json)."""
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/data/evec_pm_q12_T300_m4.mesh.evec.gz"


def test_identity_count_matches_audit_build():
    job = json.loads((ROOT / "data/provenance/manifest_final.json").read_text())["jobs"]["m4ph_audit_q12_T300_m4_e2"]
    code = [tuple(int(x) for x in re.search(r"excluded = (\d+) .*denominator check = (\d+)", s).groups())
            for s in job["metrics"]["audit_4ph_c"]]
    out = subprocess.run([sys.executable, str(ROOT / "scripts/revision/a1_identity_count.py"), str(FIXTURE),
                          "0", "0", "0", "4", "0.5", "0.5", "0.5", "1", "0.5", "0.5", "0", "1"],
                         capture_output=True, text=True, check=True).stdout
    ours = [tuple(int(x) for x in re.search(r"excluded (\d+), exact-zero denominators (\d+)", ln).groups())
            for ln in out.splitlines() if ln.startswith("target")]
    assert ours == code == [(11268, 11268), (11268, 11268), (3772, 6616)]
