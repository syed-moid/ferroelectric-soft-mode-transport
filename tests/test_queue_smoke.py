"""The compute-queue unit tests (queue/test_queue.py: validation, recovery, STOP handling, input writers), run without
a solver in a temporary queue directory. They copy stored runs before corrupting them, so they need those runs
($QUEUE_TOOLCHAIN/runs and data/raw/alamode_sto/); without them the test is skipped."""
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
Q = REPO / "queue"
TOOL = os.environ.get("QUEUE_TOOLCHAIN")
INPUTS = [Path(TOOL) / "runs/rta_on_k8_m12", Path(TOOL) / "runs/scph_on_k8"] if TOOL else []
INPUTS.append(REPO / "data/raw/alamode_sto/RTA_scph_300K.in")


@pytest.mark.skipif(not TOOL or not all(p.exists() for p in INPUTS),
                    reason="needs the stored 2026-10-07 runs and the raw ALAMODE dataset (local-only)")
def test_queue_unit_tests_pass():
    r = subprocess.run([sys.executable, "test_queue.py"], cwd=Q, capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    assert "FAIL" not in r.stdout
