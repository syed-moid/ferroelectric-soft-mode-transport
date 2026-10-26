"""Provenance headers for the two revision tables that no script generates.

- c1_secondary_sweep.csv: κ_P, κ_C and wall time collected from the rta.log files of the secondary settings sweep
  (300 K, 12³, base = own SCPH q₁ = 8³ with polarization mixing; one setting changed per run).
- zhao_values.csv: values transcribed from the full text of Zhao et al., Phys. Rev. B 104, 224304 (2021); the page and
  table or figure of each value are in the columns source_page and source_kind.
"""
from pathlib import Path

from fekappa.provenance import stamp

PJ = Path(__file__).resolve().parents[2] / "data/processed/revision"
stamp(PJ / "c1_secondary_sweep.csv", __file__,
      ["$FEKAPPA_RUNS/rta_on_k8_m12/rta.log and $FEKAPPA_RUNS/sweep_{ismear1_eps5,ismear1_eps10,ismear0_eps10,"
       "classical,na0}/rta.log (zenodo:prequeue_runs.tar.gz); values collected by hand"])
stamp(PJ / "zhao_values.csv", __file__,
      ["Zhao et al., Phys. Rev. B 104, 224304 (2021), full text; values transcribed by hand"])
print("stamped 2 tables")
