# Archived: externally calibrated scalar SCPH (rejected calibration)

This folder preserves an earlier scalar self-consistent-phonon (SCPH) calibration that
was **rejected as a validation target**. It is reproduced by
`notebooks/09_BTO_hybrid_appendix.ipynb` as a documented negative result and is otherwise
not used by the active package.

The scalar SCPH was calibrated from two literature anchors (Ehsan *et al.* 2021 two-anchor
scheme, giving `A = 208.28 meV²`). The arithmetic is correct — the self-consistency
iteration converges across the 100–650 K grid — but the BaTiO₃ dataset it was compared
against is a hybrid, Cochran-calibrated phenomenological reconstruction, **not**
first-principles renormalized phonons. The SCPH-vs-data residual therefore diverges
(`|δ| > 2800%` near the transition). The lesson preserved here is that a phenomenologically
reconstructed dataset cannot validate a first-principles self-energy framework; this
motivates the reduced soft-mode kernel developed in the manuscript.

## Contents

- `scph_calibration_snapshot.py` — the scalar-SCPH calibration code
  (`calibrate_A_from_ehsan`, `scph_iterate`, `scph_curve`, `scph_residual`,
  `scph_tc_estimate`, and the `_bose_factor` helper), kept as a standalone snapshot so the
  negative result is fully reproducible. The active package under `src/fekappa/` does not
  import or depend on this module.
