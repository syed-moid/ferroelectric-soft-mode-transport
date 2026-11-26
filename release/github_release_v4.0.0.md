Scripts, tables and raw-output archives of the paper "Polar and rotational soft modes as heat carriers and scattering
partners in cubic SrTiO₃": two anharmonic models fitted to different force data (model A, the ALAMODE SrTiO₃ example;
model B, a fit to the Materials Cloud Archive 2026.91 force data), family labels from atomic motion, and the
three-phonon linewidth of the half-current carriers split by the labels of the two scattering partners.

- Added: `scripts/v4/` (regeneration of every v4 table and number from the raw outputs; figure scripts),
  `data/processed/v4/` (summaries, report tables, number table in parts a and b), `data/provenance/v4/` (run map and
  scrubbed manifest), `tools/alamode_patches/t3_partition.diff` (linewidth partition of `anphon`, with build notes),
  `fekappa.modes` and `fekappa.partition` with tests, `docs/v4/ANALYSIS_STEPS.md`.
- Files of the 2015 study. The production force constants, Born data and numerical results of the 2015 study (T.
  Tadano and S. Tsuneyuki, Phys. Rev. B 92, 054301 (2015)) were shared privately by T. Tadano and were used in the
  calculations of the paper (the benchmark and the fit test). They are not included in this release. The results
  obtained with them are given in the paper and its Supplementary Information. Scripts that take those files as
  arguments are included. Among them are the scripts for the archived 2015 fit.
- New Zenodo archives: `v4_runs_a`, `v4_runs_b`, `v4_fits_b`, `v4_processed_a`, `v4_processed_b`; the v3 archives are
  unchanged.
- The v3 folders are renamed `scripts/v3/` and `data/processed/v3/` (`docs/PATH_MAP_v3.md`).
- Files derived from the model B dataset are under CC BY-NC 4.0, in their own folders (`DATA_LICENSE.md`).
- Package version 4.0.0.
