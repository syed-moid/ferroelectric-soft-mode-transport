# Revision analysis scripts

These scripts turn the stored solver outputs into the tables behind every number, figure and SI table of the revised
manuscript. Run them from the repository root with the package installed (`pip install -e .`).

- **Outputs:** tables go to `data/processed/revision/`, `data/processed/sto_production/` and
  `data/processed/si_v4/`. Each table in `revision/` and `si_v4/` starts with a provenance header
  (`fekappa.provenance`) naming its script, its inputs, the date, the package version and the hash of the job
  manifest.
- **Inputs:** the scrubbed job manifest `data/provenance/manifest_final.json`; tables already in the repository; and
  three kinds of input kept outside git:
  - `data/raw/alamode_sto/revision/`: the raw outputs of every queue job. They are in the Zenodo record; unpack the
    archives there.
  - `$FEKAPPA_RUNS`: the runs made before the queue existed (Zenodo `prequeue_runs.tar.gz`).
  - `$FEKAPPA_TOOLCHAIN` (`sto_example/` with the force-constant inputs) and `$FEKAPPA_DFPT` (the root of the
    lattice-shear-viscosity repository, for the harmonic cross-check only).
- **From tracked tables alone:** the manuscript numbers (`build_manuscript_numbers_v8.py`) and the figures
  (`paper/figures/make_v8_figures.py`, `make_si_v4_figures.py`) regenerate without any of these.

Order: `prod_kappa_table` → `prod_analysis` → the scripts below → `build_si_tables_v4` → `build_manuscript_numbers_v8`.

| script | inputs | output (data/processed/…) | used in |
|---|---|---|---|
| `export_manifest.py` | live queue manifest (not tracked) | `../provenance/manifest_final.json` | provenance of every job |
| `prod_kappa_table.py` | manifest | `sto_production/kappa_mesh_T.csv`, `kappa_spread_T.csv` | Fig. 1, 3, 4; §3.2; SI S3 |
| `prod_analysis.py` | manifest; raw `pm_q12` RTA and eigenvector outputs | `sto_production/mode_table_T.csv.gz`, `eta_validity_T.csv`, `kappa_sectors_T.csv`, `coherence_blocks_T.csv` | Figs. 2–4; §4, §5; SI S5, S6 |
| `d2_coherence_report.py` | raw `.kc_elem`; mode table | `revision/d2_summary.csv`, `d2_top_pairs.csv`; `sto_production/convergence_report_kC.csv` | §3.2 (r diagnostic); SI S3 |
| `sector_degeneracy_check.py` | mode table; raw eigenvectors | `revision/sector_degeneracy.csv` | §2.4; SI S5 (Table S8) |
| `e1_production.py` | raw `.kc_elem` (16³); mode table | `revision/e1_prod_block_factors.csv`, `e1_prod_factor_freezing.csv` | §5.2; SI S5 (Tables S12, S14) |
| `e1_pair_counts.py` | `e1_prod_block_factors.csv`; mode table | `revision/e1_pair_counts.csv` (per record and star-weighted), `e1_sector_sizes.csv` | §5.2; SI S5 (Tables S12, S13) |
| `e1_factor_decomposition.py` | distributed example outputs | `revision/e1_block_factors.csv`, `e1_factor_freezing.csv` | supporting analysis of the distributed dataset; not cited in the revised text |
| `e2_branch_tracking.py`, `e2_fineT_blocks.py` | distributed example outputs; `$FEKAPPA_RUNS` | `revision/e2_*.csv` | supporting analysis (branch tracking, fine T steps); not cited in the revised text |
| `c1_scph_softmodes.py`, `c1_collect.py`, `d3_sampling_isolation.py` | `$FEKAPPA_RUNS` SCPH and RTA matrix | `revision/c1_*.csv`, `d3_*.csv` | §3.1, Table 1; SI S2 (`d3_*` supporting only) |
| `d_q1_freq_shift.py`, `r_point_check.py`, `n3_q1_16_check.py` | SCPH runs (q₁ = 8³, 12³, 16³) | `revision/d_q1_freq_shift.csv`, `r_point_check.csv`, `n3_q1_16_check.csv` (incl. 500 and 800 K at q₁ = 16³) | §3.1, §3.3; SI S1, S4 (Table S7), S7b |
| `a2_harmonic_xcheck.py` | `$FEKAPPA_RUNS/harmonic_xcheck`; `$FEKAPPA_DFPT` | `revision/a2_harmonic_xcheck.csv` | SI S1 |
| `c2_common_subset.py`, `c3_damping_bookkeeping.py`, `d4_bare_subset.py` | distributed example outputs; raw 8³ outputs; `$FEKAPPA_RUNS` | `revision/c2_*.csv`, `c3_*.csv`, `d4_*.csv` | §6.4; SI S9 (Table S32, from `d4_bare_framings.csv`); the others supporting only |
| `a1_decision.py` (and `A1_T=1000`) | manifest; mode table | `revision/a1_per_mesh*.csv`, `a1_decision*.csv` | §6.2; SI S7 (R4, R5) |
| `a1_identity_count.py`, `a1_identity_count_table.py` | raw eigenvector meshes | `revision/a1_identity_counts.csv` | SI S7.8 |
| `n4_select.py`, `n4_results_table.py` | mode table; manifest | `revision/n4_selection.csv`, `n4_results.csv` | SI S7.12 (R6) |
| `eta4_margin.py` | manifest; `eta_validity_T.csv`; `n4_results.csv` | `revision/eta4_margin*.csv` | §6.2; SI S7.11 (R9) |
| `v1_validation_table.py` | manifest | `revision/v1_validation.csv`, `v1_removed_vs_mesh.csv` | §6.2; SI S7.7 (R11) |
| `s21_printed_outputs.py` | raw tolerance runs | `revision/s21_tolerance_printed.txt` (full-precision inputs of Table S21) | SI S7.9 (R7) |
| `eta_max_outside_map.py` | raw 8³ and 24³ RTA outputs | `revision/eta_max_outside_map.csv` | §4 ¶1 (8³/24³ maxima 0.154/0.164) |
| `mfp_vs_eta.py`, `wang_comparisons.py` | mode table; manifest | `revision/mfp_vs_eta.csv`, `wang_comparisons.csv` | §6.3; SI S8 |
| `exp_datasets.py`, `exp_300K_crossread.py` | `data/literature/STO/`; figure readings in the scripts | `revision/exp_*.csv` | Fig. 3; §5.1 |
| `stamp_transcribed_tables.py` | — | headers of `c1_secondary_sweep.csv` and `zhao_values.csv` (values collected or transcribed by hand) | SI S2, S10 |
| `build_si_tables_v4.py` | all of the above; `revision/prior_work_comparison.md` | `si_v4/*.csv`, `si_v4/si_numbers_v4.csv`, `si_v4/build_manifest.csv` | every SI table; the benchmark percentages of §3.1 and Table S27 |
| `build_manuscript_numbers_v8.py` | tracked tables only | `manuscript_numbers_v8.csv` | every number in the main text |

`data/processed/revision/prior_work_comparison.md` is the comparison of prior SrTiO₃ studies read from their full
texts (SI Sec. S10). The rules R1–R11 cited in the docstrings are in `docs/recorded_rules.md`.
