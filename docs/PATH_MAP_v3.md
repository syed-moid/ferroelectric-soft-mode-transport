# Paths renamed in v4.0.0

Version 4.0.0 renames two folders of v3 to neutral names. The files and their contents are unchanged except for
the path references inside them (provenance headers, imports, documentation), which name the new folders.

| v3.0.0 – v3.1.1 | v4.0.0 | files |
|---|---|---|
| `scripts/revision/` | `scripts/v3/` | 39 |
| `data/processed/revision/` | `data/processed/v3/` | 48 |

Not renamed:

- `data/raw/alamode_sto/revision/` (not tracked): the folder into which the raw-output archives of v3 unpack. The
  archives of the Zenodo record and the scripts that read the raw outputs keep this name.
- `data/provenance/manifest_final.json`: the manifest of the v3 runs. It is unchanged, so its paths and sha256
  values still describe the archives of v3.1.1, which use the old names.

The archives and the manifest of v3.1.1 therefore name the old paths; the table below maps each file.

| old path | new path |
|---|---|
| `data/processed/revision/a1_decision.csv` | `data/processed/v3/a1_decision.csv` |
| `data/processed/revision/a1_decision_T1000.csv` | `data/processed/v3/a1_decision_T1000.csv` |
| `data/processed/revision/a1_identity_counts.csv` | `data/processed/v3/a1_identity_counts.csv` |
| `data/processed/revision/a1_per_mesh.csv` | `data/processed/v3/a1_per_mesh.csv` |
| `data/processed/revision/a1_per_mesh_T1000.csv` | `data/processed/v3/a1_per_mesh_T1000.csv` |
| `data/processed/revision/a2_harmonic_xcheck.csv` | `data/processed/v3/a2_harmonic_xcheck.csv` |
| `data/processed/revision/c1_benchmark_matrix.csv` | `data/processed/v3/c1_benchmark_matrix.csv` |
| `data/processed/revision/c1_scph_softmodes.csv` | `data/processed/v3/c1_scph_softmodes.csv` |
| `data/processed/revision/c1_secondary_sweep.csv` | `data/processed/v3/c1_secondary_sweep.csv` |
| `data/processed/revision/c2_absw_bare.csv` | `data/processed/v3/c2_absw_bare.csv` |
| `data/processed/revision/c2_common_subset.csv` | `data/processed/v3/c2_common_subset.csv` |
| `data/processed/revision/c3_damping_bookkeeping.csv` | `data/processed/v3/c3_damping_bookkeeping.csv` |
| `data/processed/revision/d2_summary.csv` | `data/processed/v3/d2_summary.csv` |
| `data/processed/revision/d2_top_pairs.csv` | `data/processed/v3/d2_top_pairs.csv` |
| `data/processed/revision/d3_table_pm_mixing.csv` | `data/processed/v3/d3_table_pm_mixing.csv` |
| `data/processed/revision/d3_table_q1_sampling.csv` | `data/processed/v3/d3_table_q1_sampling.csv` |
| `data/processed/revision/d4_bare_framings.csv` | `data/processed/v3/d4_bare_framings.csv` |
| `data/processed/revision/d4_character_match_300K.csv` | `data/processed/v3/d4_character_match_300K.csv` |
| `data/processed/revision/d_q1_freq_shift.csv` | `data/processed/v3/d_q1_freq_shift.csv` |
| `data/processed/revision/dfc2_reproduction.csv` | `data/processed/v3/dfc2_reproduction.csv` |
| `data/processed/revision/e1_block_factors.csv` | `data/processed/v3/e1_block_factors.csv` |
| `data/processed/revision/e1_factor_freezing.csv` | `data/processed/v3/e1_factor_freezing.csv` |
| `data/processed/revision/e1_pair_counts.csv` | `data/processed/v3/e1_pair_counts.csv` |
| `data/processed/revision/e1_prod_block_factors.csv` | `data/processed/v3/e1_prod_block_factors.csv` |
| `data/processed/revision/e1_prod_factor_freezing.csv` | `data/processed/v3/e1_prod_factor_freezing.csv` |
| `data/processed/revision/e1_sector_sizes.csv` | `data/processed/v3/e1_sector_sizes.csv` |
| `data/processed/revision/e2_coherence_blocks_tracked.csv` | `data/processed/v3/e2_coherence_blocks_tracked.csv` |
| `data/processed/revision/e2_fineT_blocks.csv` | `data/processed/v3/e2_fineT_blocks.csv` |
| `data/processed/revision/e2_label_changes.csv` | `data/processed/v3/e2_label_changes.csv` |
| `data/processed/revision/eta4_margin.csv` | `data/processed/v3/eta4_margin.csv` |
| `data/processed/revision/eta4_margin_cases.csv` | `data/processed/v3/eta4_margin_cases.csv` |
| `data/processed/revision/eta_max_outside_map.csv` | `data/processed/v3/eta_max_outside_map.csv` |
| `data/processed/revision/exp_300K_crossread.csv` | `data/processed/v3/exp_300K_crossread.csv` |
| `data/processed/revision/exp_datasets.csv` | `data/processed/v3/exp_datasets.csv` |
| `data/processed/revision/exp_gap_vs_T.csv` | `data/processed/v3/exp_gap_vs_T.csv` |
| `data/processed/revision/mfp_vs_eta.csv` | `data/processed/v3/mfp_vs_eta.csv` |
| `data/processed/revision/n3_q1_16_check.csv` | `data/processed/v3/n3_q1_16_check.csv` |
| `data/processed/revision/n4_results.csv` | `data/processed/v3/n4_results.csv` |
| `data/processed/revision/n4_selection.csv` | `data/processed/v3/n4_selection.csv` |
| `data/processed/revision/prior_work_comparison.md` | `data/processed/v3/prior_work_comparison.md` |
| `data/processed/revision/r_point_check.csv` | `data/processed/v3/r_point_check.csv` |
| `data/processed/revision/s21_tolerance_printed.txt` | `data/processed/v3/s21_tolerance_printed.txt` |
| `data/processed/revision/s7_ratios_unrounded.csv` | `data/processed/v3/s7_ratios_unrounded.csv` |
| `data/processed/revision/sector_degeneracy.csv` | `data/processed/v3/sector_degeneracy.csv` |
| `data/processed/revision/v1_removed_vs_mesh.csv` | `data/processed/v3/v1_removed_vs_mesh.csv` |
| `data/processed/revision/v1_validation.csv` | `data/processed/v3/v1_validation.csv` |
| `data/processed/revision/wang_comparisons.csv` | `data/processed/v3/wang_comparisons.csv` |
| `data/processed/revision/zhao_values.csv` | `data/processed/v3/zhao_values.csv` |
| `scripts/revision/README.md` | `scripts/v3/README.md` |
| `scripts/revision/a1_decision.py` | `scripts/v3/a1_decision.py` |
| `scripts/revision/a1_identity_count.py` | `scripts/v3/a1_identity_count.py` |
| `scripts/revision/a1_identity_count_table.py` | `scripts/v3/a1_identity_count_table.py` |
| `scripts/revision/a2_harmonic_xcheck.py` | `scripts/v3/a2_harmonic_xcheck.py` |
| `scripts/revision/build_manuscript_numbers_v8.py` | `scripts/v3/build_manuscript_numbers_v8.py` |
| `scripts/revision/build_si_tables_v4.py` | `scripts/v3/build_si_tables_v4.py` |
| `scripts/revision/c1_collect.py` | `scripts/v3/c1_collect.py` |
| `scripts/revision/c1_scph_softmodes.py` | `scripts/v3/c1_scph_softmodes.py` |
| `scripts/revision/c2_common_subset.py` | `scripts/v3/c2_common_subset.py` |
| `scripts/revision/c3_damping_bookkeeping.py` | `scripts/v3/c3_damping_bookkeeping.py` |
| `scripts/revision/d2_coherence_report.py` | `scripts/v3/d2_coherence_report.py` |
| `scripts/revision/d3_sampling_isolation.py` | `scripts/v3/d3_sampling_isolation.py` |
| `scripts/revision/d4_bare_subset.py` | `scripts/v3/d4_bare_subset.py` |
| `scripts/revision/d_q1_freq_shift.py` | `scripts/v3/d_q1_freq_shift.py` |
| `scripts/revision/dfc2_reproduction_check.py` | `scripts/v3/dfc2_reproduction_check.py` |
| `scripts/revision/e1_factor_decomposition.py` | `scripts/v3/e1_factor_decomposition.py` |
| `scripts/revision/e1_pair_counts.py` | `scripts/v3/e1_pair_counts.py` |
| `scripts/revision/e1_production.py` | `scripts/v3/e1_production.py` |
| `scripts/revision/e2_branch_tracking.py` | `scripts/v3/e2_branch_tracking.py` |
| `scripts/revision/e2_fineT_blocks.py` | `scripts/v3/e2_fineT_blocks.py` |
| `scripts/revision/eta4_margin.py` | `scripts/v3/eta4_margin.py` |
| `scripts/revision/eta_max_outside_map.py` | `scripts/v3/eta_max_outside_map.py` |
| `scripts/revision/exp_300K_crossread.py` | `scripts/v3/exp_300K_crossread.py` |
| `scripts/revision/exp_datasets.py` | `scripts/v3/exp_datasets.py` |
| `scripts/revision/export_manifest.py` | `scripts/v3/export_manifest.py` |
| `scripts/revision/mfp_vs_eta.py` | `scripts/v3/mfp_vs_eta.py` |
| `scripts/revision/n3_q1_16_check.py` | `scripts/v3/n3_q1_16_check.py` |
| `scripts/revision/n4_results_table.py` | `scripts/v3/n4_results_table.py` |
| `scripts/revision/n4_select.py` | `scripts/v3/n4_select.py` |
| `scripts/revision/prod_analysis.py` | `scripts/v3/prod_analysis.py` |
| `scripts/revision/prod_kappa_table.py` | `scripts/v3/prod_kappa_table.py` |
| `scripts/revision/r_point_check.py` | `scripts/v3/r_point_check.py` |
| `scripts/revision/s21_printed_outputs.py` | `scripts/v3/s21_printed_outputs.py` |
| `scripts/revision/s7_ratios_single_rounding.py` | `scripts/v3/s7_ratios_single_rounding.py` |
| `scripts/revision/sector_degeneracy_check.py` | `scripts/v3/sector_degeneracy_check.py` |
| `scripts/revision/stamp_transcribed_tables.py` | `scripts/v3/stamp_transcribed_tables.py` |
| `scripts/revision/v1_validation_table.py` | `scripts/v3/v1_validation_table.py` |
| `scripts/revision/wang_comparisons.py` | `scripts/v3/wang_comparisons.py` |
