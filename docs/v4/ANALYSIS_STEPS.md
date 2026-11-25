# Analysis steps of v4

<!-- DRAFT: a map of the steps, their scripts and runs. No interpretation. Steps marked [HELD]
need the files of the 2015 study. Run names: data/provenance/v4/run_map.csv. -->

Each step names what it computes, the script in `scripts/v4/`, the runs it reads and the outputs in
`data/processed/v4/` (a = model A, b = model B). Rules marked "fixed before" were written down before the
calculation they govern; the dated records are part of the author's working notes, not of this package.

| step | what it computes | script(s) | runs | outputs |
|---|---|---|---|---|
| B1 | fits of model B (H1: harmonic terms fixed; H2: all orders) and force errors on the training and test sets; model A on B's test set | `b0_inventory.py`, `b1_convert_dfset.py`, `b1_h1_fc2_conversion.py`, `b1_h1_fc2_verify.py`, `b1_force_eval.py`, `figures/make_fig1_parity.py` | fits (ALM) | Table 1, Table S7, Fig. 1 |
| B3 | SCPH soft modes of H1 and H2 on q₁ = 2³, 4³, 8³, 12³; comparison with the source study's bands | `scph_softmodes.py`, `b3_source_compare.py` | `B_scph_q*`, `H2_scph_q*` | `b/b3_scph_softmodes.csv` |
| T1 | transport of model B at 300 K (κ<sub>P</sub>, κ<sub>C</sub>, R mode, coherence blocks); A next to B | `a3_transport_analysis.py`, `b4_model_comparison.py`, `kp_eta_distribution.py` | `B_rta_q12_m16_T300`, `B_rta_q8_m12_T300`, `B_evec_q12_m16_T300`, `A_*` | `b/t1/`, `b/b4_model_comparison.csv`, `{a,b}/kp_eta_300K_m16*` |
| A1–A4 | model A: regeneration of the distributed fit, penalty scan (accepted interval, endpoints a198 and a600), reduced-basis refit (1,053 parameters), image-convention test | `a1b_compare_fc.py`, `a2_compare.py`, `a3_extract.py`, `a3_transport_analysis.py` | `A_bench_*`, `A198_*`, `A600_*` | `a/benchmark/`, `a/t1/charphase_*` |
| T2 | family census on the full mesh (polar score weighted by cos²(πq) per axis; rotational score by sin²(πq_b) sin²(πq_c)); frozen labels (four definitions) | `t2_family_census.py`, `t2_freeze_labels.py`, `t2_census_tables.py` | `*_evec_*`, `*_harm_*` | `{a,b}/labels/`, `{a,b}/census/` (archives) |
| T3 | patched `anphon`: per carrier, the three-phonon linewidth binned by the labels of the two partners; checks (Σ bins = total; total = `.result`; outputs byte-identical to the base build; label invariance) | patch `tools/alamode_patches/t3_partition.diff`; `t3_reconstruct.py`, `t3_symcheck.py`, `t3_label_invariance.py`, `t3_neutral_labels.py`, `t4_encoding_check.py` | `*_partition_*` | checks |
| T4 | carriers (whole degenerate groups ranked by κ<sub>P</sub> until ≥ 50 %, fixed before), their linewidth by number of soft partners and by partner-family pair, κ<sub>P</sub> by family, modes with η > 0.1; 8³ sensitivity to TRISYM | `t4_analysis.py`, `t4_trisym_sensitivity.py`, `t4_report.py`, `figures/make_t4_figures.py` | `{A,B}_partition_q12_m16_T300`, `*_m8_T300_trisym{1,0}` | `{a,b}/t4/`, Table 5, Fig. 4, Fig. S3 |
| T5 | partners within frequency bands (edges 105 / 145 / 160 cm⁻¹ from the terciles and maximum of model A's rotational frequencies, fixed before): share per 1,000 available modes, rotational against non-rotational | `t5_control_labels.py`, `t5_control_analysis.py`, `t5_report.py` | `*_partition_q12_m16_T300_bands` | `{a,b}/t5/`, Table 6, Fig. 5 |
| T5b | the lowest band divided at 80 and 95 cm⁻¹ (fixed before); mean frequencies; regions of the zone of the partners below 80 cm⁻¹ | same, `v2_nonrot_below80.py`, `rot_below_105.py` | `*_partition_q12_m16_T300_lowbands` | `{a,b}/t5b/`, Table 7, Fig. S2 |
| T6 | the T4 rows on 20³ (a share is mesh-stable if 16³ and 20³ differ by ≤ 2 percentage points, fixed before) | `t6_report.py` | `*_partition_q12_m20_T300` | `{a,b}/t4/m20/`, Tables 5, S8–S12 |
| T7 | the T4 rows at 600 K (labels by the same rule at 600 K) | `t7_report.py` | `*_partition_q12_m16_T600` | `{a,b}/t4/T600/`, Tables S13–S17 |
| T8 | full sum (TRISYM = 0) against the symmetry-reduced sum on 16³ | `t8_widths.py`, `t8_report.py` | `*_partition_q12_m16_T300_fullsum` | `{a,b}/t4/trisym0/`, Tables S30–S32 |
| P0 [HELD] | checks of the files of the 2015 study against the distributed example | `p0_archive_files.py` | — | held |
| P1 [HELD] | the archived fit in the steps from model A to the configuration of the 2015 study (fit, non-analytic term, Born data), q₁ 2³–12³ | `make_inputs_archived_fit.py`, `scph_softmodes.py` | `P_scph_*` | Table 2, Fig. 2 |
| P2 [HELD] | κ at 300 K in those steps against the 2015 value; mode comparison with the 2015 per-mode file (positional pairing) | `p2_analysis.py` | `P_rta_*_q8_m12_T300` | Table 2 |
| P3 [HELD] | the archived fit under the settings of model A: T2–T5b repeated; comparison with model A by limits fixed before (± 3 percentage points; ratio classes ≥ 2.0, 1.3–2.0, < 1.3) | `p3_sampling.py`, `p3_report.py` | `P_*_settingsA`, `P_partition_*` | Table 8, Tables S33–S38 |
| P4 [HELD] | forces of the archived fit on model A's 40 configurations and model B's test set | `b1_force_eval.py --xml` | — | Table 1, Table S7 |
| P5 [HELD] | κ(T) in the configuration of the 2015 study, 200–1000 K | `p5_kappa_T.py`, `figures/make_kappa_T_2015.py` | `P_rta_e_q8_m12_T*` | Table S39, Fig. S4 |
| P6 [HELD] | where the difference at 300 K sits: the 2015 per-mode lifetimes (and heat capacities at the 2015 frequencies) in the present conductivity sum, present group velocities | `p6_kappa_split.py` | `P_rta_e_q8_m12_T300` | Sec. 3.2, SI Sec. S12.1, Table S40 |
| GA | graphical abstract: κ<sub>P</sub> by family, linewidth with at least one / no soft partner, per-mode ratios by band, read from the T4, T5 and T5b outputs; `--without-archived` for models A and B only; [HELD] with the archived fit | `figures/make_graphical_abstract_v3.py` | — | graphical abstract |

The driver `regenerate_v4.py` runs B3, T1, T2, T4–T8 and the counts in this order from the stored raw outputs;
`regenerate_v4_held.py` runs P1–P6 [HELD]; `build_tables_v4.py` and `build_manuscript_numbers_v4.py` write the tables
and the number table.
