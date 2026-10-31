# Recorded rules and predictions

Each rule or prediction below was written down at the stated time (UTC), relative to the runs it governs. The rule
text is quoted as recorded. Job ids are those of `data/provenance/manifest_final.json`; the tables hold the numbers.
Widths Γ are half-widths (HWHM) in cm⁻¹.

| id | subject | written (UTC) | relative to the runs | outcome | reported in |
|---|---|---|---|---|---|
| R1 | SCPH grid q₁: 8³ or 12³ | before 2026-10-09T21:31Z | before | criteria not met; q₁ = 12³ adopted | SI S1, S4 (Table S7) |
| R2 | κ<sub>C</sub> at 24³, 300 K | 2026-10-14T14:05Z | before | prediction failed | SI S3 (Table S5, Fig. S3) |
| R3 | acceptance tests of the four-phonon exclusion build | 2026-10-15T17:25Z | after the first 4³ results | both replacement tests passed | SI S7.10 (Table S22, Fig. S6) |
| R4 | four-phonon rule, 300 K mesh sequence | 2026-10-15T17:25Z | before | all three soft modes in the first class | SI S7.11 (Tables S24, S25) |
| R5 | four-phonon rule, 1000 K | 2026-10-16T16:16Z | before | all five modes in the first class | SI S7.11 (Tables S24, S25) |
| R6 | heat-carrier sample | 2026-10-16T16:16Z | before | descriptive (no pass/fail) | SI S7.12 (Table S26) |
| R7 | identity-tolerance sensitivity | 2026-10-16 | before | all changes ≤ 1.9 %; no further run triggered | SI S7.9 (Table S21) |
| R8 | κ<sub>C</sub> at 20³, 1000 K | 2026-10-18T23:04:44Z | before | prediction held | SI S3 (Table S5, Fig. S3) |
| R9 | additive form of the four-phonon margin | 2026-10-18 | after the results of R4–R6 | additive η₄ in the text; ratio form in the SI | §6.2; SI S7.11 |
| R10 | labelling of exactly degenerate subspaces | 2026-10-22 | after the production labels | whole-subspace labels adopted | §2.4; SI S5 (Table S8) |
| R11 | implementation checks of the exclusion | 2026-10-21/22 | after the results of R4–R7 | all checks passed | §6.2; SI S7.7 (Tables S18, S19) |

R1, R2 and R4–R8 were written before their runs. R3 and R9–R11 were written after the results named in the table.

## R1. SCPH grid q₁

**Runs:** `scph_pm_q8`, `scph_pm_q12`; `rta_pm_q8_T250_m12`, `rta_pm_q8_T300_m12`, `rta_pm_q12_T250_m12`,
`rta_pm_q12_T300_m12`. The rule was applied by the queue job `d0_decision` at 2026-10-10T01:29:58Z. The q₁ = 8³ SCPH
run and two of these RTA runs were made before the queue started and were imported into it.

> Adopt q1 = 8³ if, at both 250 K and 300 K, going 8³ → 12³ moves κ<sub>P</sub> by < 3% AND the R and M frequencies by < 5%;
> otherwise adopt q1 = 12³. Transport mesh fixed at 12³; PM on; tetrahedron; NONANALYTIC = 3; KAPPA_COHERENT = 2.

**Outcome.** Criteria not met at either temperature; q₁ = 12³ adopted.

| T (K) | q₁ | κ<sub>P</sub> | κ<sub>C</sub> | Ω Γ-TO1 | Ω M | Ω R (cm⁻¹) |
|---|---|---|---|---|---|---|
| 250 | 8³ | 10.1181 | 0.2589 | 128.24 | 78.28 | 20.12 |
| 250 | 12³ | 9.7065 | 0.2736 | 127.83 | 77.17 | 14.65 |
| 250 | Δ (%) | −4.07 | +5.68 | −0.32 | −1.42 | −27.22 |
| 300 | 8³ | 9.2384 | 0.2790 | 134.04 | 80.82 | 28.93 |
| 300 | 12³ | 9.1112 | 0.2772 | 133.69 | 79.86 | 25.82 |
| 300 | Δ (%) | −1.38 | −0.65 | −0.26 | −1.19 | −10.73 |

Tables: `data/processed/revision/d_q1_freq_shift.csv`, `data/processed/si_v4/s1_scph_convergence.csv`.

## R2. κ<sub>C</sub> at 24³, 300 K

**Runs:** `rta_pm_q12_T300_m24` (started 2026-10-14T20:27Z).

> At 24³ the κ<sub>C</sub>-weighted median r is ≈ 0.7. |Δκ<sub>C</sub>(20³ → 24³)| / κ<sub>C</sub>(20³) is smaller than the 300 K 16³ → 20³
> change (|−13.1%|) and of the order of the 500 K 16³ → 20³ change (−7.2%), i.e. 3–11% in magnitude.
> Outcome rule: "holds" if |Δ| < 13.1%; "holds as stated" if 3% ≤ |Δ| ≤ 11%; "fails" if |Δ| ≥ 13.1%. The sign is not
> predicted.

**Outcome.** Failed.
- Δκ<sub>C</sub>(20³ → 24³) = (0.5411 − 0.4417)/0.4417 = +22.5 %.
- The median-r part held: r = 0.70 at 24³.
- The 300 K sequence 0.508 (16³) → 0.442 (20³) → 0.541 (24³) is not monotonic.
- The prediction compared r of the finer mesh at 300 K with r of the coarser mesh at 500 K. It fails under either
  convention: with r of the coarser mesh the 300 K 20³ → 24³ pair has r = 0.59. Later comparisons use r of the coarser
  mesh of each pair.

Tables: `data/processed/revision/d2_summary.csv`, `data/processed/sto_production/kappa_mesh_T.csv`.

## R3. Acceptance tests of the four-phonon exclusion build

Two tests were written before any run of the build: (1) Γ₄ with the exclusion reproduces the fitted intercept B of
Γ₄ = a/ε + B within the fit's 0.8 %; (2) Γ₄ with the exclusion is independent of ε. On 2026-10-15T17:25Z, after the
4³ results of both tests were known, test (1) was replaced:

> Gate 1 as written ("patched Γ₄ reproduces the fitted B within the fit's 0.8%") is replaced by two direct tests:
> (1) the independent count of the excluded set matches the patched code exactly (all three modes); (2) the amount
> removed is 99.7–102.5% of the fitted a/ε. Gate 2 (ε-independence) is accepted at the decision-relevant factor-2
> precision, not at 10%.

Reason recorded with the change: the fitted B combines the identity term with the other 1/ε-like terms of the 4³
mesh, so it cannot be matched to 0.8 %.

**Results (4³, 300 K).**
- Test (1) as first written failed for all three modes.
- Test (2): passed for Γ-TO1 at 10 % and failed for R (22 %) and M (33 %); passed for all three at a factor of 2.
- Excluded counts, code against the independent count: 11268 / 11268 (Γ-TO1), 11268 / 11268 (R), 3772 / 3772 (M).
- Removed amount: 99.7 % (R) to 102–103 % (Γ-TO1, M) of the fitted a/ε.

Tables: `data/processed/si_v4/s7_gate_4cubed.csv`, `data/processed/revision/a1_identity_counts.csv`.

## R4. Four-phonon rule, 300 K mesh sequence

**Runs:** `m4ph_audit_q12_T300_m{4,5,6,8}_e*` and `evec_pm_q12_T300_m{4,6,8}` (started 2026-10-15T17:39Z), with the
rule build `alamode_4ph_audit` (`tools/alamode_patches/BUILD.md`).

> Γ₄,thr = the Γ₄ that makes η₃₊₄ = 0.5 at 300 K, from production Γ₃ (16³) and Ω: Γ-TO1 64.6, M 36.5, R 8.7.
> Stability, both required: S1: Γ₄(8³, ε = 5) is within a factor of 2 of Γ₄(6³, ε = 5); S2: the 6³ ε-scan
> (ε = 2, 5, 10) spans no more than a factor of 2.
> Outcome, per mode passing stability, with Γ₄,max the largest Γ₄ over {6³ at ε = 2, 5, 10; 8³ at ε = 5}:
> Γ₄,max ≤ ½ Γ₄,thr → first class; ½ Γ₄,thr < Γ₄,max < Γ₄,thr → second class; Γ₄,max ≥ Γ₄,thr → third class.
> A mode failing S1 or S2: "Γ₄ not determined at feasible meshes". The 4³ and 5³ jobs are reported but do not enter
> the rule. R-mode Γ₄ and η are reported as "in this model".

**Outcome.** All three modes passed S1 and S2 and fell in the first class.

| mode | Γ₄(6³; ε = 2 / 5 / 10) | Γ₄(8³, ε = 5) | S1 | S2 | Γ₄,max / Γ₄,thr |
|---|---|---|---|---|---|
| Γ-TO1 | 2.587 / 2.504 / 2.466 | 2.310 | 0.92 | 1.05 | 0.040 |
| M | 1.754 / 2.100 / 2.335 | 1.662 | 0.79 | 1.33 | 0.064 |
| R | 1.692 / 1.979 / 2.092 | 2.111 | 1.07 | 1.24 | 0.243 (in this model) |

Tables: `data/processed/revision/a1_decision.csv`, `a1_per_mesh.csv`.

## R5. Four-phonon rule, 1000 K

**Runs:** `m4ph_audit_q12_T1000_m{6,8}_e*` and `evec_pm_q12_T1000_m{6,8}` (started 2026-10-16T16:17Z).
- Descriptive and outside the rule: `m4ph_auditdesc_q12_T1000_m8_e*`.

> The R4 rule applied unchanged at 1000 K, with Γ₄,thr = ½Ω − Γ₃ from production Γ₃ (16³) and Ω at 1000 K:
> Γ-TO1 85.6, M 53.2, R 39.1, X branch 7: 108.7, R branch 13: 232.1.

**Outcome.** All five modes passed S1 and S2 and fell in the first class. Γ₄,max / Γ₄,thr = 0.052 (Γ-TO1), 0.110
(M), 0.143 (R, in this model), 0.083 (X b7), 0.034 (R b13). Tables: `data/processed/revision/a1_decision_T1000.csv`,
`a1_per_mesh_T1000.csv`.

## R6. Heat-carrier sample

**Runs:** `m4ph_n4_q12_T{300,800}_m{6,8}_e*` (started 2026-10-17T13:26Z).

> Descriptive; no pass/fail; no κ<sub>P</sub> recomputation from these modes. Two independent samples: q on the 8³ mesh,
> ranked from the 16³ table, and q on the 6³ mesh, ranked from the 12³ table; per-mode κ<sub>P</sub> contribution at 300 K;
> classes by the character rule; the 2 largest per class, at most one mode per q-point, ties by lower q index then
> branch. The modes are run unchanged at 800 K.

**Outcome.**
- Γ₄/Γ₃ = 0.12–1.08 at 300 K and 0.15–1.11 at 800 K.
- Lifetime factor 1/(1 + Γ₄/Γ₃) = 0.48–0.90 at 300 K and 0.47–0.87 at 800 K.
- The sample carries 1.2–1.3 % of κ<sub>P</sub>.

Tables: `data/processed/revision/n4_selection.csv`, `n4_results.csv`.

## R7. Identity-tolerance sensitivity

**Runs:** `m4ph_audittol6_q12_T300_m6_e{2,5}` with the 1e-6 Ry build `alamode_4ph_audit_tol1e6`.

> Rerun 6³ at ε = 2 and 5 for the three modes with the identity-frequency tolerance at 1e-6 Ry instead of 1e-9 Ry.
> The outcome of R4 stays the 1e-9 result; the 1e-6 values are a sensitivity check. If the 6³ difference exceeds 10%
> of Γ₄ for any mode, also run 8³ ε = 5 with the 1e-6 build.

**Outcome.** The largest change is 1.9 % (R, ε = 2), and every change is a decrease. The 10 % condition was not
reached, so no 8³ run followed. Tables: `data/processed/si_v4/s7_tolerance.csv`; printed outputs in
`data/processed/revision/s21_tolerance_printed.txt`.

## R8. κ<sub>C</sub> at 20³, 1000 K

**Runs:** `rta_pm_q12_T1000_m20`, `evec_pm_q12_T1000_m20` (started 2026-10-18T23:05Z).

> |Δκ<sub>C</sub>(16³→20³)| ≤ 2.5 % at 1000 K, with κ<sub>C</sub>(16³, 1000 K) = 0.51473 W m⁻¹ K⁻¹. Holds if
> 0.50186 ≤ κ<sub>C</sub>(20³) ≤ 0.52760; fails otherwise. Basis: r(16³, 1000 K) = 1.13, above the r ≥ 0.88 range where
> |Δκ<sub>C</sub>(16³→20³)| ≤ 2.2%.

**Outcome.** Held: κ<sub>C</sub>(20³, 1000 K) = 0.50306 W m⁻¹ K⁻¹, Δκ<sub>C</sub> = −2.27 % (pair sums). Tables:
`data/processed/revision/d2_summary.csv`, `data/processed/sto_production/kappa_mesh_T.csv`.

## R9. Additive form of the four-phonon margin

Written after the results of R4–R6. The ratio Γ₄/Γ₃ reached 2.91 for Γ-TO1 at 1000 K and 1.11 in the heat-carrier
sample, so a statement based on "Γ₄/Γ₃ ≤ 1 in every sampled mode" no longer applied.
- The text uses the additive form η₃₊₄ = η₃ + η₄, with η₄ = Γ₄/Ω.
- The ratio form is kept as an SI table.

No pass/fail test is attached. Over the 71 mesh/ε cases the largest η₄ is 0.0818 (R, 300 K, 8³, ε = 5). Tables:
`data/processed/revision/eta4_margin.csv`, `eta4_margin_cases.csv`.

## R10. Labelling of exactly degenerate subspaces

Written after the production labels were known. Inside an exactly degenerate subspace the per-eigenvector label
depends on the basis returned by the diagonalizer.
- Subspaces with |ΔΩ| < 10⁻³ cm⁻¹ are labelled as a whole, by their summed projections (`sector_char`).
- The per-eigenvector label is kept as `sector_char_eig`.

The change affects 6 % of such subspaces, carrying at most 0.16 % of κ<sub>P</sub>, and moves the hard/hard share of κ<sub>C</sub> by
at most 0.6 percentage points. Table: `data/processed/revision/sector_degeneracy.csv`.

## R11. Implementation checks of the exclusion

Written after the results of R4–R7. These checks verify the implementation; they test no prediction.

**Runs:** `m4ph_trace_q12_T300_m4_e5`, `m4ph_v1omit_q12_T300_m6_e5`, `m4ph_pm_q12_T300_m6_e5`,
`m4ph_v2strict_q12_T300_m{6,8}_e*`, `m4ph_v3dump_q12_T300_m6_e5` and `m4ph_v3dump_q12_T300_m6_e5_M`.

**Criterion, set before the strict-index runs finished:**
- the omitted-sum check passes;
- the enlarged − strict difference is itself a finite-mesh term (it falls with N_k and goes as 1/ε);
- the sign convention is consistent.

**Outcome.** All three conditions were met.
- The omitted sum reproduces unmodified − patched to printed precision.
- ε × (enlarged − strict) is constant to 3 % at 6³, and the difference falls by a factor of 2.7–5.5 from 6³ to 8³.
- The trace of one quartet agrees with the code to ≤ 2 × 10⁻¹⁵ relative.
- The full-precision dump explains the 6³ count difference of R4 as the 7-digit precision of the eigenvalue file.

Tables: `data/processed/revision/v1_validation.csv`, `v1_removed_vs_mesh.csv`.
