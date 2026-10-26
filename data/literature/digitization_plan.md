# Literature digitization plan — v1

For each source below, the steps are:

1. Open the PDF in `data/<filename>.pdf`
2. Navigate to the figure(s) listed
3. Use [WebPlotDigitizer](https://apps.automeris.io/wpd/) to extract data points
4. Save extracted points to the CSV path listed, with the schema given
5. After all CSVs are filled, run `pytest tests/test_literature_inputs.py`

If a figure does not contain the requested data, leave the CSV empty
(header only) and add a note to `data/literature/extraction_log.md`
explaining why.

---

## BTO sources

### BTO-1. Ehsan et al 2021 (Phys. Rev. B 103, 094108)

- File: `data/Ehsan_BTO.pdf`
- Already extracted (from prior chat): two anchor points
  - 0 K: bare imaginary ω₀ = 166i cm⁻¹ → ω₀² = −424 meV²
  - 500 K: stabilized Ω = 170 cm⁻¹ → Ω(500 K) = 21.08 meV; Ω² = 444 meV²
- These two points are pre-filled in
  `data/literature/BTO/ehsan_2021/anchors.csv` (no digitization needed).
- **Ω(T) curve:** scan the paper for any T-dependent SCP Ω(T) curves for
  the cubic phase. If present, digitize into
  `data/literature/BTO/ehsan_2021/scp_omega_T.csv` with schema
  `{T_K, Omega_meV, Omega_sq_meV2}`. **Status:** Ehsan 2021 has no Ω(T)
  curve — Fig. 4 plots U(T) and TS(T); Fig. 5 plots F(T). Soft-mode info
  from this paper is exhausted by the two anchor points. (See
  `extraction_log.md` for full note.) The `scp_omega_T.csv` template
  remains as an empty optional placeholder.
- **Fig. 5 — free energies F_cubic(T) and F_tetra(T)** [optional, added
  v1 extension]
  - CSV: `data/literature/BTO/ehsan_2021/free_energy_T.csv`
  - Schema: `{T_K, F_cubic_meV, F_tetra_meV, source_figure}`
  - Purpose: second independent F(T) crossing for cubic-tetragonal T_c,
    parallel to Tadano 2022 Fig. 4.
- **Table I — structural properties (cubic + tetragonal, PBEsol)** [added
  v1 extension]
  - CSV: `data/literature/BTO/ehsan_2021/structural_properties.csv`
  - Schema: `{phase, a_angstrom, c_angstrom, c_over_a, V_cell_A3, atom,
    Wyckoff_x, Wyckoff_y, Wyckoff_z, source}`
  - Purpose: cross-check the project's `a = 4.00 Å, V_cell = 64.0 Å³` lock
    against Ehsan's PBEsol native cell; quantify the lattice gap that
    Tadano 2022's `dT_c/da ≈ 8500 K/Å` rescaling applies to.
- **Fig. 2 — cubic phonon band structure, nonthermalized (0 K) vs
  thermalized (500 K)** [added v1 extension]
  - CSV: `data/literature/BTO/ehsan_2021/cubic_dispersion.csv`
  - Schema: `{T_K, condition, branch_label, q_label, q_inv_angstrom,
    omega_cm_inv}`
  - `condition`: `nonthermalized` (0 K bare) vs `thermalized` (500 K).
  - Purpose: soft-TO dispersion near Γ along Γ-X (and optionally Γ-R)
    gives an *independent Λ² estimate* for BTO at the Ehsan calibration
    condition.
- **Table II — Born effective charges Z\* and ε_∞** [optional, medium
  priority]
  - CSV: `data/literature/BTO/ehsan_2021/born_charges_epsilon_infty.csv`
  - Schema: `{phase, atom, Z_xx, Z_yy, Z_zz, epsilon_inf_xx,
    epsilon_inf_yy, epsilon_inf_zz, source}`
  - Purpose: LO-TO splitting and polar character of the soft mode if/when
    the manuscript discusses polarization dynamics or IR oscillator
    strengths.

### BTO-2. Masuki/Tadano 2022 (arXiv:2205.08789v2)

- File: `data/Tadano_BTO.pdf`
- **Fig. 4** — free energy ℱ(T) for four phases (cubic, tetra, ortho, rhombo),
  T ∈ [350, 650] K
  - CSV: `data/literature/BTO/tadano_2022/free_energy_T.csv`
  - Schema: `{T_K, F_cubic_meV, F_tetra_meV, F_ortho_meV, F_rhombo_meV}`
  - Purpose: cross-validation of T_c via free-energy crossings
- **Fig. 15** — SCP phonon dispersions at T = 200, 600, 1000 K, at +1%
  lattice constant
  - CSV: `data/literature/BTO/tadano_2022/scp_dispersion_plus1pct.csv`
  - Schema: `{T_K, branch_label, q_label, omega_cm_inv}`
  - Purpose: lattice-corrected cross-validation only.
  - **NOTE: at +1% lattice — must be noted in all downstream uses.**
  - Extract specifically the Γ-point soft-mode value at each T. Other
    branches optional.
- **Table I** — transition temperatures
  - CSV: `data/literature/BTO/tadano_2022/transition_temperatures.csv`
  - Schema: `{transition, T_c_calc_K, T_c_exp_K}`
  - Rows: cubic-tetra, tetra-ortho, ortho-rhombo
  - Body text: `dT_c/da ≈ 8500 K/Å`; lattice-corrected cubic-tetra
    T_c = 511 K (record as note in `extraction_log.md`).
- **Fig. 7 — total κ(T) + experimental overlays** [scale-overlay only]
  - CSV: `data/literature/BTO/tadano_2022/kappa_total_T.csv`
  - Schema: `{T_K, kappa_W_per_m_K, kappa_uncertainty, method, source}`
  - `method` ∈ {`Tadano_calc`, `exp_Muta`, `exp_Popuri`, ...}.
  - **Use as scale-overlay only**, not as a κ-validation target
    (v1 explicit non-goal).

### BTO-3. Shirane et al 1970 (Phys. Rev. B 2, 3651)

- Already partially extracted from v1 into `data/raw/` as
  `ShiraneFig4a.csv` and `ShiraneFig4b.csv` (energy-loss spectra at two
  (q, e) configurations, 295 K, single domain).
- **Phase caveat:** Shirane 1970 measurements are all at 295 K, which is
  *below* T_c ≈ 403 K, so BTO is in the **tetragonal phase**. Numbers used
  for the cubic-phase soft mode in the manuscript inherit this
  approximation; flag in any downstream notebook.
- **Fig. 3 (e‖[110]) + Fig. 6 (e‖[100]) — soft-TO dispersions** [added
  v1 extension; supersedes earlier placeholder
  `polar_TO_dispersion_295K.csv`, which was empty and has been removed]
  - CSV: `data/literature/BTO/shirane_1970/cubic_dispersion.csv`
  - Schema: `{q_label, q_inv_angstrom, omega_meV, omega_uncertainty_meV,
    branch, polarization_e, wavevector_q, phase, T_K, source_figure, note}`
  - **CSV-quoting note:** `q_label` values like `[0,0.1,0]` contain commas
    and must be enclosed in double quotes (`"[0,0.1,0]"`). The current
    file is already correctly quoted via the standard `csv` module.
  - Purpose: reproducible Λ² fit from both Σ branch (Fig. 3) and Δ branch
    (Fig. 6) — gives one Λ² with uncertainty, or distinct Λ²_∥ / Λ²_⊥ if
    the two directions disagree. Replaces the inherited 5100 meV²·Å² with
    a properly sourced number.
- `data/literature/BTO/shirane_1970/Lambda_sq.csv` is pre-filled with the
  central value 5100 meV²·Å² already established in prior chat work; the
  Λ² fit from `cubic_dispersion.csv` is intended to either corroborate
  this number or replace it.
  - Schema: `{Lambda_sq_meV2_A2, Lambda_sq_uncertainty_meV2_A2, T_K, reference_q_range}`

---

## STO sources

### STO-1. Tadano-Tsuneyuki 2015 (Phys. Rev. B 92, 054301)

- File: `data/Tadano_STO.pdf`
- **This is the primary STO benchmark — extract carefully.**
- **Fig. 5 — Ω(T) for two Γ₁₅ modes** (TO1 = FE soft mode; TO2 the next
  Γ₁₅) with SCPH±PM, TDEP, and Servoin/Yamada-Shirane experimental
  overlays
  - CSV: `data/literature/STO/tadano_2015/scph_omega_T.csv`
  - Schema (extended v1): `{T_K, mode, method, Omega_meV,
    Omega_sq_meV2, source_figure}`
  - `mode` ∈ {`TO1`, `TO2`}; `method` ∈ {`SCPH_with_PM`,
    `SCPH_without_PM`, `TDEP`, `exp_Servoin`, `exp_YamadaShirane`}.
- **Fig. 4 — Cubic STO dispersion at 300 K** (SCPH + harmonic +
  experimental)
  - CSV: `data/literature/STO/tadano_2015/scph_dispersion.csv`
  - Schema (extended v1): `{T_K, q_label, q_inv_angstrom, omega_meV,
    branch, method, source_figure}`
  - `method` ∈ {`SCPH`, `harmonic`, `experimental`}.
  - Primary purpose: constrain Λ²_STO from the soft-TO dispersion near Γ.
- **Table I — soft-mode ω(300 K) vs q₁-grid + harmonic** [medium]
  - CSV: `data/literature/STO/tadano_2015/scph_convergence_300K.csv`
  - Schema: `{T_K, q_label, q_inv_angstrom, omega_meV, branch, method,
    source_figure}` (`branch` packs mode + grid; e.g.
    `FE_8x8x8`, `AFD_4x4x4`, `M3_harmonic`).
- **Fig. 7 — total κ(T) for c-STO + Muta/Popuri experimental** [optional,
  scale overlay only; not a validation target]
  - CSV: `data/literature/STO/tadano_2015/kappa_total_T.csv`
  - Schema: `{T_K, kappa_W_per_m_K, kappa_uncertainty, method, source}`
  - `method` ∈ {`Tadano_calc`, `exp_Muta`, `exp_Popuri`}.
  - **Use as scale-overlay only**, not as a quantitative κ validation
    target (v1 explicit non-goal).
- Skip Fig. 2 (LASSO IFC magnitudes), Fig. 3 (AIMD vs LASSO scatter), and
  Fig. 6 (R₂₅ AFD mode T-dep) — method internals or non-FE soft mode.

### STO-2. Verdi et al 2023 (Phys. Rev. Materials 7, L030801)

- File: `data/Verdi_STO.pdf`
- **Fig. 4(a) — TO1 Γ-point Ω(T)** (MLFF+SSCHA shifted-by-39.8 cm⁻¹ +
  INS [Ref. 39] + √-fit)
  - CSV: `data/literature/STO/verdi_2023/sscha_omega_T.csv`
  - Schema (extended v1): `{T_K, Omega_meV, Omega_sq_meV2, method,
    shift_meV, source_figure}`
  - `method` ∈ {`MLFF_SSCHA`, `DFT_SSCHA`, `INS_Ref39`, `sqrt_fit`}.
  - **Critical caveat:** the Fig. 4(a) MLFF+SSCHA curve is shifted DOWN
    by 39.8 cm⁻¹ (= 4.93 meV) to match the experimental low-T value
    (per Verdi caption). Record this as `shift_meV = 4.93` for the
    shifted curve, `shift_meV = 0` for INS. The unshifted MLFF+SSCHA
    lives in Fig. S1 of the Supporting Information.
- **Fig. 3 + Fig. 5 — phonon dispersion at 0 K (Fig. 3) and 0/300 K
  (Fig. 5)** with harmonic + DFT_SSCHA + MLFF_SSCHA + INS experimental
  - CSV: `data/literature/STO/verdi_2023/dispersion.csv`
  - Schema: `{T_K, q_label, q_inv_angstrom, omega_cm_inv, branch, method,
    supercell, source_figure}`
  - `method` ∈ {`DFT_harmonic`, `DFT_SSCHA`, `MLFF_SSCHA`, `experimental`};
    `supercell` ∈ {`2x2x2`, `3x3x3`}.
- **Fig. 4(b) — 1/ε(T)** (MLFF+SSCHA + Curie-Weiss fit) [low-medium]
  - CSV: `data/literature/STO/verdi_2023/inverse_dielectric_T.csv`
  - Schema: `{T_K, inverse_epsilon, method, source_figure}`
  - `method` ∈ {`MLFF_SSCHA`, `Curie_Weiss_fit`}.
  - LST-related to Ω²_TO1; serves as an internal consistency check on the
    Fig. 4(a) Ω(T) extraction.
- Skip Fig. 1 (schematic phase diagrams, no data) and Fig. 2 (energy
  landscape vs O displacement — defer to a future turn if the manuscript
  ever discusses the anharmonic potential shape).

### STO-3. Fauqué et al. 2022 — *Mesoscopic fluctuating domains in strontium titanate* (Phys. Rev. B 106, L140301; arXiv:2203.15495)

- File: `data/Fauque_STO.pdf`. **HIGHEST-impact STO source.**
- **Fig. 1(d) — DHO ω(T) at Γ, R, M, 4–1400 K** [HIGH]
  - CSV: `data/literature/STO/fauque_2022/dho_fits_T.csv`
  - Schema: `{T_K, q_label, q_point_rlu, mode, omega_meV, omega_uncertainty_meV, gamma_meV, gamma_uncertainty_meV, method, source_figure}`
  - `q_label` ∈ {`Gamma`, `R`, `M`, `near_Gamma_TA`, ...};
    `q_point_rlu` string (e.g. `"(2,0,0)"`, `"(1.5,0.5,0.5)"`, `"(2,-0.05,-0.05)"`);
    `mode` ∈ {`TO1`, `TA`, `R_AFD`, ...};
    `method` ∈ {`INS_DHO`, `TDEP`}.
  - **Critical verification while extracting:** does Fig. 1(d) (or any SI
    panel) plot DHO Γ(T) explicitly? If yes — this is the central STO
    Γ(T) source for the κ kernel. If only ω(T) is plotted, fill
    `gamma_meV` with `NaN` and consider email-requesting fit parameters.
- **Figs. 2(c) and 3(d)** — TA-branch DHO ω(T) at finite Q [LOW-MEDIUM,
  same CSV]: merge into `dho_fits_T.csv` with `mode = "TA"` and the
  corresponding `q_point_rlu`.
- **Fig. 1(c) — TDEP dispersion at 300 K along M-Γ-R** [MEDIUM]
  - CSV: `data/literature/STO/fauque_2022/tdep_dispersion_300K.csv`
  - Schema (same as Tadano scph_dispersion):
    `{T_K, q_label, q_inv_angstrom, omega_meV, branch, method, source_figure}`
  - Independent Λ²_STO source at 300 K.
- Skip Figs. 1(a)–(b), 2(a)–(b), 3(a)–(c), 3(e), 4 (sketches, raw scans,
  ultrasound comparison, free-energy cartoons).

### STO-4. Schmidt & Spaldin 2025 — STO ¹⁶O + ¹⁸O SSCHA (arXiv:2508.10735)

- File: `data/Schmidt_STO.pdf`.
- **Fig. 3 top-left — Ω_Eu(T) for ST¹⁶O and ST¹⁸O + experimental** [HIGH]
  - CSV: `data/literature/STO/schmidt_2025/sscha_omega_T.csv`
  - Schema: `{T_K, isotope, Omega_meV, Omega_sq_meV2, method, source_figure}`
  - `isotope` ∈ {`16O`, `18O`, `experimental`};
    `method` ∈ {`SSCHA_MLIP`, `exp_Ref50`}.
  - **Imaginary frequencies tabulated as negative `Omega_meV`** (Schmidt
    convention).
  - Third independent SSCHA Ω(T) source for STO triangulation; ST¹⁸O is
    the isotope cross-check.
- **Fig. 2 — SSCHA energy vs polarization (double well, 0 K, 16O + 18O,
  RPA-MLIP overlay)** [MEDIUM]
  - CSV: `data/literature/STO/schmidt_2025/energy_polarization.csv`
  - Schema: `{T_K, polarization_C_per_m2, E_meV_per_fu, isotope, method, source_figure}`
  - Verify polarization unit in PDF (C/m² assumed; could be µC/cm²).
- Defer Fig. 3 top-right, bottom-left, Fig. 4 (E(T), polarization(T), 2D
  Eu-mode maps) to a future turn if §3.3 anharmonic-potential discussion
  needs them. Skip Fig. 1 (eigenvector sketches).

### KTO-1. Meier, Mingo, van Roekeghem 2022 — KTaO₃ SSCHA+MLIP (arXiv:2206.08296)

- File: `data/Meier_KTO.pdf`. **New material** (KTaO₃). Strategic addition
  for the Section 4 parametric κ_soft/κ₀ equation: parallels STO Verdi
  2023 in scheme; together with Ranalli 2023 (on the reference shelf,
  not yet in `data/`) gives KTaO₃ two-source triangulation.
- **Fig. 3(a) — Ω²_TO(T): SSCHA+MLIP vs Shirane [10] experimental** [HIGH]
  - CSV: `data/literature/KTO/meier_2022/omega_T.csv`
  - Schema: `{T_K, Omega_meV, Omega_sq_meV2, method, source_figure}`
  - `method` ∈ {`SSCHA_MLIP`, `exp_Shirane_Ref10`, `Barrett_fit`}.
- **Fig. 3(b) — ε(T): calc + Fujishita [9] exp + Barrett fit** [MEDIUM]
  - CSV: `data/literature/KTO/meier_2022/inverse_dielectric_T.csv`
  - Schema (same as Verdi inverse_dielectric_T):
    `{T_K, inverse_epsilon, method, source_figure}`
  - `method` ∈ {`SSCHA_MLIP`, `exp_Fujishita_Ref9`, `Barrett_fit`}.
  - LST cross-check on KTaO₃ analogous to STO Verdi 4(b).
- **Fig. 2 — Temperature-dependent phonon spectrum of KTaO₃** [MEDIUM]
  - CSV: `data/literature/KTO/meier_2022/dispersion_T.csv`
  - Schema: `{T_K, q_label, q_inv_angstrom, omega_meV, branch, method, source_figure}`
  - Multi-T dispersion → multi-T Λ²_KTO sampling.
- Skip Fig. 1 (MLIP validation), Fig. 4 (lattice-parameter sweep at 0 K),
  Fig. 5 (epitaxial strain effects), Fig. 6 (electrostrictive tensors) —
  off-scope for the bulk unstrained kernel.

### KTO-2. Ranalli et al. 2023 — KTaO₃ MLFF+SSCHA (Adv. Funct. Mater.; `data/Ranalli_2023_KTO.pdf`)

- **Second KTaO₃ source** — strengthens the third material (was Meier-only).
  Same methodology lineage as Verdi STO (Carla Verdi co-authors both;
  MLFF+SSCHA, shift-to-match-experiment), so the extraction transfers from
  the Verdi STO pattern almost directly.
- **UNIT CONFIRMED (figure viewed):** Fig. 4(a) y-axis is **cm⁻¹** (not meV).
  Convert to meV in the CSV via ×0.1239842; log native cm⁻¹ in this plan.
  Fig. 3 (full dispersion) is also cm⁻¹.
- **Fig. 4(a) — TO1 Γ-point Ω(T), 0–300 K [HIGH; primary target]**
  - CSV: `data/literature/KTO/ranalli_2023/omega_T.csv`
  - Schema (identical to Verdi `sscha_omega_T`): `{T_K, Omega_meV,
    Omega_sq_meV2, method, shift_meV, source_figure}`.
  - `method` ∈ {`MLFF_SSCHA_bubble` (green △, main result — the curve NB 07b
    triangulates against Meier), `INS_Ref39` (black ▽ experimental anchor;
    NB: *this PDF's caption cites Ref [37]* — verify ref number),
    `sqrt_fit` (blue line, ω²_TO ∝ C(T−Tc); shows the <30 K plateau
    breakdown — secondary)}.
  - **Shift flag (parallels Verdi 4.93 / Schmidt negative-Ω):** the computed
    frequencies are plotted **shifted DOWN by 34.71 cm⁻¹ = 4.30 meV** to
    match the exp low-T soft mode (green & black coincide at ~20 cm⁻¹ at
    low T). Capture the curve **as plotted (shifted)**; record
    `shift_meV = 4.30` for `MLFF_SSCHA_bubble` + `sqrt_fit`, `0` for
    `INS_Ref39`. Do NOT pre-correct in the CSV (NB 05b plots raw + shift
    annotation; subtract back for any Cochran fit later).
  - `Omega_sq_meV2` = sign(Ω)·Ω² in meV² (for the LST cross-check).
  - `omega_uncertainty` not in schema (figure shows no error bars).
  - **phase = cubic throughout** (KTaO₃ is cubic at all T) — recorded here in
    provenance, not as a column, to keep the schema identical to Verdi/Meier
    (mechanical loader + `_SCHEMAS` cross-check). Sentinels to expect:
    low-T plateau ≈ 20 cm⁻¹ (≈ 2.48 meV); 300 K ≈ 72 cm⁻¹ (≈ 8.9 meV) for
    MLFF_SSCHA; exp ≈ 83 cm⁻¹ (≈ 10.3 meV) at 300 K.
- **Fig. 4(b) — 1/ε(T) + Curie-Weiss [MEDIUM; extract if clean]**
  - CSV: `data/literature/KTO/ranalli_2023/inverse_dielectric_T.csv`
  - Schema (identical to Verdi/Meier `inverse_dielectric_T`): `{T_K,
    inverse_epsilon, method, source_figure}`.
  - `method` ∈ {`Calculated` (computed from the shifted Ω, per caption),
    `Curie_Weiss` (blue reference line)}.
  - Purpose: KTaO₃ LST cross-check (1/ε vs Ω²) parallel to the STO one
    (R²=0.9979); shows the Barrett low-T saturation. If 4(b) is cramped and
    the points don't separate cleanly from the line, **skip** rather than
    produce a noisy extraction (NB 07b LST is a nice-to-have, not load-bearing).
- **HOLD on Λ²:** Fig. 3 is a full calculated dispersion (Γ-X-M-Γ-R, cm⁻¹)
  with a soft-TO branch near Γ — a *potential* small-q Λ² source, but it's a
  computed full-BZ band, not discrete small-q INS. Do NOT create a
  `dispersion.csv` yet; decide after NB 05b shows what's needed. KTaO₃'s Λ²
  may be bracketed rather than fit (fine for an interior-point anchor).

### STO-5. Akimov et al. 2000 — SrTiO₃ thin film E-field Raman [SKIP]

- File: `data/akimov_STO.PDF`.
- **Not used.** Thin-film STO with applied electric field; the cubic-to-
  tetragonal transition is at 120 K in the film vs 105 K bulk, so the
  zero-field curve is not a clean bulk reference. The bulk overlays in
  Fig. 2(a) (stars from Ref [8], crosses from Ref [15]) are likely
  Servoin / Yamada-Shirane, both already digitized inside the Tadano 2015
  Fig. 5 overlay (as `exp_Servoin` and `exp_YamadaShirane`). No marginal
  bulk-κ value. No folder or CSVs created.

### STO-6. Yamada-Shirane 1969 (J. Phys. Soc. Jpn. 26, 396)

- **PDF OBTAINED: `data/Shirane_1969_STO.PDF`.** The original
  paper supersedes the Shirane 1970 Fig. 1 fallback.
- **Table II — q=0 soft-mode Ω(T), 4.5–293 K [PRE-FILLED, no
  digitization].** The 19 (T, ℏω, ±) values are tabulated in the paper
  text and were transcribed directly into
  `data/literature/STO/yamada_shirane_1969/omega_T.csv` with
  `condition_direction = q_0_TO_soft`, `q_inv_angstrom = 0.00`. This is
  the **gold-standard experimental STO FE soft-mode Ω(T)** and is the
  primary INS validation curve for the STO benchmark notebook.
  - Schema: `{material, condition_direction, q_inv_angstrom, T_K,
    omega_meV, omega_uncertainty_meV}` (unchanged from v1).
- **Fig. 6 — intrinsic width Γ(T) of the soft mode at q=0 [HIGH; DIGITIZE].**
  - CSV: `data/literature/STO/yamada_shirane_1969/gamma_T.csv`
  - Schema: `{T_K, gamma_meV, gamma_uncertainty_meV}`
  - This is the **experimental STO Γ(T)** — the quantity that lets
    `eta.eta_bracket` be replaced by a fitted `eta.eta_fit`. Γ values are
    plotted, not tabulated; prior-chat reference linewidth 0.22 → 0.85 meV
    from 20 → 250 K. Fig. 5 annotates the per-temperature ℏΓ best-fit
    values used to build Fig. 6, so cross-check Fig. 5 labels against the
    Fig. 6 curve while digitizing.
- **Fig. 2 — ω(q) dispersion along [100] at several T [OPTIONAL MEDIUM].**
  Multi-T dispersion; the existing omega_T.csv already holds the 4.5 K
  [100]/[110] dispersion (Figs. 2/3). Add other-T [100] curves with
  `condition_direction = q_100_e_010` and the matching `T_K` if a
  multi-T Λ²_STO(T) sampling is wanted later.
- **Skip:** Fig. 4 (Ω(T) at q=0 — redundant with Table II; Fleury-Worlock
  Raman + LST overlay only), Figs. 7/8 (relative intensities), Table I
  (structure factors).

---

## BTO sources (v1 library acquisitions)

### BTO-4. Harada, Axe, Shirane 1971 — *Soft Modes in Cubic BaTiO₃* (Phys. Rev. B 4, 155)

- **PDF OBTAINED: `data/Shirane_1971_BTO.PDF`.** This is the
  **cubic-phase (T > T_c ≈ 403 K)** counterpart to the tetragonal-phase
  Shirane 1970 (295 K). It removes the tetragonal-phase caveat currently
  attached to the BTO cubic soft mode. **Verify the measurement
  temperature/phase of each datum from the clean PDF** (the paper targets
  the cubic phase; record `phase` and `T_K` per row).
- **Fig. 5 — TA + soft-TO dispersion along [100] [HIGH; DIGITIZE].** The
  optic branch is given as quasiharmonic (Ω_TO)² = (Ω_TO0)² + A q² — the
  `A` coefficient IS Λ²_cubic.
  - CSV: `data/literature/BTO/harada_axe_shirane_1971/cubic_dispersion.csv`
  - Schema (mirrors Shirane 1970): `{q_label, q_inv_angstrom, omega_meV,
    omega_uncertainty_meV, branch, polarization_e, wavevector_q, phase,
    T_K, source_figure, note}`
  - **CSV-quoting note:** `q_label` tuples with commas (e.g. `"[0,0.1,0]"`)
    MUST be double-quoted (see the Fauqué/Shirane lesson).
- **Fig. 8 — low-lying phonon dispersion along principal symmetry
  directions [MEDIUM].** Fold extra branches/directions into the same
  `cubic_dispersion.csv` (`source_figure = Fig_8`).
- **Λ²_cubic [optional; derive or read].**
  - CSV: `data/literature/BTO/harada_axe_shirane_1971/Lambda_sq.csv`
  - Schema (mirrors Shirane 1970): `{Lambda_sq_meV2_A2,
    Lambda_sq_uncertainty_meV2_A2, T_K, reference_q_range}`
  - Fill either by reading `A` off the Fig. 5 caption/fit or by refitting
    the digitized `cubic_dispersion.csv`. Compare to the tetragonal 5100
    meV²·Å².
- **Table I — mode frequencies, damping, coupling vs q [HIGH; DIGITIZE].**
  - CSV: `data/literature/BTO/harada_axe_shirane_1971/mode_parameters.csv`
  - Schema (long format, one row per mode): `{q_label, q_inv_angstrom,
    T_K, phase, mode, omega_meV, omega_uncertainty_meV, gamma_meV,
    gamma_uncertainty_meV, structure_factor, coupling_param, source}`
  - `mode` ∈ {`TO`, `TA`}. Gives the **cubic-phase soft-TO damping Γ** —
    the BTO analogue of the YS-1969 Fig. 6 STO Γ(T).
- **Fig. 7 — soft-TO anisotropy near q=0 [OPTIONAL MEDIUM].** (Ω)² =
  (Ω₀)² + λ(q) q² as a function of the angular orientation of q.
  - CSV: `data/literature/BTO/harada_axe_shirane_1971/anisotropy.csv`
  - Schema: `{angle_deg, q_inv_angstrom, omega_meV, omega_uncertainty_meV,
    T_K, phase, source_figure, note}`
  - Purpose: directional Λ²_∥ / Λ²_⊥ anisotropy.
- **Skip:** Figs. 1–4, 6, 9 (scan geometry, energy profiles, TA/TO
  interference methodology, resolution demonstration) — no soft-mode data
  not already captured by Table I / Figs. 5,7,8.

---

## Reviewed-and-skipped

### STO-7. Maity et al. 2025 — *Soft phonon and the central peak at the cubic-to-tetragonal transition in SrTiO₃* (PRL)

- **PDF: `data/Maity_2025_STO.pdf`. Status: REVIEWED AND SKIPPED for the
  soft-mode pipeline.** This paper studies the **antiferrodistortive (AFD)
  R₂₅ octahedral-rotation soft mode** at Q_R = (3/2, 3/2, 5/2) and the
  associated central peak at the ~105 K cubic→tetragonal transition. That
  is **not the ferroelectric Γ soft mode** the κ pipeline is built on — it
  is the same non-FE mode class explicitly skipped in Tadano 2015 (Fig. 6,
  "R₂₅ AFD mode T-dep — non-FE soft mode"). No soft-mode CSVs created.
  - Fig. 2(c) (AFD ω(T) + Γ(T) at Q_R), Figs. 3–4 (central-peak intensity,
    correlation lengths): out of scope (AFD, not FE).
- **OPTIONAL — Fig. 1(e) κ(T) scale-overlay only.** STO κ(T) for two
  defect-engineered samples (EFG, VER) + pristine STO [Ref 21]. κ overlays
  are a Turn-4 **non-goal** (scale-overlay, not a validation target), and
  these are defect samples rather than pristine. **No file created pending
  open.** If wanted, suggested schema parallel to the existing
  `kappa_total_T.csv` files: `{T_K, kappa_W_per_m_K, kappa_uncertainty,
  sample, method, source}` with `sample` ∈ {`EFG`, `VER`, `pristine_Ref21`}.

---

## Provenance table

After all extractions, the table in `data/literature/provenance_table.md`
should reflect the actual status of every quantity. Update the `Status`
column as files are filled.
