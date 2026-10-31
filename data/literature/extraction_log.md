# Literature extraction log — v1

Updated as sources are digitized. Each entry records what was extracted,
which figure(s), how many points, any anomalies, and the WebPlotDigitizer
calibration parameters so each extraction is reproducible.

---

## BTO-1. Ehsan et al 2021 (PRB 103, 094108)

**Status:** anchors pre-filled; Fig. 5 + Table I + Fig. 2 + Table II
extracted; `scp_omega_T.csv` confirmed n/a.

**Pre-filled anchors** (no digitization, `BTO/ehsan_2021/anchors.csv`):
- 0 K bare imaginary: 166i cm⁻¹ → ω₀² = -424 meV² (Lock 1).
- 500 K thermalized-displacement SCP: Ω = 170 cm⁻¹ = 21.08 meV;
  Ω² = 444 meV² (Lock 2).

**Fig. 5 — F_cubic(T) and F_tetra(T)** (`free_energy_T.csv`):
- 20 (T, F_cubic) points across 300–600 K; F_tetra at 11 points across
  300–455 K (F_tetra is only plotted up to the crossing; above the
  crossing, the tetragonal phase is unstable so Ehsan does not show
  F_tetra — those cells are `NaN`).
- **F_cubic = F_tetra ≈ 0 at T ≈ 455 K** → Ehsan's cubic-tetra crossing.
- This is the BTO analog of Tadano 2022 Fig. 4; the two crossings can be
  compared as an independent cross-validation of the calculated T_c.

**`scp_omega_T.csv` — n/a:** Ehsan 2021 contains no Ω<sub>soft</sub>(T) curve.
Fig. 4 = vibrational U(T) and TS(T); Fig. 5 = free energies F_cubic(T)
and F_tetra(T). Soft-mode frequency information from this paper is
exhausted by the two anchor points already captured in `anchors.csv`. The
template stays in place as an empty optional CSV.

**Table I — structural properties (`structural_properties.csv`):**
- 10 rows: 5 atoms × 2 phases (cubic at expt T=473 K; tetragonal at expt
  T=293 K). Cubic a = 4.013 Å (V = 64.63 Å³); tetragonal a = 4.002 Å,
  c = 4.028 Å (c/a = 1.007, V = 64.51 Å³).
- **Cross-check vs. project lock (a = 4.00 Å, V = 64.0 Å³):** Ehsan's
  cubic a is 0.013 Å higher (+0.32%), V_cell is 1.0% larger. Small enough
  that the lock's transport-kernel volume prefactor is not materially
  off; large enough that lattice corrections matter for Tadano-2022 +1%
  comparisons (see Tadano notes below).

**Fig. 2 — cubic phonon band structure (`cubic_dispersion.csv`):**
- 16 rows across both `condition = nonthermalized (0 K)` and
  `thermalized (500 K)`, for Highest Optical, Ferroelectric Soft Mode,
  and Acoustic TA/LA branches at Γ, X, M, R.
- **Soft-mode at Γ matches the two anchors:** 0 K nonthermalized
  ω = −166 cm⁻¹ (imaginary, matches the bare 166i lock); 500 K
  thermalized ω = +170 cm⁻¹ (matches the 21.08 meV anchor).
- Soft-mode at X: 0 K = −130 cm⁻¹, 500 K = +160 cm⁻¹.
- **v1 unit conversion (applied):** the `q_inv_angstrom` column was
  originally populated with reduced-coordinate magnitudes (0, 0.5, 0.707,
  0.866 = |Γ|, |X|, |M|, |R| in r.l.u.). Detected by NB 05b pre-flight
  gate. Converted in-place by multiplying with `2π/a ≈ 1.5657 Å⁻¹` using
  the Ehsan Table I cubic `a = 4.013 Å`; the column now genuinely holds
  Å⁻¹ values: Γ=0, X=0.7829, M=1.1070, R=1.3559.
- **v1 methodological decision:** these finite-q points are at the
  **zone boundary** (X, M, R), well outside the small-q regime where the
  long-wavelength expansion `Ω²(q) = Ω₀² + Λ²q²` holds. They are
  therefore **NOT used for the BTO Λ² refit** in NB 06. The apparent
  "negative slope" at 500 K thermalized (Ω(X)=160 < Ω(Γ)=170 cm⁻¹) is
  almost certainly SCP soft-branch anti-crossing with another mode along
  Γ→X, not a real Λ² < 0. The dispersion data retains standalone value as
  a shape illustration of the Ehsan SCP scheme but is omitted from the
  long-wavelength kernel inputs.

**Table II — Born charges and ε_∞ (`born_charges_epsilon_infty.csv`):**
- 7 rows: cubic (Ba, Ti, O1) + tetragonal (Ba, Ti, O1, O2).
- Cubic Z*(Ti) = 7.26 (large, drives polar character of the soft mode).
- Cubic ε_∞ = 6.00 (isotropic); tetragonal ε_∞_xx = 5.74, ε_∞_zz = 5.38
  (anisotropic, as expected for the polarized phase).

---

## BTO-2. Masuki/Tadano 2022 (arXiv:2205.08789v2)

**Status:** Fig. 4 + Fig. 15 + Table I all extracted.

**Fig. 4 — free energy ℱ(T), four phases (`free_energy_T.csv`):**
- 9 rows across 350–580 K covering F_cubic, F_tetra, F_ortho, F_rhombo.
- Cubic-tetra crossing visible around 410 K in the digitization.

**Fig. 15 — SCP dispersion at +1% lattice (`scp_dispersion_plus1pct.csv`):**
- 28 rows: Highest Optical + Soft Optical + Acoustic TA/LA at 200, 600,
  1000 K, mostly at Γ, X, M, R, and a Γ<sub>mid</sub> sampling point.
- **Soft optical at Γ<sub>mid</sub>:** 40 cm⁻¹ (200 K), 140 cm⁻¹ (600 K), 175 cm⁻¹
  (1000 K). In meV: 5.0 / 17.4 / 21.7. **Critical caveat: at +1% lattice.**
  Without rescaling via `dT_c/da ≈ 8500 K/Å` these are NOT directly
  comparable to Ehsan's native-lattice 500 K anchor.

**Table I — transition temperatures (`transition_temperatures.csv`):**
- 3 rows: cubic-tetra (606 K calc / 390 K exp), tetra-ortho (509 / 270),
  ortho-rhombo (411 / 180).
- Body text: `dT_c/da ≈ 8500 K/Å`; lattice-corrected cubic-tetra
  T_c = 511 K (worth recording when these values get used downstream).

**Fig. 7 — total κ(T) + experimental (`kappa_total_T.csv`):**
- 28 rows: Tadano calculated curve (200-800 K, monotonically decreasing
  ~10.3 → 5.4 W/m·K) plus Popuri and Muta experimental data (with
  uncertainties on the experimental points).
- **Use as scale-overlay only**, not as κ-validation target (v1
  explicit non-goal).

---

## BTO-3. Shirane et al 1970 (PRB 2, 3651)

**Status:** Λ² central value pre-filled; Fig. 3 + Fig. 6 dispersion data
extracted; old empty placeholder `polar_TO_dispersion_295K.csv` removed.

**Phase caveat (all rows):** Shirane 1970 measurements are at 295 K,
which is *below* T_c ≈ 403 K, so BTO is in the **tetragonal phase**. The
Λ² fit inherits this approximation; flag in any downstream use that the
Shirane-derived Λ² is a tetragonal-phase value applied to the cubic-phase
soft-mode treatment.

**Pre-filled `Lambda_sq.csv`:** 5100 meV²·Å² central value from an earlier
reference fit (no uncertainty stored).

**Figs. 3 & 6 — soft-TO dispersion (`cubic_dispersion.csv`):**
- 21 rows (12 from Fig. 6 e‖[100]; 9 from Fig. 3 e‖[110]).
- TO points from Fig. 6 (Δ branch, q‖[010] or q‖[001]) span q = 0.39–0.63 Å⁻¹
  with ω = 5.5–7.8 meV (with ~1–2 meV error bars where provided).
- TO points from Fig. 3 (Σ branch, q‖[1̄10]) span q = 0.11–0.22 Å⁻¹ with
  ω = 9.0–17.0 meV — much steeper dispersion than Fig. 6 (consistent with
  a longitudinal Coulomb contribution along that direction).
- TA points are recorded for context (lower curves).
- **CSV-quoting note (resolved):** `q_label` values like `[0,0.1,0]`
  contain commas and were quoted with the standard csv module so the
  schema's 11 columns parse correctly.
- **Action for v1:** fit Λ² separately from the Fig. 3 and Fig. 6 TO
  branches; if the two values agree within uncertainty, fold into a
  single Λ² with that uncertainty; if they disagree (Λ²_∥ vs Λ²_⊥),
  report both in `Lambda_sq.csv` with `reference_q_range` annotations.

---

## STO-1. Tadano-Tsuneyuki 2015 (PRB 92, 054301)

**Status:** Fig. 5 + Fig. 4 + Table I extracted; Fig. 7 (κ overlay)
optional, still pending.

**Fig. 5 — Ω(T) for two Γ₁₅ modes (`scph_omega_T.csv`):**
- 30 rows: TO1 (FE soft mode) + TO2 across four methods (SCPH_with_PM,
  SCPH_without_PM, TDEP, exp_Servoin, exp_YamadaShirane).
- TO1 SCPH_with_PM: 13.64 meV at 0 K → 22.44 meV at 1000 K (monotone).
- TO1 exp_YamadaShirane spans 15-250 K (2.23 → 9.67 meV); exp_Servoin
  spans 300-920 K (10.91 → 18.60 meV). These two experimental sets cover
  STO complementary T ranges — together they bracket the SCPH curve.
- Schema extended v1 to `{T_K, mode, method, Omega_meV,
  Omega_sq_meV2, source_figure}`.

**Fig. 4 — Cubic STO dispersion at 300 K (`scph_dispersion.csv`):**
- 14 rows at 300 K: Highest_Optical (~99 meV at Γ and M, ~92 meV at X),
  Soft_Optical at Γ<sub>mid</sub> = 12.4 meV (SCPH; matches exp), Structural soft
  R-point = 6.2 meV SCPH (vs -11.2 meV harmonic, i.e., imaginary harmonic
  R₂₅ mode stabilized by SCPH).
- Methods covered: SCPH, harmonic, experimental.
- Schema extended v1 to add `method` and `source_figure` columns.

**Table I — Soft-mode convergence at 300 K (`scph_convergence_300K.csv`):**
- 21 rows: FE (Γ), AFD (R), M₃ (M) soft modes at q₁ grids 2×2×2 through
  12×12×12, plus harmonic reference for each.
- FE soft converges from 17.85 (2×2×2) to 16.74 meV (10×10×10, 12×12×12).
  Harmonic value is -7.19 meV (imaginary).
- AFD soft converges from 8.55 (2×2×2) to 4.34 meV (12×12×12); harmonic
  -9.42 meV. Tighter q-grid pulls AFD down significantly.
- Same schema as `scph_dispersion.csv`; `branch` column packs mode + grid
  (`FE_8x8x8`, `AFD_harmonic`, etc.).

**Fig. 7 — total κ(T) and Muta/Popuri exp (`kappa_total_T.csv`):**
- PENDING (optional). Schema `{T_K, kappa_W_per_m_K, kappa_uncertainty,
  method, source}`. **Use as scale-overlay only**, not as a κ-validation
  target (v1 explicit non-goal).

---

## STO-2. Verdi et al 2023 (PRM 7, L030801)

**Status:** Fig. 4(a) + Figs. 3/5 dispersion + Fig. 4(b) all extracted.

**Fig. 4(a) — TO1 Γ Ω(T) (`sscha_omega_T.csv`):**
- 28 rows across three method curves:
  - MLFF_SSCHA at 11 T points (0-300 K). **shift_meV = 4.93** (= 39.8 cm⁻¹
    downward shift per Verdi caption to match low-T experimental value).
  - INS_Ref39 at 11 T points (5-300 K), `shift_meV = 0`.
  - sqrt_fit at 7 T points (2-300 K), `shift_meV = 4.93` (inherits the
    MLFF shift since it's a fit through the shifted MLFF curve).
- At 300 K: MLFF (shifted) 8.93 meV vs INS 10.17 meV vs sqrt_fit (shifted)
  9.17 meV. Disagreements stay within the shift magnitude.
- Schema extended v1 to add `shift_meV` column for traceability.
  Unshifted MLFF+SSCHA values live in Verdi Fig. S1; not yet captured.

**Figs. 3 + 5 — phonon dispersion (`dispersion.csv`):**
- 22 rows.
- Fig. 3 (0 K, 2×2×2 supercell): DFT_harmonic shows soft TO1 at Γ
  imaginary (-160 cm⁻¹); DFT_SSCHA and MLFF_SSCHA stabilize it to
  +40 cm⁻¹. Confirms the bubble-term renormalization. Highest optical
  ~530 cm⁻¹ at X, ~890 cm⁻¹ at M, ~940 cm⁻¹ at R.
- Fig. 5 (0 K and 300 K, 3×3×3 supercell): MLFF_SSCHA static dispersion +
  experimental overlay. 0 K Γ soft optical = 22 cm⁻¹ MLFF_SSCHA vs 0 cm⁻¹
  exp; 300 K Γ soft optical = 115 cm⁻¹ MLFF_SSCHA vs 85 cm⁻¹ exp.

**Fig. 4(b) — 1/ε(T) (`inverse_dielectric_T.csv`):**
- 26 rows: 19 MLFF_SSCHA T points (10-300 K) + 7 Curie_Weiss_fit T points.
- 1/ε ranges from 0.0006 (10 K) to 0.0243 (300 K) for MLFF_SSCHA.
  Curie-Weiss line agrees to within ~1% across the full range.
- Useful for LST cross-check: the Cochran relation 1/ε ∝ Ω²_TO1 should
  hold; comparing slope vs the Ω(T) data in `sscha_omega_T.csv` checks
  internal consistency.

---

## STO-3. Fauqué et al. 2022 — Mesoscopic fluctuating domains in strontium titanate (Phys. Rev. B 106, L140301; arXiv:2203.15495)

**Status:** PENDING (templates created v1). **HIGHEST-impact STO source.**

- `dho_fits_T.csv` — TODO. Target ≥ 15 (T, ω) points at Γ across 4–1400 K
  plus the R and M points from Fig. 1(d). Merge in Figs. 2(c) and 3(d)
  TA-branch points with `mode = "TA"` and the appropriate `q_point_rlu`.
  **Critical verification while extracting:** does Fig. 1(d) or any SI
  panel plot DHO Γ(T) explicitly? If yes → this becomes the central STO
  Γ(T) source for the κ kernel and replaces the T¹/T² bracketing scheme
  for `eta.eta_T_STO`. If only ω is plotted, fill `gamma_meV` with `NaN`
  and consider email-requesting the DHO fit parameters from the authors.
- `tdep_dispersion_300K.csv` — TODO. Target ≥ 10 (q, ω) points along
  M-Γ-R from Fig. 1(c) for soft TO and other relevant branches.

---

## STO-4. Schmidt & Spaldin 2025 — STO ¹⁶O + ¹⁸O SSCHA (arXiv:2508.10735)

**Status:** PENDING (templates created v1).

- `sscha_omega_T.csv` — TODO. Target ≥ 8 (T, ω) points per isotope across
  the plotted T range, plus any experimental points from Ref [50] with
  `isotope = "experimental"`. **Imaginary frequencies tabulated as
  negative `Omega_meV`** per Schmidt convention.
- `energy_polarization.csv` — TODO (optional). Verify the polarization
  axis unit in the PDF (C/m² assumed; could be µC/cm²).

---

## KTO-1. Meier et al. 2022 — KTaO₃ SSCHA+MLIP (arXiv:2206.08296)

**Status:** PENDING (templates created v1). **New material** —
Section 4 parametric κ<sub>soft</sub>/κ₀ fit third pillar (alongside BTO + STO).

- `omega_T.csv` — TODO. Target ≥ 10 (T, ω²) points from Fig. 3(a),
  separated by `method` (SSCHA_MLIP calc vs Shirane exp).
- `inverse_dielectric_T.csv` — TODO. Fig. 3(b) calc + Fujishita + Barrett.
- `dispersion_T.csv` — TODO. Fig. 2 multi-T dispersion. Pick the soft TO
  branch at ≥ 3 temperatures for Λ²_KTO sampling.

**Lattice constants for KTaO₃** (locks needed for the kernel):
cubic perovskite KTaO₃ at $a \approx 3.989$ Å → $V_{\rm cell} \approx 63.5$
Å³. **Action**: confirm against Meier 2022 Table or body text before NB 06
takes the value as a lock.

---

## STO-5. Akimov et al. 2000 — STO thin-film E-field Raman

**Status:** REVIEWED AND SKIPPED. Thin-film STO under electric field;
cubic-to-tetragonal at 120 K (vs 105 K bulk). The zero-field Fig. 2(a)
curve is not a clean bulk reference. The bulk overlays in Fig. 2(a)
(stars = Ref [8], crosses = Ref [15]) are already digitized inside the
Tadano 2015 Fig. 5 overlay as `exp_Servoin` and `exp_YamadaShirane`. No
folder or CSVs created. No marginal value for the bulk κ pipeline.

---

## STO-6. Yamada-Shirane 1969 (JPSJ 26, 396) — and the Shirane 1970 Fig. 1 fallback

**Status:** PDF NOT IN `data/`. Partial fallback available: Shirane 1970
Fig. 1 redraws the Yamada-Shirane STO 4 K dispersion (the BTO panel in
Fig. 1 is schematic Hüller theory and is NOT a measurement — don't
digitize it).

- `omega_T.csv` — schema updated to 6 columns `{material,
  condition_direction, q_inv_angstrom, T_K, omega_meV,
  omega_uncertainty_meV}` to accommodate ω(q, T) dispersion data rather
  than ω(T) at fixed q. TODO: digitize the STO panel of Shirane 1970
  Fig. 1 (or the original Yamada-Shirane 1969 paper if obtained).
- `gamma_T.csv` — schema unchanged. TODO if Yamada-Shirane PDF obtained;
  earlier reference: linewidth 0.22 → 0.85 meV from 20 K to 250 K.

---

## v1 library acquisitions (Yamada-Shirane 1969, Harada-Axe-Shirane 1971, Maity 2025)

### STO-6. Yamada-Shirane 1969 (`data/Shirane_1969_STO.PDF`)

**Status:** Table II PRE-FILLED into `omega_T.csv`; Fig. 6 Γ(T) pending.

- **Table II → `STO/yamada_shirane_1969/omega_T.csv`** (pre-filled, no
  WPD): 19 q=0 soft-mode points, 4.5–293 K, transcribed verbatim from the
  paper text (`condition_direction = q_0_TO_soft`, `q_inv_angstrom = 0.00`,
  uncertainties = the tabulated ±). Internally consistent with the
  pre-existing 4.5 K dispersion rows (ω(q=0, 4.5 K) = 1.9 meV in both).
  This is now the gold-standard experimental STO Ω(T) validation curve.
- **Fig. 6 → `gamma_T.csv`** (TODO HIGH): experimental STO Γ(T) at q=0;
  plotted not tabulated — digitize. Fig. 5 annotates the per-T ℏΓ best-fit
  values feeding Fig. 6 (cross-check). Replaces `eta.eta_bracket` with a
  fitted `eta.eta_fit` once filled.
- Skipped: Fig. 4 (redundant with Table II), Figs. 7/8 (intensities),
  Table I (structure factors).

### BTO-4. Harada-Axe-Shirane 1971 (`data/Shirane_1971_BTO.PDF`)

**Status:** folder + 4 templates created; all pending digitization.
Cubic-phase (T > T_c) BTO — removes the tetragonal-phase caveat on the
Shirane 1970 (295 K) soft mode. VERIFY per-row phase/T from the clean PDF.

- **Fig. 5 + Fig. 8 → `cubic_dispersion.csv`** (TODO HIGH): soft-TO + TA
  dispersion along [100] (Fig. 5) and full low-lying dispersion (Fig. 8).
  Optic branch is quasiharmonic (Ω<sub>TO</sub>)² = (Ω<sub>TO0</sub>)² + A q²; `A` = Λ²_cubic.
- **Table I → `mode_parameters.csv`** (TODO HIGH): long format, mode ∈
  {TO, TA}; gives cubic-phase soft-TO damping Γ + ω + coupling vs q.
- **Lambda_sq.csv** (optional): fill from Fig. 5 `A` or refit
  `cubic_dispersion.csv`; compare to tetragonal 5100 meV²·Å².
- **Fig. 7 → `anisotropy.csv`** (optional): directional Λ²_∥/Λ²_⊥.
- Skipped: Figs. 1–4, 6, 9 (geometry / profiles / interference / resolution).

### STO-7. Maity et al. 2025 (`data/Maity_2025_STO.pdf`) — REVIEWED AND SKIPPED

AFD R₂₅ octahedral-rotation soft mode at Q_R = (3/2,3/2,5/2) + central
peak at the ~105 K cubic→tetragonal transition — NOT the ferroelectric Γ
soft mode (same non-FE class skipped in Tadano 2015 Fig. 6). No soft-mode
CSVs created. Only Fig. 1(e) κ(T) (defect EFG/VER + pristine [21]) is
marginally relevant, and κ overlays are a Turn-4 non-goal — optional
`kappa_defect_T.csv` deferred.

---

## v1 Wave 2 — library-source verification (HAS 1971 + YS 1969)

Verified against the rendered PDF figures (poppler), not the CSV alone.

### HAS 1971 (cubic BaTiO3) — physics review

- **Soft TO is OVERDAMPED (paper finding, verified).** Table I (image) +
  p.159 text: the DHO mode-damping Γ<sub>TO</sub> (Eq. 2, the renormalized/
  quasiharmonic damping — NOT a FWHM/HWHM or column artifact) satisfies the
  overdamping condition Γ<sub>TO</sub> ≳ 2 Ω<sub>TO</sub> at every (q, T). Empirical
  η = Γ<sub>TO</sub>/Ω<sub>TO</sub> ≈ 2.3 (neutron, finite q); cited 4–6 at q=0 (optical) and 9
  at q=0.08 a*. Example (q=0.2 a*, 450 °C): Ω<sub>TO</sub> = 15±3, Γ<sub>TO</sub> = 34±6 meV →
  η = 2.27. **Direct evidence cubic BTO is overdamped; the underdamped
  soft-mode framing is STO-specific** (→ NB 06 BTO methodological note,
  alongside the Hlinka 2008 two-mode complication).
- **CSV CORRECTION (mode_parameters.csv): phase + temperature.** Table I's
  header prints "T (°K)" with 150/250/450, which had led to ortho/rhombo
  phase tags. But Figs. 4/5/6/8 all label 150/250/450 **°C**, the text gives
  "Ω₀=(3±1) meV at 150 °C", and the Ω hardening (9.8→12.7→15.0) is the
  classic cubic soft mode above T_c=130 °C. → all rows are **cubic**, at
  T = 423/523/723 K (the "°K" header is a misprint). Fig. 5 independently
  labels its TA data "150 °C / 230 °C", confirming the °C convention (and
  that cubic_dispersion.csv's 503 K open circles = 230 °C are correct).
- **Cubic Λ² (Lambda_sq.csv updated 859 → 972±50).** Fig. 5 fits the soft TO
  to Ω²(q)=Ω₀²+A q² with Ω₀=(3±1) meV; p.160 gives **A(100)=(972±50)**,
  A₁(110)=A(111)=(4750±300) meV²Å². The soft-TO sits at q=0/0.313/0.478 Å⁻¹
  (ζ≤0.3) — small-q, usable for the fit, NOT zone-boundary-only like Ehsan.
  The Fig-5 3-point digitization gave 859 (~12% low); the text value 972 is
  authoritative. **This is the first CUBIC-phase BTO Λ²**, resolving the NB
  05b finding that all prior BTO Λ² sources were tetragonal. The tetragonal
  Shirane 5100 is NOT in conflict: it matches the steep A(110)/A(111)=4750
  direction — the gap is the soft-TO anisotropy (Fig. 7), so "Λ²_BTO" must
  carry a direction in NB 06.
- **Fig. 7 angle convention.** Abscissa = direction along arcs
  [111]–[100]–[110] ([100] at the Λ minimum). The CSV signed encoding
  ([100]=0°, [111]=−54.74°, [110]=+45°) faithfully parameterizes it;
  sanity bound (−90, 90) is a valid superset (data span −54.74° to +45°).

### YS 1969 (cubic SrTiO3) — Γ(T) verification

- **Fig. 6 Γ(T) VERIFIED (image).** All 8 gamma_T.csv points match: 0.24
  (20 K) → 0.84 meV (250 K), monotonic; symbols • E₀=35 meV, ■ E₀=13 meV.
  Canonical low-T experimental STO soft-mode linewidth → lets `eta.eta_fit`
  exist alongside `eta.eta_bracket` (NB 08). STO is strongly **underdamped**:
  η = Γ/Ω ≈ 0.087 at 250 K (Γ=0.84, Ω=9.6 from Table II) — contrast cubic
  BTO η≈2.3. Table II q=0 Ω(T) (pre-filled v1) verified consistent.

### Γ sanity-range bound finalized

`gamma_meV` upper bound set to **50 meV** (was the reactive 60): HAS Table I
max Γ<sub>TO</sub>=34±6 (→40); optical Γ/Ω=9 at q=0.08 a* with Ω₀≈3 implies ~44 meV at
small q. 50 brackets the overdamped cubic-BTO regime with margin.

---

## KTaO₃ wave — Ranalli 2023 (second KTaO₃ source)

`data/Ranalli_2023_KTO.pdf`. Loaders + NB 05b verification added; NB 06/07/eta/
kernel untouched (parallel data-verification wave). Both CSVs digitized by user.

- **`omega_T.csv` (Fig. 4(a), native cm⁻¹ → meV).** `method` ∈ {`MLFF+SSCHA`,
  `expt`, `fit_paraelectric`}. The computed `MLFF+SSCHA` curve is plotted
  **shifted DOWN by 34.71 cm⁻¹ = 4.30 meV** to match the exp low-T soft mode;
  recorded as `shift_meV = -4.30` (negative = down-shift; same convention
  family as Verdi −4.93 and Schmidt negative-Ω). `expt`/`fit_paraelectric`
  unshifted. Recover unshifted Cochran Ω via Ω − shift_meV (+4.30). phase =
  cubic throughout (KTaO₃ cubic at all T; kept in provenance, not a column).
- **`inverse_dielectric_T.csv` (Fig. 4(b)).** `method` ∈ {`MLFF+SSCHA`, `expt`,
  `Curie-Weiss_fit`}. Feeds the KTaO₃ LST 1/ε vs Ω² cross-check (NB 07b).
- **Two-source consistency (verified, NB 05b Fig. F):** the INS experimental
  anchors AGREE across sources (Ranalli expt 300 K = 10.29 vs Meier Shirane
  10.70 meV, Δ=0.41). The computed curves carry a ~1.8 meV spread (Ranalli
  shifted MLFF 9.05 vs Meier Calc 10.86 @300 K; low-T plateau Ranalli 2.48 vs
  Meier 3.33) — a genuine inter-source shift-convention/method spread, **noted
  not asserted-equal**. QPE low-T plateau confirmed (Ω(5 K)=Ω(20 K)=2.48 meV).
- **Λ² HELD:** Fig. 3 is a full calculated dispersion (Γ-X-M-Γ-R, cm⁻¹) — a
  potential near-Γ Λ² source, but Γ-point-only in the digitized data like
  Meier; KTaO₃ Λ² to be bracketed, not fit. No `dispersion.csv` created.

---

## Vogt wave — Vogt 1995 (PRB 51, 9426): experimental Ω(T)+γ(T) for BOTH QPEs

`data/Vogt_1995_KTO.pdf` (+ `Vogt_1988_KTO.pdf`). Loaders + NB 05b verification;
NB 06/07/eta/kernel untouched. Implements Amendment-2 §5. **Amendment 1's
KTaO₃ Γ-source case-analysis (a/b/c) and "bracket/interior anchor" framing are
RETIRED** — superseded by the resolutions below.

### γ-convention VERDICT (the gating item)
Vogt's Fig. 5/6 fit (Sec. III) = "a damped harmonic oscillator response **and a
δ function** at the center" → the **central peak is separated**; the plotted γ is
the intrinsic zone-center DHO damping (kernel-appropriate). The plotted symbol is
lowercase **γ** (standard DHO damping constant), distinct from the theoretical
Eq. (13) self-energy **Γ** (which is HWHM). Reading γ as the standard DHO constant
⇒ **γ = FWHM ⇒ Γ<sub>HWHM</sub> = γ/2**. Confirmed by the numbers (STO 300 K: Ω=10.54,
γ=3.04 → η=0.144 on the FWHM basis; the HWHM reading would give 0.288). **Residual
factor-of-2 systematic NOTED** (if Vogt's fit used his Eq.-13 Γ=HWHM, use γ
directly and η doubles); Refs 29-30 (PRB 34,410 (1986); 41,1184 (1990)) — not in
repo — would settle it. Both readings keep STO & KTaO₃ underdamped (η<0.5) and
clustered; only the absolute η calibration carries the factor. `eta.py` will apply
Γ<sub>HWHM</sub> = γ/2 by default (later wave); loaders return plotted γ as-is.

### Schema normalization (to the Amendment-2 §3 locked schema)
- KTO omega: method `experimental_KTaO3`→`exp_HyperRaman_Vogt1995`;
  `exp_HyperRaman_Vogt` (Table III 300 K=10.54)→`exp_HyperRaman_Vogt1988`.
- STO omega: split `experimental_SrTiO3_{A2u,Eu,F1u}` → `method=exp_HyperRaman_Vogt1995`
  + **`phase_mode`** column (`tetragonal_A2u`/`tetragonal_Eu` <105 K; `cubic_F1u` ≥110 K).
- **Two-variant gamma schema (intentional):** STO-Vogt gamma 4-col (with phase_mode),
  KTaO₃-Vogt + YS gamma 3-col. NOT reconciled.
- omega schema uses `Omega_uncertainty_meV` (not Omega_sq); `shift_meV = 0` (exp).

### Findings (NB 05b verified)
- **η CLUSTER, not span:** STO η≈0.144, KTaO₃ η≈0.148 (FWHM basis; |Δ|=0.004) —
  two underdamped QPEs cluster; BTO η≈2.3 sits alone overdamped; η-axis interior
  (~0.3–2.3) is EMPTY → Section 4 interior is extrapolated, not anchored.
- **Vogt = PRIMARY η source** (true q=0, central peak separated) for both QPEs →
  eta_fit. YS 1969 = documented cross-check.
- **YS-vs-Vogt STO γ spread:** raw 3.39× at 250 K → ~1.7× after FWHM→HWHM /2;
  documented physics (YS finite-q + central-peak-contaminated vs Vogt true-q=0,
  central peak separated), NOT a data error — do not reconcile to equality.
- **Vogt tie-breaks Meier–Ranalli:** Vogt exp Ω(300 K)=10.48 ≈ Meier 10.86 ≫
  Ranalli shifted-MLFF 9.05 → Ranalli's −4.30 meV downshift is slightly
  over-applied above ~150 K.

### Flags / open items
- **SPOT-CHECK (KTaO₃ Fig. 5 damping):** digitized high-T γ reads slightly HIGH —
  295 K γ=3.10 meV=25 cm⁻¹ vs the Fig. 5 right axis ~20-22 cm⁻¹ (152 K point also
  non-monotonic, likely scatter). Impact: η<sub>KTO</sub> ≈ 0.24 (figure) vs 0.30 (CSV);
  **clustering unaffected**. NOT auto-adjusted (experimental values; re-digitization
  of the KTaO₃ damping tail recommended as a focused follow-up).
- **STO tetragonal A₂ᵤ/E_u γ:** Vogt fit the sub-105 K split with the SAME γ
  (p.8051), but the digitized A₂ᵤ/E_u γ rows differ (e.g. 0.12 vs 0.08 at 10 K).
  Sub-105 K is discussion-only (cubic η uses `cubic_F1u` ≥110 K), so low-impact;
  marked for optional re-digitization, not changed.
- `KTO/vogt_1995/model_params.md` is METADATA (Barrett-fit description), no loader.

---

## Vogt wave — RE-VERIFICATION of corrected KTaO₃ γ tail (supersedes the two flags above)

The KTaO₃ `vogt_1995/gamma_T.csv` tail was re-digitized and independently
re-verified against Vogt 1995 Fig. 5 (rendered figure). **Both prior flags are
now resolved**; uncertainty convention settled.

- **KTaO₃ γ tail corrected & VERIFIED vs Fig. 5 (within ~1 cm⁻¹):** 295 K
  3.10 → **2.60 meV = 21.0 cm⁻¹** (figure top damping point ~21–22 cm⁻¹);
  5 K = 0.19 meV = 1.5 cm⁻¹ (figure ~1.2–1.6 cm⁻¹ via the right-axis mapping);
  60–295 K monotonic. The earlier high-T pixel-shift is gone.
- **140–170 K SMOOTHING (noted):** Fig. 5 shows *genuine scatter* there
  (~±1–1.5 cm⁻¹); the corrected CSV is monotonic (13.1/13.6/14.2 cm⁻¹). The
  smoothing is **sub-error-bar** (within the ±1 cm⁻¹ = 0.12 meV bar) but the raw
  figure is not monotonic — noted, not "fixed".
- **Uncertainty convention (settled, both QPEs):** ±1 cm⁻¹ ≈ **0.12 meV** base
  (stat) for T ≥ 60 K; **0.17 meV** for T < 60 K = quadrature stat⊕sys
  √(0.12² + 0.12²), the extra 0.12 meV (≈1 cm⁻¹) being the **defect-induced
  systematic floor** Vogt notes for "below about 60 K" (p.8051: even nominally
  pure samples show soft-mode distortion around symmetry-breaking defects). The
  60 K boundary matches Vogt's wording verbatim. Applied **consistently to all
  four Vogt files** (KTaO₃ + STO, γ + Ω).
- **Cross-file consistency:** STO Vogt γ and Ω now use the same 0.12/0.17
  convention (were 0.05/0.08) — desync fixed.
- **1988 cross-ref point exception:** KTaO₃ omega's `exp_HyperRaman_Vogt1988`
  Table III point (300 K = 10.54 meV) legitimately keeps its own ±0.37 meV
  (≈3 cm⁻¹); the 0.12/0.17 convention is 1995-data only. (NB 05b assertion
  filters to 1995 rows.)
- **Flag-2 RESOLVED:** STO sub-105 K A₂ᵤ/E_u γ now share one value at each T
  (10 K both 0.10, 50 K 0.28, 100 K 0.62), per Vogt's same-damping fit (p.8051).
- **η impact:** η<sub>KTO</sub>(300 K) = (γ/2)/Ω = (2.60/2)/10.48 = **0.124** (was 0.148);
  still clusters with η<sub>STO</sub> = 0.144 (|Δ| = 0.020 < 0.05). Convention applied at
  the consumption layer (Figure J / eta.py), **NOT baked into the CSV** (CSV
  stores the plotted FWHM γ as digitized).

---

## Vogt wave — FINAL: dense WebPlotDigitizer KTaO₃ files authoritative (supersedes prior addenda)

The KTaO₃ `vogt_1995/{omega_T,gamma_T}.csv` were re-done as **dense WPD
extractions** (~39 pts each, native non-round T) and are now authoritative;
they **supersede the round-number / 0.84-regime-split version entirely**.
Re-verified against Fig. 5 (rendered):

- **Boundary (the gating check) PASSES — no scaling artifact.** 81.59 K →
  0.868 meV (7.0 cm⁻¹), 106.61 K → 1.201 meV (9.7 cm⁻¹): a smooth 7→9.7 cm⁻¹
  rise matching Fig. 5. The dense file was **never scaled**, so the flat
  83–105 K junction of the old 0.84-split version is gone. (Corrects the prior
  "0.84 axis-recalibration" framing — there is none here.)
- **Ω re-digitized too** (corrects the earlier "Ω unchanged at 10.48"): the
  1995 curve now tops at 294.36 K → **10.15 meV** (5.19 K → 2.44 plateau);
  there is no 1995 300 K point. The **1988 Table III** point (300 K = 10.54 ±
  **0.37**) is retained as an independent 300 K cross-check — *consistent within
  the 1988 ±0.37* (not the old "0.06 match").
- **γ-file checks (vs figure):** (a) 11.5 K = 0.093 meV dips below 5.85 K
  (0.147) — kept (sub-0.17-bar scatter; figure scatters at low T). (b)
  295.99/297.89 K both 2.570 = **single point double-picked → one dropped**
  (figure shows one isolated point above the 20 gridline). (c) 55.62/56.56 and
  245.35/249.10 — distinct values, two real points, kept. (bonus) 267.29/276.79
  both 2.361 — figure shows two distinct dots at one height → real plateau, kept.
  γ file: 40 → **39 rows** after the drop.
- **Uncertainty:** sub-60 K = 0.17 meV = √(0.12²+0.12²) (stat ±1 cm⁻¹ ⊕ sys
  ±1 cm⁻¹ defect floor, Vogt p.8051 "below about 60 K"); ≥60 K = 0.12; applied
  to all dense Vogt files (55.62/56.56 K both < 60 → 0.17). 1988 point keeps
  ±0.37 (1995-only convention).
- **η (dense):** η<sub>KTO</sub> = (γ/2)/Ω = (2.57/2)/10.15 = **0.127** (HWHM basis),
  still clusters with η<sub>STO</sub> = 0.144 (|Δ| = 0.018). Convention applied at the
  consumption layer (Figure J / eta.py), not in the CSV.
- STO Vogt files unchanged (round-T): A₂ᵤ/E_u share one γ (Flag-2), sub-60 K
  0.17. NB 05b gate PASS 22/22.

---

## Vogt wave — STO dense correction pass (supersedes the round STO Vogt files)

The STO `vogt_1995/{omega_T,gamma_T}.csv` are now **dense WebPlotDigitizer
extractions** (omega 47 rows, gamma 36), authoritative; they supersede the
round-number STO Vogt version. Four corrections applied + re-verified vs Fig. 6:

1. **Ω phase re-tag:** 3 sub-105 K points (85.89/93.58/97.77 K) were tagged
   `cubic_F1u` but lie below T_a = 105 K (the A₂ᵤ/E_u split is merely
   *unresolved* above ~85 K, not cubic). Re-tagged → `transition_merged`, so
   `cubic_F1u` now starts ≥110 K in BOTH files (boundary consistent with γ);
   the cubic-η filter is phase-clean.
2. **Uncertainty boundary** moved from ~53.7 K to **60 K** (γ 53.69/59.91,
   Ω 53.78 → 0.17), matching Vogt's "below about 60 K" and the KTaO₃ files.
3. **γ scale — THIRD independent figure read (decided, not reconciled):** three
   reads of Fig. 6 damping at 300 K — dense WPD **19.4 cm⁻¹**, eye-read
   **~21–22 cm⁻¹ (±2)**, round hand-read **24.5 cm⁻¹**. The eye-read rejects
   the round but sits ~2 cm⁻¹ above the dense. **Decision: keep the dense
   (marker-center-accurate) value; NOT corrected toward the round.** Both QPE
   dense extractions read lower than their round hand-reads (round likely read
   to scatter-top, not marker-center).
4. **Non-monotonic γ points** (223.89→1.80, 248.63→2.05 meV): confirmed real
   plotted scatter in Fig. 6 (~±1 cm⁻¹, within the 0.12 meV bar) → kept.

**η<sub>STO</sub> (LOCKED) = 0.109** (HWHM basis, Γ<sub>HWHM</sub>=γ/2 over Ω), **central; carry
~±0.02 reading uncertainty** (0.109 dense central → ~0.125 eye-read upper; dense
adopted as marker-center-accurate). **Cluster is robust across the whole range:**
η<sub>STO</sub> ≈ 0.11–0.13 vs η<sub>KTO</sub> ≈ 0.127, both ≪ 0.5 (underdamped), BTO ≈ 2.3 alone.

Knock-ons: (a) **YS-vs-Vogt STO damping spread shrank to ~2.4× raw** (was ~3.3×
on the round γ; ~1.2× residual after FWHM→HWHM /2) — qualitative point unchanged
(Vogt > YS; q-resolution + central peak), magnitude smaller (note in NB 07 §3.3).
(b) **Ω F1u top** 299 K → 10.96 meV (η denominator). (c) Fixed a column-order
swap in the dense omega (`phase_mode`/`source_figure`); made NB/test T-hardcodes
dense-robust; Flag-2 re-expressed as "single shared-γ sub-105 K branch" (the
dense gamma has no separate E_u row — the shared-γ fit is one curve).
