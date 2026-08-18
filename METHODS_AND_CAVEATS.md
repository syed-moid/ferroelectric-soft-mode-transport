# Methods, checks, and caveats

**Scope.** Turn 4B built the literature-benchmark data stack and the reduced
soft-mode transport kernel for the three target materials — SrTiO₃ (STO),
KTaO₃ (KTO), BaTiO₃ (BTO) — from digitized published data. This document is the
single record of *what was checked, what passed, what fired and how it was
resolved, what was retired, and what is assumed or deferred*. It is written for a
referee or a future maintainer who needs the audit trail without re-reading every
notebook.

**Closing state (v1, literature-only analysis).** `origin/main = b5699d2`. Test suite
**261 passed / 3 skipped** (`fekappa.eta` 9/9, `fekappa.kappa_kernel` 17/17); `ruff` clean.
All notebook gates PASS. No publisher PDFs tracked; all figures gitignored (regenerate from
notebooks); the one curated notebook-stage interface (`section4_eta_kappa_triples.csv`) is
tracked.

**Revision state (v2, computed SrTiO₃ core, 2026-08-18).** Test suite **297 passed /
2 skipped**; `ruff` clean (`src`, `tests`). New gates all PASS: tidy-table row counts and
locked benchmark numbers (κ(300 K) = 9.5283 W/mK; 0/525 vs 11/525 unstable; τ = 16.19/17.54
ps at 54.9 cm⁻¹); mode-sum closure vs the solver κ ≤ 0.3 % at all 21 temperatures; five-way
partition exact; Λ² fits R² ≥ 0.999 (0–400 K); kernel reduces to the v1 literature path when
fed literature inputs; HWHM-once regression extended to the computed dataset. The 300 K
production run regenerates bit-for-bit after a toolchain rebuild. Standing caveat C8 (below)
records the fixed-cell SCPH Ω₀ systematic and its handling.

---

## 1. Flags that FIRED and were RESOLVED

| # | Flag | Where | Resolution | Status |
|---|------|-------|------------|--------|
| F1 | **BTO Λ² movement > 10% gate** | NB 06 | Cubic-phase refit gave Λ²[100]=818 vs locked paper-text A(100)=972 (−15.8%, tripped the 10% gate). Resolved **by sign-off**: keep **972 central**; 818 corroborates *within* the ±20% envelope [778, 1166]. | Resolved (sign-off) |
| F2 | **KTaO₃ Ω-source / shift tie-break** | NB 07b | Ranalli MLFF+SSCHA plotted −34.71 cm⁻¹ (−4.30 meV) shifted → 9.05 meV @300 K vs Vogt 10.15 / Meier 10.86. Resolved **by Vogt case-a**: Vogt-exp 10.15 **and** Ranalli-*expt* 10.29 agree (Δ=0.14 meV) — two independent experiments — against the over-shifted MLFF; shift over-applied above ~150 K. | Resolved (data) |
| F3 | **γ-convention (FWHM vs HWHM)** | Vogt wave | DHO line-shape + magnitudes → Vogt plotted γ is **FWHM**. Locked: Γ_HWHM = γ/2, applied **once at the consumption layer** (`eta.py`); CSV stores γ as-plotted; YS γ used directly as HWHM. Regression-tested against double-application. | Resolved (locked) |
| F4 | **HAS Table I °K/°C misprint** | Wave 2 | Table-I "°K" header is a misprint (figures say °C); rows initially mis-tagged ortho/rhombo at 150/250 → re-verified against figures as **all cubic at 423/523/723 K**. | Resolved (vs figures) |
| F5 | **CSV unquoted-comma column splits** | Wave 1 / 2 | Fauqué `dho_fits_T.csv` `q_point_rlu` (e.g. `(2,0,0)`) and HAS `cubic_dispersion.csv` `[0,0,0]` split across columns → re-quoted. | Resolved (fixed) |
| F6 | **`gamma_meV` sanity bound too tight** | inputs tests | 20 → 60 → finalized **50 meV** (from HAS Table I Γ_TO=34±6). `angle_deg` (0,180)→(−90,90). | Resolved (bound set) |
| F7 | **STO dense Vogt cross-assignment** | STO dense | 3 sub-105 K Ω points mis-tagged `cubic_F1u` → re-tagged `transition_merged` (cubic_F1u ≥110 K, matching γ); uncertainty boundary 53.7→60 K; column-order swap fixed; Flag-2 rewritten as a single shared-γ sub-105 K branch. | Resolved (corrected) |
| F8 | **STO γ-scale third read** | STO dense | @300 K: dense WPD 19.4 / eye ~21–22±2 / round 24.5 cm⁻¹. Decided **keep dense (marker-center)**, *not* corrected toward the round hand value → η_STO=0.109 (see caveat C1). | Resolved (kept dense) |

---

## 2. Flags / framings RETIRED (superseded)

| Retired item | Superseded by |
|--------------|---------------|
| **KTaO₃-Γ-absent** — Amendment-1 case-a/b/c *bracket* framing (no KTaO₃ experimental Ω(T)) | Vogt 1995 experimental hyper-Raman Ω(T)+γ(T) acquired → Vogt case-a (F2) |
| **γ-convention-unconfirmed** — "factor-of-2 unresolved" placeholder | F3 FWHM verdict locked |
| **η as a *span*** (η_STO…η_BTO continuum) | η as a **cluster**: two underdamped QPE (≈0.11–0.13) + one overdamped FE (≈2.3), empty interior |
| **BTO Λ² from digitized 859 / scalar Λ²** | HAS 1971 cubic-phase **directional** Λ² (soft [100]=972, stiff [110]/[111]=4750) |

---

## 3. Checks that PASSED (routine gates)

| Notebook / module | Gate | Key results |
|-------------------|------|-------------|
| NB 05b pre-flight | 22/22 | full CSV inventory + schema cross-check |
| NB 06 BTO | PASS | directional Λ² anisotropy ~4.9× (972 vs 4750); overdamping η=2.24–2.39 (mean 2.29) |
| NB 07 STO | PASS | 6-paper/7-curve Ω_TO1(T) triangulation; cross-method q=0 anchors agree **15.4% ≤ 30%**; LST R²=**0.9978**; η_STO=**0.109** |
| NB 07b KTaO₃ | PASS | Meier+Ranalli+Vogt triangulation; η_KTO=**0.126** (<0.5); LST (Ranalli internal) R²=**0.9934**; Λ² bracketed |
| NB 08 cluster | 9/9 | cluster \|Δη\|=**0.016**; interior (0.5,2.0) empty; κ_soft ≪ total (~0.1 vs ~10) |
| NB 09 appendix | 7/7 | negative result reproduced (δ=**+2806%** @400 K) + fully disclaimed |
| `fekappa.eta` | 9/9 | η_STO=0.1096, η_KTO=0.1257, cluster, γ/2-applied-once regression, YS-direct-HWHM |
| `fekappa.kappa_kernel` | 17/17 | K_soft units, n_B/C limits, Debye q_max, κ∝1/η, overdamped guard, Ω²<0 skip, dict bracket |

---

## 4. Standing CAVEATS (carry forward)

- **C1 — η_STO = 0.109 ± 0.02** (reading uncertainty). Dense marker-center value;
  an eye-read upper bound is ~0.125. Both keep STO firmly underdamped and within the
  QPE cluster. *(See F8.)*
- **C2 — Λ²_STO = 889 meV²Å² is an EFFECTIVE branch curvature, not a small-q fit.**
  Verdi resolves only Γ/X/M/R (no points in (0, 0.15 Å⁻¹)); 889 = Verdi 3×3×3
  MLFF_SSCHA Γ→X-TO1 curvature, which spans **outside** the small-q regime
  Ω²=Ω₀²+Λ²q² (the same reason Ehsan zone-boundary points were rejected for BTO).
  STO ±30%; the 861↔889 (2×2×2 vs 3×3×3) grid spread sits inside it.
- **C3 — KTaO₃ Λ² is an STO-proxy ±50%.** No KTaO₃ finite-q soft dispersion exists;
  all KTaO₃ sources are q=0. Wider envelope than STO by design.
- **C4 — Section-4 fit interior is EXTRAPOLATED, not anchored.** The η-axis interior
  (0.5 < η < 2.0) is empty: the dataset is two underdamped QPE + one overdamped FE.
  Any Section-4 relation through that interior is interpolation between disjoint
  regimes, not data-anchored.
- **C5 — BTO κ_soft is ORDER-OF-MAGNITUDE ONLY (overdamped).** η_BTO≈2.27 (HAS,
  verified) ⇒ the K_P quasiparticle picture (well-defined v_q, τ_q) fails; the kernel
  emits `OverdampedModeWarning`. BTO Ω₀ is taken from the HAS cubic TO at the
  nearest-Γ point (q≈0.313 Å⁻¹) — disclosed, negligible against the order-of-magnitude
  framing. The real BTO value needs the Wigner-coherent term (deferred, see §5).
- **C6 — Vogt γ/2 is applied exactly once**, at the `eta.py` consumption layer. CSVs
  store γ as-plotted (FWHM); never pre-halve, never pass an already-HWHM value with
  `convention="FWHM"`. Guarded by a regression test.
- **C7 — Two-variant γ schema is intentional** (STO 4-col with `phase_mode`, KTaO₃
  3-col) — not to be reconciled.
- **C8 — Ω₀ fixed-cell SCPH systematic (computed STO dataset).** The SCPH-coupled
  ALAMODE dataset (`data/processed/sto_alamode/`) carries the known fixed-volume
  SCPH overestimate of the SrTiO₃ soft mode: computed Ω₀(300 K) = 167.1 cm⁻¹ vs
  Vogt's hyper-Raman 88.4 cm⁻¹ (ratio 1.89, growing toward low T; ratios tabulated,
  never tuned away). The dataset reproduces the fixed-cell ("without PM") SCPH curve
  of Tadano & Tsuneyuki, Phys. Rev. B **92**, 054301 (2015) — ≤1.5% at the digitized
  0/400 K points (`tadano_fig5_match.csv`); the gap to experiment is discussed in
  that paper itself and in Verdi *et al.*, Phys. Rev. Materials **7**, L030801
  (2023) and Schmidt & Spaldin, arXiv:2508.10735. Closing it requires the
  cubic-bubble real-part shift, thermal expansion, and/or an SSCHA-grade treatment
  — out of scope here. Consequences handled as follows: (i) the zone-wide η map
  stays on computed inputs, with the caveat that the soft-branch Ω overestimate
  biases η *downward* there (narrow-width conclusions optimistic for the soft
  sector); a sensitivity table rescaling soft-branch Ω by the Vogt/computed ratio
  is emitted (`eta_validity_sensitivity.csv`); (ii) the soft kernel is evaluated
  both fully computed and literature-anchored (Vogt/Verdi Ω₀(T) with computed Γ₀,
  Λ²); (iii) κ(T) validation is **unaffected** — κ is acoustic-dominated and the
  9.53 vs 8.7 W/mK (300 K) agreement stands independently of the soft-mode Ω₀
  systematic.

---

## 5. DEFERRED (future work)

- **Wigner coherent κ_C** — the wave-like term, required for the overdamped BTO soft
  mode (where K_P is invalid) and for interband-spacing-comparable-to-linewidth cases.
- **Full SCPH self-consistency loop** — `fekappa.scph.scph_iterate` is a stub; the
  aborted scalar-SCPH attempt is preserved under `archive/scph_calibration/`
  (and is what NB 09 reproduces). A real loop awaits new first-principles force constants.
- **The Section-4 fit itself** — this work produced only the **anchor triples**
  (`data/processed/section4_eta_kappa_triples.csv`); fitting the κ(η) relation is future work.
- **κ_bte_rta / κ_total integrators** (`fekappa.kappa`) — still stubs.

---

## 6. Provenance

Per-result provenance — which source, notebook, and convention underlies each value — is
given in the tables above, in the notebooks (`notebooks/`), and in the digitized-source
records under `data/literature/` (`provenance_table.md`, `extraction_log.md`,
`digitization_plan.md`).
