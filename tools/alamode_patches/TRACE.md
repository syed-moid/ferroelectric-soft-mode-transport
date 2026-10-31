# Diagram (c): sign, prefactor, units and vertex convention in ALAMODE c5b015c

This note traces the second-order quartic ("sunset", diagram (c)) self-energy of `anphon` from the source to the
printed linewidth, and checks one quartet by hand against a build that prints it. Line numbers refer to the
unmodified c5b015c tree.

## Files

| file | purpose |
|---|---|
| `apply_variants.py` | writes the trace build and the three validation builds into clean c5b015c copies |
| `trace_one_quartet.py` | recomputes one printed quartet by hand and compares every intermediate quantity |
| `diff_excluded_sets.py` | compares the excluded set dumped at double precision with the count from `.mesh.evec` |
| `*.diff` | the source changes of each build against c5b015c |

Builds (each its own copy of c5b015c and its own build directory; same compiler, flags and libraries as the
production build):

| tree | change |
|---|---|
| `alamode_4ph_trace` | one added print of a single fixed quartet; the sum is not changed |
| `alamode_4ph_v1omit` | the enlarged (frequency) exclusion of the audit build, plus an accumulator of the omitted terms |
| `alamode_4ph_v2strict` | strict-index exclusion, plus accumulators of the omitted terms and of the enlarged-but-not-strict terms |
| `alamode_4ph_v3dump` | the enlarged exclusion, plus a dump of every excluded quartet and every kept near-miss at full precision |

## What the code computes

`anphon/selfenergy.cpp`, `Selfenergy::selfenergy_c` (from line 296):

- **Frequency prescription.** Line 326: `omega_shift = omega + im * epsilon`, i.e. ω̃ = ω + iε with ε > 0, the
  retarded prescription.
- **Frequency sums.** Lines 358–369: `omega_sum[j] = L(Δ_j)` with L(Δ) = 1/(ω̃ − Δ) − 1/(ω̃ + Δ) and
  Δ₀ = ω₁ + ω₂ + ω₃, Δ₁ = ω₁ + ω₂ − ω₃, Δ₂ = −ω₁ + ω₂ + ω₃, Δ₃ = ω₁ − ω₂ + ω₃.
- **Thermal coefficients.** Lines 382–386: C₀ = n₁n₂ + n₂n₃ + n₃n₁ + n₁ + n₂ + n₃ + 1, C₁ = n₃n₁ + n₂n₃ + n₃ − n₁n₂,
  C₂ = n₁n₂ + n₃n₁ + n₁ − n₂n₃, C₃ = n₂n₃ + n₁n₂ + n₂ − n₃n₁. Here n = 1/(exp(ω/k_BT) − 1), with ω in Ry
  (`thermodynamics.cpp`, `fB`, line 65). `fB` returns 0 for ω < 10⁻⁸ Ry.
- **Prefactor.** Line 394: `factor = -1 / (nk^2 * 2^5 * 3)` = −1/(96 N_k²).
- **Sums.** (k₁, k₂) run over the full mesh, k₃ = q − k₁ − k₂, and (s₁, s₂, s₃) over all ordered branch triples.

The code quantity is therefore

  Σ_code(ω) = −(1/(96 N_k²)) Σ_{k₁k₂} Σ_{s₁s₂s₃} |V⁽⁴⁾_code|² Σ_j C_j L(Δ_j).

**Vertex.** `anharmonic_core.cpp`, `AnharmonicCore::V4` (lines 492–545):

  V⁽⁴⁾_code = Σ e e e e Φ⁽⁴⁾(q) / √(M M M M) / √(ω₀ω₁ω₂ω₃).

- Φ⁽⁴⁾(q) is the lattice Fourier sum of the quartic force constants per primitive cell, with no 1/N factor
  (`calc_phi4_reciprocal`).
- The vertex returns 0 if any of the four frequencies is below 10⁻⁸ Ry (line 511). Unstable modes therefore
  never enter as partners.
- Energies are in Ry with ħ = 1.

## Output

`anphon/mode_analysis.cpp`, `print_selfenergy`:

- **Printed width.** Line 473 writes `in_kayser(2.0 * self_c.imag())` under the header
  "2*Gamma (FWHM)" (line 462).
- **Printed shift.** For the bubble, line 528 writes `in_kayser(-self_a.real())`, and the printed shifted
  frequency is ω − Re Σ_code.

So the code convention is: half-width Γ = +Im Σ_code, and shift = −Re Σ_code. Therefore

  **Σ_code = −Σ^R**,

where Σ^R is the retarded self-energy of the main text: Eq. (8), D^R ∝ 1/((ω + iδ)² − Ω² − 2ΩΣ^R), pole shift
+Re Σ^R, and Γ = −Im Σ^R, Eq. (17). The bubble (`selfenergy_a`, factor +1/(16 N_k), line 222, with the two
denominators in the opposite order) follows the same convention.

**Units.**

- `in_kayser(x) = x · Ry_to_kayser`, with Ry_to_kayser = Hz_to_kayser/time_ry ≈ 109737.3. This is an angular
  frequency (in Ry, ħ = 1) expressed in cm⁻¹ as Γ̂ = Γ<sub>angular</sub>/(2πc).
- The input ε is given in cm⁻¹ and converted with the inverse factor (`integration.cpp`, line 65). It is the
  half-width of the Lorentzian that replaces each δ function: Im[1/(x + iε)] = −ε/(x² + ε²).

**Lifetime.** The transport code uses τ = 1/(2Γ) (`conductivity.cpp`, lines 430 and 450:
`0.5 / damp_tmp`). Γ is a half-width (HWHM), and 2Γ is the population relaxation rate.

**Queue.** The queue's `validate_mode4ph` divides the printed FWHM by 2 and stores Γ₄ as a HWHM in cm⁻¹.

## Sign of each quartet

With x real and ε > 0, Im[1/(x + iε)] = −ε/(x² + ε²). Hence

  Im L(Δ) = −ε/((ω − Δ)² + ε²) + ε/((ω + Δ)² + ε²) = −4εωΔ / {[(ω − Δ)² + ε²][(ω + Δ)² + ε²]},

which for ω > 0 has the sign of −Δ.

- **C₀.** Δ₀ > 0 and C₀ > 0.
- **C₁ to C₃.** Since 1 + n_i = n_i e^{x_i} with x_i = ω_i/k_BT, C₁ = n₁n₂n₃(e^{x₁+x₂} − e^{x₃}), which has the
  sign of Δ₁. C₂ and C₃ follow by relabelling.
- **Every term.** Hence C_j Im L(Δ_j) ≤ 0 for every j and every quartet, for T > 0 and stable frequencies.
- **Result.** With K = |V⁽⁴⁾_code|²/(96 N_k²) ≥ 0, Im Σ_code = −K Σ_j C_j Im L(Δ_j) ≥ 0, and the printed
  Γ₄ = +Im Σ_code ≥ 0.

In the retarded convention Σ^R = −Σ_code = +K Σ_j C_j L(Δ_j) and Γ₄ = −Im Σ^R ≥ 0. Both statements describe the
same numbers. Writing the code's negative prefactor and the definition Γ₄ = −Im Σ in the same equation is the
inconsistency to avoid.

## Vertex normalization against the manuscript (Appendix A)

Appendix A is that of the submitted manuscript (v7.3). The revised manuscript has no appendix; its SI Sec. S7.1 writes
diagram (c) with the code vertex and the 1/96 prefactor.

- **Diagram (c).** Appendix A, Eq. (A6), defines V_n = Φ_λ Π_j √(ħ/(2Ω_j)). The code vertex omits the factor
  (ħ/2)^{n/2}. For n = 4, V⁽⁴⁾_A = V⁽⁴⁾_code/4 (ħ = 1), so |V⁽⁴⁾_code|² = 16 |V⁽⁴⁾_A|² and

    1/(96 N_k²) |V⁽⁴⁾_code|² = 1/(6 N_k²) |V⁽⁴⁾_A|².

  1/6 = 1/3! is the symmetry factor of the sunset diagram, which has three equivalent internal lines summed over
  ordered branch triples.
- **Bubble.** For n = 3, |V⁽³⁾_code|² = 8 |V⁽³⁾_A|², so the bubble prefactor 1/(16 N_k) = 1/(2 N_k) with the
  Appendix A vertex.
- **Momentum sums.** The 1/N_k (bubble) and 1/N_k² (sunset) factors normalize the free wave-vector sums. The
  vertices are intensive (per primitive cell). Appendix A does not state how its vertex scales with N; an
  equation using the 1/96 prefactor must say that V⁽⁴⁾ is the code vertex, without (ħ/2)².

## One quartet by hand

Trace build, 4³ mesh, ε = 5 cm⁻¹, 300 K, external mode Γ-TO1 (q = 0, branch 4), quartet ik₁ = 5, ik₂ = 22
(ik₃ = 57), branches 5, 10, 14 (1-based). `trace_one_quartet.py` on the run's `rta.log`:

| quantity | by hand vs printed (relative difference) |
|---|---|
| n₁, n₂, n₃ | 0 |
| S₀–S₃ = L(Δ_j) | ≤ 2 × 10⁻¹⁵ |
| closed form of Im L | ≤ 5 × 10⁻¹⁶ |
| product forms of C₁–C₃ | ≤ 3 × 10⁻¹⁶ |
| \|V⁽⁴⁾\|² = Re² + Im² | 0 |
| prefactor −1/(96 · 64²) | 0 |
| ε | 5 cm⁻¹ exactly |

- **Signs.** C_j Im L(Δ_j) = −0.50, −3.55, −0.73, −3.63 (all ≤ 0).
- **Contribution.** The quartet's contribution to Σ_code has Im = +8.19 × 10⁻¹⁸ Ry (+8.99 × 10⁻¹³ cm⁻¹).
- **Unchanged sum.** The run prints 2Γ₄ = 38.304 cm⁻¹ (Γ₄ = 19.152 cm⁻¹). This equals the unmodified production
  binary at the same settings, so the added print does not change the sum.

## Conclusion

The code supports the convention "code quantity with the negative prefactor, Γ₄ = +Im Σ_code". This is
equivalent to the retarded convention Σ^R = −Σ_code with prefactor +1/(96 N_k²) and Γ₄ = −Im Σ^R.

- **Printed width.** Γ₄ is a half-width in cm⁻¹, and τ = 1/(2Γ).
- **Prefactor.** 1/96 refers to the code vertex. With the Appendix A vertex it is 1/6.
- **Changes.** No sign or factor in the code needs changing. Only the printed equation has to carry a consistent
  label.
