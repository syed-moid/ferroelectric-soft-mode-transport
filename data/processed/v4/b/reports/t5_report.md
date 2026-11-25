# T5 report: frequency-matched control (21-Nov-2026)

the analysis plan; criteria and band edges fixed before the labels in docs/v4/ANALYSIS_STEPS.md, step T5. Descriptive; no
pass criterion. Tables only; the author writes the interpretation.

## 1. Runs and checks

Node, patched anphon a0a5d30c…, OMP 12, TRISYM = 1, 16³, production inputs (rta.in identical to the T4 runs);
`t3/run_t5.sh`, log `queue/incoming/node_b/t5/run_t5.log`; outputs pulled to `queue/incoming/node_b/t5/`,
sha256 of every file checked. Label tables `t5_control/t5_labels_{A,B}_m16.txt` (label = 4 × family + band),
invariance check 0 violations (`analysis/t5_label_invariance.txt`).

| run | label table sha256 | time | max RSS | Σ bins = total, per carrier group (max rel.) | total = `.result` (max rel.) | outputs vs reference |
|---|---|---|---|---|---|---|

Carriers: the T4 carriers (primary c0.5, 16³, frozen 50 % rule): A 203, B 211 irreducible modes; κ_P weights equal to T4 (checked, rel. 10⁻⁸).

## 2. Definitions

- S_c: κ_P-weighted mean over the carriers of the share of the carrier's group-averaged linewidth from processes
  with at least one partner in class c (a process with both partners in c counted once). Shares of the total width
  (decay + absorption); of the decay part (decay processes / decay width); of the absorption part (absorption
  processes / absorption width).
- N_c: number of modes in c on the full 16³ mesh (star-weighted; 61 440 modes).
- s_c = S_c / N_c, given in **% per 1000 modes** (S_c in % × 1000 / N_c).
- Per band: rot against non-rot (neither ∪ pol; at least one partner in the union) and against neither alone;
  ratio = s_rot / s_nonrot (s_rot / s_neither).
- Bands (cm⁻¹): 0 Ω < 105; 1 105 ≤ Ω < 145; 2 145 ≤ Ω < 160; 3 Ω ≥ 160.

## 3. Share of the carriers' linewidth by partner class (family × band)

### Model A

| family | band (cm⁻¹) | N_c | S_c total (%) | of which decay (%) | of which absorption (%) | S_c of decay part (%) | S_c of absorption part (%) | s_c total | s_c decay | s_c absorption |
|---|---|---|---|---|---|---|---|---|---|---|
| neither | Ω < 105 | 2147 | 15.46 | 8.52 | 6.94 | 51.96 | 8.75 | 7.199 | 24.199 | 4.076 |
| neither | 105 ≤ Ω < 145 | 10204 | 38.60 | 14.51 | 24.09 | 42.66 | 29.99 | 3.783 | 4.181 | 2.940 |
| neither | 145 ≤ Ω < 160 | 441 | 1.82 | 0.46 | 1.36 | 1.32 | 1.64 | 4.137 | 3.004 | 3.727 |
| neither | Ω ≥ 160 | 46189 | 78.09 | 4.45 | 73.64 | 11.54 | 95.93 | 1.691 | 0.250 | 2.077 |
| pol | Ω < 105 | 0 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | — | — | — |
| pol | 105 ≤ Ω < 145 | 55 | 0.67 | 0.09 | 0.58 | 0.26 | 0.68 | 12.198 | 4.640 | 12.284 |
| pol | 145 ≤ Ω < 160 | 264 | 1.66 | 0.26 | 1.40 | 0.68 | 1.57 | 6.294 | 2.584 | 5.949 |
| pol | Ω ≥ 160 | 64 | 0.34 | 0.06 | 0.29 | 0.15 | 0.32 | 5.386 | 2.305 | 5.060 |
| rot | Ω < 105 | 660 | 17.51 | 11.48 | 6.04 | 48.45 | 7.33 | 26.536 | 73.412 | 11.100 |
| rot | 105 ≤ Ω < 145 | 804 | 3.04 | 1.21 | 1.83 | 3.42 | 2.26 | 3.785 | 4.259 | 2.817 |
| rot | 145 ≤ Ω < 160 | 612 | 1.65 | 0.50 | 1.15 | 1.33 | 1.41 | 2.690 | 2.167 | 2.300 |
| rot | Ω ≥ 160 | 0 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | — | — | — |

s_c in % per 1000 modes. Classes overlap (one process can have partners in two classes), so the S_c do not
sum to 100 %.

### Model B

| family | band (cm⁻¹) | N_c | S_c total (%) | of which decay (%) | of which absorption (%) | S_c of decay part (%) | S_c of absorption part (%) | s_c total | s_c decay | s_c absorption |
|---|---|---|---|---|---|---|---|---|---|---|
| neither | Ω < 105 | 2003 | 12.80 | 6.88 | 5.92 | 44.44 | 7.60 | 6.392 | 22.188 | 3.795 |
| neither | 105 ≤ Ω < 145 | 10276 | 38.48 | 14.98 | 23.51 | 42.95 | 29.34 | 3.745 | 4.179 | 2.855 |
| neither | 145 ≤ Ω < 160 | 276 | 1.13 | 0.25 | 0.88 | 0.68 | 1.06 | 4.084 | 2.468 | 3.835 |
| neither | Ω ≥ 160 | 46438 | 77.89 | 5.45 | 72.45 | 13.33 | 95.42 | 1.677 | 0.287 | 2.055 |
| pol | Ω < 105 | 0 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | — | — | — |
| pol | 105 ≤ Ω < 145 | 67 | 0.78 | 0.10 | 0.68 | 0.27 | 0.80 | 11.687 | 4.073 | 11.928 |
| pol | 145 ≤ Ω < 160 | 252 | 1.84 | 0.28 | 1.56 | 0.75 | 1.74 | 7.310 | 2.956 | 6.915 |
| pol | Ω ≥ 160 | 64 | 0.37 | 0.06 | 0.31 | 0.14 | 0.35 | 5.711 | 2.223 | 5.542 |
| rot | Ω < 105 | 660 | 19.08 | 11.91 | 7.16 | 51.41 | 8.76 | 28.902 | 77.895 | 13.267 |
| rot | 105 ≤ Ω < 145 | 816 | 3.57 | 1.49 | 2.08 | 3.83 | 2.61 | 4.371 | 4.688 | 3.196 |
| rot | 145 ≤ Ω < 160 | 588 | 1.77 | 0.58 | 1.20 | 1.49 | 1.51 | 3.014 | 2.540 | 2.567 |
| rot | Ω ≥ 160 | 0 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | — | — | — |

s_c in % per 1000 modes. Classes overlap (one process can have partners in two classes), so the S_c do not
sum to 100 %.

## 4. Per band: rot partners against non-rot partners, per mode

### Decay + absorption

| band (cm⁻¹) | model | N rot | s rot | N non-rot | s non-rot | s rot / s non-rot | N neither | s neither | s rot / s neither |
|---|---|---|---|---|---|---|---|---|---|
| Ω < 105 | A | 660 | 26.536 | 2147 | 7.199 | 3.69 | 2147 | 7.199 | 3.69 |
| Ω < 105 | B | 660 | 28.902 | 2003 | 6.392 | 4.52 | 2003 | 6.392 | 4.52 |
| 105 ≤ Ω < 145 | A | 804 | 3.785 | 10259 | 3.822 | 0.99 | 10204 | 3.783 | 1.00 |
| 105 ≤ Ω < 145 | B | 816 | 4.371 | 10343 | 3.790 | 1.15 | 10276 | 3.745 | 1.17 |
| 145 ≤ Ω < 160 | A | 612 | 2.690 | 705 | 4.943 | 0.54 | 441 | 4.137 | 0.65 |
| 145 ≤ Ω < 160 | B | 588 | 3.014 | 528 | 5.622 | 0.54 | 276 | 4.084 | 0.74 |
| Ω ≥ 160 | A | 0 | — | 46253 | 1.692 | — | 46189 | 1.691 | — |
| Ω ≥ 160 | B | 0 | — | 46502 | 1.678 | — | 46438 | 1.677 | — |

### Decay (share of the decay part)

| band (cm⁻¹) | model | N rot | s rot | N non-rot | s non-rot | s rot / s non-rot | N neither | s neither | s rot / s neither |
|---|---|---|---|---|---|---|---|---|---|
| Ω < 105 | A | 660 | 73.412 | 2147 | 24.199 | 3.03 | 2147 | 24.199 | 3.03 |
| Ω < 105 | B | 660 | 77.895 | 2003 | 22.188 | 3.51 | 2003 | 22.188 | 3.51 |
| 105 ≤ Ω < 145 | A | 804 | 4.259 | 10259 | 4.165 | 1.02 | 10204 | 4.181 | 1.02 |
| 105 ≤ Ω < 145 | B | 816 | 4.688 | 10343 | 4.161 | 1.13 | 10276 | 4.179 | 1.12 |
| 145 ≤ Ω < 160 | A | 612 | 2.167 | 705 | 2.843 | 0.76 | 441 | 3.004 | 0.72 |
| 145 ≤ Ω < 160 | B | 588 | 2.540 | 528 | 2.696 | 0.94 | 276 | 2.468 | 1.03 |
| Ω ≥ 160 | A | 0 | — | 46253 | 0.252 | — | 46189 | 0.250 | — |
| Ω ≥ 160 | B | 0 | — | 46502 | 0.289 | — | 46438 | 0.287 | — |

### Absorption (share of the absorption part)

| band (cm⁻¹) | model | N rot | s rot | N non-rot | s non-rot | s rot / s non-rot | N neither | s neither | s rot / s neither |
|---|---|---|---|---|---|---|---|---|---|
| Ω < 105 | A | 660 | 11.100 | 2147 | 4.076 | 2.72 | 2147 | 4.076 | 2.72 |
| Ω < 105 | B | 660 | 13.267 | 2003 | 3.795 | 3.50 | 2003 | 3.795 | 3.50 |
| 105 ≤ Ω < 145 | A | 804 | 2.817 | 10259 | 2.990 | 0.94 | 10204 | 2.940 | 0.96 |
| 105 ≤ Ω < 145 | B | 816 | 3.196 | 10343 | 2.914 | 1.10 | 10276 | 2.855 | 1.12 |
| 145 ≤ Ω < 160 | A | 612 | 2.300 | 705 | 4.559 | 0.50 | 441 | 3.727 | 0.62 |
| 145 ≤ Ω < 160 | B | 588 | 2.567 | 528 | 5.305 | 0.48 | 276 | 3.835 | 0.67 |
| Ω ≥ 160 | A | 0 | — | 46253 | 2.077 | — | 46189 | 2.077 | — |
| Ω ≥ 160 | B | 0 | — | 46502 | 2.054 | — | 46438 | 2.055 | — |

s in % per 1000 modes; — where the class is empty. Data: `analysis/t5/t5_{A,B}.json`.
