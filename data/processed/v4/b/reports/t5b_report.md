# T5b report: control inside the lowest band ()

the analysis plan; band edges fixed before the labels in docs/v4/ANALYSIS_STEPS.md, step T5b (terciles of model A's rot
frequencies below 105 cm⁻¹, rounded to 5 cm⁻¹; bands above 105 merged). Same quantities as T5
(the T5 report tables (Markdown)) plus the frequency of each class. Descriptive; no pass criterion. Tables only;
the author writes the interpretation.

## 1. Runs and checks

Node, patched anphon a0a5d30c…, OMP 12, TRISYM = 1, 16³, 300 K, production inputs (rta.in identical to the
T4 runs); `t3/run_t5b.sh`, log `queue/incoming/node_b/t5b/run_t5b.log`; outputs pulled to
`queue/incoming/node_b/t5b/`, sha256 of every file checked. Label tables `t5b_control/t5b_labels_{A,B}_m16.txt`
(label = 4 × family + band), invariance check 0 violations (`analysis/t5b_label_invariance.txt`).

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
- Bands (cm⁻¹): 0 Ω < 80; 1 80 ≤ Ω < 95; 2 95 ≤ Ω < 105; 3 Ω ≥ 105.

## 3. Share of the carriers' linewidth by partner class (family × band)

### Model A

| family | band (cm⁻¹) | N_c | S_c total (%) | of which decay (%) | of which absorption (%) | S_c of decay part (%) | S_c of absorption part (%) | s_c total | s_c decay | s_c absorption |
|---|---|---|---|---|---|---|---|---|---|---|
| neither | Ω < 80 | 453 | 4.51 | 2.59 | 1.92 | 32.76 | 2.42 | 9.955 | 72.322 | 5.343 |
| neither | 80 ≤ Ω < 95 | 658 | 4.77 | 2.63 | 2.14 | 10.13 | 2.69 | 7.255 | 15.402 | 4.088 |
| neither | 95 ≤ Ω < 105 | 1036 | 6.51 | 3.62 | 2.88 | 11.47 | 3.65 | 6.280 | 11.073 | 3.522 |
| neither | Ω ≥ 105 | 56834 | 93.30 | 17.68 | 75.62 | 51.06 | 97.92 | 1.642 | 0.898 | 1.723 |
| pol | Ω < 80 | 0 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | — | — | — |
| pol | 80 ≤ Ω < 95 | 0 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | — | — | — |
| pol | 95 ≤ Ω < 105 | 0 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | — | — | — |
| pol | Ω ≥ 105 | 383 | 2.68 | 0.41 | 2.27 | 1.08 | 2.57 | 6.990 | 2.833 | 6.710 |
| rot | Ω < 80 | 240 | 13.57 | 9.21 | 4.36 | 40.68 | 5.21 | 56.547 | 169.505 | 21.709 |
| rot | 80 ≤ Ω < 95 | 180 | 2.35 | 1.32 | 1.03 | 4.32 | 1.26 | 13.036 | 23.983 | 6.973 |
| rot | 95 ≤ Ω < 105 | 240 | 2.30 | 1.07 | 1.23 | 3.95 | 1.44 | 9.578 | 16.457 | 6.009 |
| rot | Ω ≥ 105 | 1416 | 4.59 | 1.65 | 2.94 | 4.59 | 3.63 | 3.241 | 3.239 | 2.566 |

s_c in % per 1000 modes. Classes overlap (one process can have partners in two classes), so the S_c do not
sum to 100 %.

### Model B

| family | band (cm⁻¹) | N_c | S_c total (%) | of which decay (%) | of which absorption (%) | S_c of decay part (%) | S_c of absorption part (%) | s_c total | s_c decay | s_c absorption |
|---|---|---|---|---|---|---|---|---|---|---|
| neither | Ω < 80 | 393 | 3.58 | 2.14 | 1.44 | 27.13 | 1.84 | 9.117 | 69.021 | 4.671 |
| neither | 80 ≤ Ω < 95 | 598 | 3.93 | 2.08 | 1.85 | 9.58 | 2.39 | 6.568 | 16.024 | 4.003 |
| neither | 95 ≤ Ω < 105 | 1012 | 5.56 | 2.91 | 2.64 | 9.66 | 3.38 | 5.492 | 9.548 | 3.341 |
| neither | Ω ≥ 105 | 56990 | 93.34 | 18.63 | 74.71 | 52.10 | 97.70 | 1.638 | 0.914 | 1.714 |
| pol | Ω < 80 | 0 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | — | — | — |
| pol | 80 ≤ Ω < 95 | 0 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | — | — | — |
| pol | 95 ≤ Ω < 105 | 0 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | — | — | — |
| pol | Ω ≥ 105 | 383 | 2.99 | 0.44 | 2.55 | 1.16 | 2.90 | 7.808 | 3.029 | 7.562 |
| rot | Ω < 80 | 264 | 15.22 | 9.78 | 5.44 | 44.09 | 6.53 | 57.653 | 167.004 | 24.743 |
| rot | 80 ≤ Ω < 95 | 180 | 2.36 | 1.24 | 1.11 | 4.29 | 1.38 | 13.087 | 23.827 | 7.656 |
| rot | 95 ≤ Ω < 105 | 216 | 2.25 | 0.99 | 1.27 | 3.40 | 1.50 | 10.434 | 15.753 | 6.955 |
| rot | Ω ≥ 105 | 1404 | 5.25 | 2.00 | 3.24 | 5.16 | 4.08 | 3.737 | 3.676 | 2.908 |

s_c in % per 1000 modes. Classes overlap (one process can have partners in two classes), so the S_c do not
sum to 100 %.

## 4. Per band: rot partners against non-rot partners, per mode

### Decay + absorption

| band (cm⁻¹) | model | N rot | s rot | N non-rot | s non-rot | s rot / s non-rot | N neither | s neither | s rot / s neither |
|---|---|---|---|---|---|---|---|---|---|
| Ω < 80 | A | 240 | 56.547 | 453 | 9.955 | 5.68 | 453 | 9.955 | 5.68 |
| Ω < 80 | B | 264 | 57.653 | 393 | 9.117 | 6.32 | 393 | 9.117 | 6.32 |
| 80 ≤ Ω < 95 | A | 180 | 13.036 | 658 | 7.255 | 1.80 | 658 | 7.255 | 1.80 |
| 80 ≤ Ω < 95 | B | 180 | 13.087 | 598 | 6.568 | 1.99 | 598 | 6.568 | 1.99 |
| 95 ≤ Ω < 105 | A | 240 | 9.578 | 1036 | 6.280 | 1.53 | 1036 | 6.280 | 1.53 |
| 95 ≤ Ω < 105 | B | 216 | 10.434 | 1012 | 5.492 | 1.90 | 1012 | 5.492 | 1.90 |
| Ω ≥ 105 | A | 1416 | 3.241 | 57217 | 1.641 | 1.98 | 56834 | 1.642 | 1.97 |
| Ω ≥ 105 | B | 1404 | 3.737 | 57373 | 1.636 | 2.28 | 56990 | 1.638 | 2.28 |

### Decay (share of the decay part)

| band (cm⁻¹) | model | N rot | s rot | N non-rot | s non-rot | s rot / s non-rot | N neither | s neither | s rot / s neither |
|---|---|---|---|---|---|---|---|---|---|
| Ω < 80 | A | 240 | 169.505 | 453 | 72.322 | 2.34 | 453 | 72.322 | 2.34 |
| Ω < 80 | B | 264 | 167.004 | 393 | 69.021 | 2.42 | 393 | 69.021 | 2.42 |
| 80 ≤ Ω < 95 | A | 180 | 23.983 | 658 | 15.402 | 1.56 | 658 | 15.402 | 1.56 |
| 80 ≤ Ω < 95 | B | 180 | 23.827 | 598 | 16.024 | 1.49 | 598 | 16.024 | 1.49 |
| 95 ≤ Ω < 105 | A | 240 | 16.457 | 1036 | 11.073 | 1.49 | 1036 | 11.073 | 1.49 |
| 95 ≤ Ω < 105 | B | 216 | 15.753 | 1012 | 9.548 | 1.65 | 1012 | 9.548 | 1.65 |
| Ω ≥ 105 | A | 1416 | 3.239 | 57217 | 0.896 | 3.62 | 56834 | 0.898 | 3.61 |
| Ω ≥ 105 | B | 1404 | 3.676 | 57373 | 0.910 | 4.04 | 56990 | 0.914 | 4.02 |

### Absorption (share of the absorption part)

| band (cm⁻¹) | model | N rot | s rot | N non-rot | s non-rot | s rot / s non-rot | N neither | s neither | s rot / s neither |
|---|---|---|---|---|---|---|---|---|---|
| Ω < 80 | A | 240 | 21.709 | 453 | 5.343 | 4.06 | 453 | 5.343 | 4.06 |
| Ω < 80 | B | 264 | 24.743 | 393 | 4.671 | 5.30 | 393 | 4.671 | 5.30 |
| 80 ≤ Ω < 95 | A | 180 | 6.973 | 658 | 4.088 | 1.71 | 658 | 4.088 | 1.71 |
| 80 ≤ Ω < 95 | B | 180 | 7.656 | 598 | 4.003 | 1.91 | 598 | 4.003 | 1.91 |
| 95 ≤ Ω < 105 | A | 240 | 6.009 | 1036 | 3.522 | 1.71 | 1036 | 3.522 | 1.71 |
| 95 ≤ Ω < 105 | B | 216 | 6.955 | 1012 | 3.341 | 2.08 | 1012 | 3.341 | 2.08 |
| Ω ≥ 105 | A | 1416 | 2.566 | 57217 | 1.720 | 1.49 | 56834 | 1.723 | 1.49 |
| Ω ≥ 105 | B | 1404 | 2.908 | 57373 | 1.712 | 1.70 | 56990 | 1.714 | 1.70 |

s in % per 1000 modes; — where the class is empty. Data: `analysis/t5b/t5b_{A,B}.json`.

## 5. Frequency of each class (full 16³ mesh, unweighted over the class's modes; cm⁻¹)

| model | family | band (cm⁻¹) | N | Ω mean | Ω min | Ω max |
|---|---|---|---|---|---|---|
| A | neither | Ω < 80 | 453 | 63.13 | 0.00 | 79.41 |
| A | neither | 80 ≤ Ω < 95 | 658 | 88.45 | 80.24 | 94.46 |
| A | neither | 95 ≤ Ω < 105 | 1036 | 100.70 | 95.21 | 104.87 |
| A | neither | Ω ≥ 105 | 56834 | 337.10 | 105.40 | 830.66 |
| A | pol | Ω < 80 | 0 | — | — | — |
| A | pol | 80 ≤ Ω < 95 | 0 | — | — | — |
| A | pol | 95 ≤ Ω < 105 | 0 | — | — | — |
| A | pol | Ω ≥ 105 | 383 | 153.27 | 133.69 | 165.63 |
| A | rot | Ω < 80 | 240 | 64.34 | 25.82 | 79.86 |
| A | rot | 80 ≤ Ω < 95 | 180 | 88.11 | 80.32 | 93.67 |
| A | rot | 95 ≤ Ω < 105 | 240 | 100.33 | 95.02 | 104.78 |
| A | rot | Ω ≥ 105 | 1416 | 137.19 | 106.88 | 158.41 |
| B | neither | Ω < 80 | 393 | 61.55 | 0.00 | 75.31 |
| B | neither | 80 ≤ Ω < 95 | 598 | 87.58 | 80.11 | 93.52 |
| B | neither | 95 ≤ Ω < 105 | 1012 | 100.40 | 95.09 | 104.63 |
| B | neither | Ω ≥ 105 | 56990 | 337.27 | 105.11 | 832.69 |
| B | pol | Ω < 80 | 0 | — | — | — |
| B | pol | 80 ≤ Ω < 95 | 0 | — | — | — |
| B | pol | 95 ≤ Ω < 105 | 0 | — | — | — |
| B | pol | Ω ≥ 105 | 383 | 152.70 | 133.42 | 166.81 |
| B | rot | Ω < 80 | 264 | 64.15 | 22.85 | 78.44 |
| B | rot | 80 ≤ Ω < 95 | 180 | 88.50 | 83.16 | 93.26 |
| B | rot | 95 ≤ Ω < 105 | 216 | 99.85 | 96.08 | 103.99 |
| B | rot | Ω ≥ 105 | 1404 | 136.32 | 105.90 | 157.56 |

Per-band unions (as in §4):

| model | partners | band (cm⁻¹) | N | Ω mean | Ω min | Ω max |
|---|---|---|---|---|---|---|
| A | rot | Ω < 80 | 240 | 64.34 | 25.82 | 79.86 |
| A | non-rot (neither ∪ pol) | Ω < 80 | 453 | 63.13 | 0.00 | 79.41 |
| A | neither | Ω < 80 | 453 | 63.13 | 0.00 | 79.41 |
| A | rot | 80 ≤ Ω < 95 | 180 | 88.11 | 80.32 | 93.67 |
| A | non-rot (neither ∪ pol) | 80 ≤ Ω < 95 | 658 | 88.45 | 80.24 | 94.46 |
| A | neither | 80 ≤ Ω < 95 | 658 | 88.45 | 80.24 | 94.46 |
| A | rot | 95 ≤ Ω < 105 | 240 | 100.33 | 95.02 | 104.78 |
| A | non-rot (neither ∪ pol) | 95 ≤ Ω < 105 | 1036 | 100.70 | 95.21 | 104.87 |
| A | neither | 95 ≤ Ω < 105 | 1036 | 100.70 | 95.21 | 104.87 |
| A | rot | Ω ≥ 105 | 1416 | 137.19 | 106.88 | 158.41 |
| A | non-rot (neither ∪ pol) | Ω ≥ 105 | 57217 | 335.87 | 105.40 | 830.66 |
| A | neither | Ω ≥ 105 | 56834 | 337.10 | 105.40 | 830.66 |
| B | rot | Ω < 80 | 264 | 64.15 | 22.85 | 78.44 |
| B | non-rot (neither ∪ pol) | Ω < 80 | 393 | 61.55 | 0.00 | 75.31 |
| B | neither | Ω < 80 | 393 | 61.55 | 0.00 | 75.31 |
| B | rot | 80 ≤ Ω < 95 | 180 | 88.50 | 83.16 | 93.26 |
| B | non-rot (neither ∪ pol) | 80 ≤ Ω < 95 | 598 | 87.58 | 80.11 | 93.52 |
| B | neither | 80 ≤ Ω < 95 | 598 | 87.58 | 80.11 | 93.52 |
| B | rot | 95 ≤ Ω < 105 | 216 | 99.85 | 96.08 | 103.99 |
| B | non-rot (neither ∪ pol) | 95 ≤ Ω < 105 | 1012 | 100.40 | 95.09 | 104.63 |
| B | neither | 95 ≤ Ω < 105 | 1012 | 100.40 | 95.09 | 104.63 |
| B | rot | Ω ≥ 105 | 1404 | 136.32 | 105.90 | 157.56 |
| B | non-rot (neither ∪ pol) | Ω ≥ 105 | 57373 | 336.04 | 105.11 | 832.69 |
| B | neither | Ω ≥ 105 | 56990 | 337.27 | 105.11 | 832.69 |
