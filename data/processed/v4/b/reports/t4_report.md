## 2. Results (opened here first)

All shares are κ_P-weighted means over the carriers (each carrier mode weighted by its κ_P contribution), of the
carrier's group-averaged three-phonon linewidth (decay + absorption). Soft = pol, rot or mixed. 300 K, 16³, q₁ 12³,
TRISYM = 1 (production widths).

### 2a. Carriers (frozen 50 % rule)

| | model A | model B |
|---|---|---|
| κ_P mode sum (W/mK) | 8.8907 | 8.6850 |
| carrier modes (irreducible) / degenerate groups | 203 / 198 | 211 / 205 |
| share of κ_P carried by the carriers | 50.03 % | 50.05 % |
| carriers by family, primary_c0.5 (neither / pol / rot / mixed) | 154 / 0 / 49 / 0 | 162 / 0 / 49 / 0 |
| carriers by family, secondary_c0.5 (neither / pol / rot / mixed) | 153 / 1 / 48 / 1 | 161 / 1 / 48 / 1 |
| carriers by family, primary_c0.4 (neither / pol / rot / mixed) | 131 / 0 / 72 / 0 | 139 / 0 / 72 / 0 |
| carriers by family, primary_c0.6 (neither / pol / rot / mixed) | 165 / 0 / 38 / 0 | 174 / 0 / 37 / 0 |

### 2b. Linewidth of the carriers by number of soft partners

| definition | model | neither soft (%) | one soft (%) | both soft (%) | absorption part (%) |
|---|---|---|---|---|---|
| primary_c0.5 | A | 76.13 | 21.36 | 2.51 | 77.69 |
| primary_c0.5 | B | 73.71 | 23.37 | 2.91 | 76.99 |
| secondary_c0.5 | A | 68.85 | 28.11 | 3.04 | 77.69 |
| secondary_c0.5 | B | 66.87 | 29.61 | 3.52 | 76.99 |
| primary_c0.4 | A | 73.23 | 23.79 | 2.97 | 77.69 |
| primary_c0.4 | B | 70.52 | 25.95 | 3.53 | 76.99 |
| primary_c0.6 | A | 78.03 | 19.93 | 2.04 | 77.69 |
| primary_c0.6 | B | 75.80 | 21.82 | 2.37 | 76.99 |

Medians over the carriers (primary c0.5), neither / one / both soft: A 77.36 / 18.29 / 0.16 %; B 74.31 / 20.12 / 0.19 %.

### 2c. Linewidth of the carriers by partner family pair (primary c0.5; %)

| partner pair | model A | model B |
|---|---|---|
| neither/neither | 76.13 | 73.71 |
| neither/pol | 2.62 | 2.91 |
| neither/rot | 18.74 | 20.46 |
| pol/pol | 0.00 | 0.00 |
| pol/rot | 0.06 | 0.08 |
| rot/rot | 2.45 | 2.84 |

Pairs absent from the table are zero in both models. pol/pol: A −1.6 × 10⁻⁴ %, B 3.1 × 10⁻⁴ % (shown as 0.00).
Negative bins: 3 per model, all pol/pol absorption parts of carriers; largest magnitude 8.9 × 10⁻⁵ (A) and
1.2 × 10⁻⁹ (B) of the carrier's width. The absorption term −(f1 − f2)(W1 − W2) is not sign-definite pair by pair
when the two partners are close in frequency (the tetrahedron weights W1, W2 overlap); docs/v4/ANALYSIS_STEPS.md, step T3 §1
corrected. Other definitions: `analysis/t4/t4_{A,B}_summary.json`.

### 2d. κ_P carried by each family, next to κ_P carried by modes with η > 0.1 (all modes; %)

| definition | model | neither | pol | rot | mixed | η > 0.1 (group-averaged) |
|---|---|---|---|---|---|---|
| primary_c0.5 | A | 85.96 | 0.45 | 13.59 | 0.00 | 0.00185 |
| primary_c0.5 | B | 86.01 | 0.51 | 13.49 | 0.00 | 0.00198 |
| secondary_c0.5 | A | 83.97 | 2.44 | 13.36 | 0.23 | 0.00185 |
| secondary_c0.5 | B | 83.79 | 2.73 | 13.26 | 0.23 | 0.00198 |
| primary_c0.4 | A | 79.88 | 0.68 | 19.44 | 0.00 | 0.00185 |
| primary_c0.4 | B | 79.77 | 0.81 | 19.42 | 0.00 | 0.00198 |
| primary_c0.6 | A | 89.33 | 0.25 | 10.42 | 0.00 | 0.00185 |
| primary_c0.6 | B | 89.49 | 0.28 | 10.23 | 0.00 | 0.00198 |

### 2e. Sensitivity to TRISYM, 8³

Same labels and settings, TRISYM = 1 against TRISYM = 0; each run's own 50 % carriers. Difference in percentage
points (TRISYM 1 − TRISYM 0).

| model | κ_P (W/mK) TRISYM 1 / 0 | carriers 1 / 0 (common) | neither / one / both soft, primary c0.5, TRISYM 1 (%) | same, TRISYM 0 (%) | largest abs. difference over all shares and definitions (pp) | at |
|---|---|---|---|---|---|---|
| A | 7.5376 / 7.5374 | 45 / 46 (44) | 71.24 / 25.93 / 2.83 | 70.91 / 26.10 / 2.99 | 0.326 | primary_c0.5: nsoft0 |
| B | 7.2799 / 7.2847 | 47 / 47 (46) | 67.78 / 28.38 / 3.84 | 68.27 / 27.99 / 3.74 | 0.570 | primary_c0.6: nsoft0 |

Escalation rule (> 2 pp): **not triggered**; the 16³ TRISYM = 0 run goes on the list of the analysis plan

