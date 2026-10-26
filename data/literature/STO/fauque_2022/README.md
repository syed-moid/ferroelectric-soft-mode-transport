# Fauqué et al. 2022 (arXiv:2203.15495): digitized data

Reference: B. Fauqué, P. Bourges, A. Subedi, K. Behnia, B. Baptiste, B. Roessli, T. Fennell, S. Raymond, and P. Steffens,
*Mesoscopic fluctuating domains in strontium titanate*, Phys. Rev. B 106, L140301 (2022),
https://doi.org/10.1103/PhysRevB.106.L140301 (arXiv:2203.15495). Title checked against Crossref on 2026-10-22; earlier notes
and the submitted manuscript gave it wrongly as "Mesoscopic tunneling in strontium titanate".

`dho_fits_T.csv` holds mode energies read from Fig. 1(d) (INS DHO fits, open circles = experimental; squares = TDEP). `tdep_dispersion_300K.csv` is from Fig. 1(c).

The CSV files carry no comment lines, because `src/fekappa/literature_io.py` reads them with a plain `pandas.read_csv` and checks the column schema. Provenance notes go here instead.

## Correction of 2026-10-15 (R point, `mode == "AFD"`, `method == "experimental"`)

| T_K | old omega_meV | new omega_meV | note |
|---|---|---|---|
| 200 | 3.5 | 4.69 | re-read centroid at T ≈ 203 K; T_K kept at the nominal 200 used for every other row at this temperature |
| 300 | 5.5 | 6.61 | the old value was the TDEP square at 300 K (still listed in its own `TDEP` row), not the experimental circle |

- Method: Fig. 1(d) rendered at 600 dpi and pixel-calibrated on the axis ticks; marker centroids located on the render.
- Centroid precision ±0.03 meV. The operative uncertainty is the marker half-size, ±0.62 meV (±5.0 cm⁻¹). The `omega_uncertainty_meV` column stays 0.0 for consistency with the other rows, which were not re-read; read it as "not recorded", not as exact.
- The other R rows were checked in the same reading and agree within 0.1 meV; they were not changed.

## Revision of 2026-10-22 (R point, `mode == "AFD"`, `method == "experimental"`)

| T_K | old omega_meV | new omega_meV | note |
|---|---|---|---|
| 600 | (no row) | 9.85 | experimental circle at 600 K, added; it was read on 2026-10-14 but not entered |
| 1000 | 11.7 | 11.6 | re-read centroid; the old value was an earlier, coarser reading |

- Source: Fauqué et al., Phys. Rev. B 106, L140301 (2022), Fig. 1(d), open circles (INS DHO fits). Read on 2026-10-14 from a
  600 dpi render, pixel-calibrated on the axis ticks, in the same reading as the 2026-10-15 correction above.
- Uncertainty: centroid precision ±0.03 meV; the operative uncertainty is the marker half-size, ±0.62 meV (±5.0 cm⁻¹).
- The 1000 K change (0.1 meV) is well inside that uncertainty. Earlier values are kept in the version history of this file and
  in the tables above.
- Loader check: `load_fauque2022_omega_T()` reads the file with the schema unchanged; the R rows go from 9 to 10.
