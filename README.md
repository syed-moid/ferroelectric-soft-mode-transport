<!-- DRAFT for the author. Facts from the manuscript with their source tables; the author writes and edits the text.
[HELD] marks scripts that run on the files of the 2015 study when these are given as arguments; those files, and
everything derived from them, are not part of this release. -->

# Polar and rotational soft modes as heat carriers and scattering partners in cubic SrTiO₃

Computational study of cubic SrTiO₃ at 300 K (and 600 K): self-consistent phonon (SCPH) frequencies and
eigenvectors, three-phonon linewidths in the relaxation-time approximation, and the particle (κ<sub>P</sub>) and coherence
(κ<sub>C</sub>) conductivities, for two anharmonic models fitted to different force data. Modes are labelled polar,
rotational or neither by their atomic motion, with a rule fixed before any result was resolved by family. The
three-phonon linewidth of the modes that carry half of κ<sub>P</sub> is split by the labels of the two scattering partners
(patched `anphon`, `tools/alamode_patches/t3_partition.diff`).

- **Model A**: the SrTiO₃ example distributed with ALAMODE (commit c5b015c), force constants as distributed.
- **Model B**: a fit (H1) to the independent 300 K force data of the Materials Cloud Archive record 2026.91.
- **[HELD]** The archived 2015 fit (production force constants of the 2015 study, provided by its author).

Version 4.0.0 holds the scripts, tables and number table of this paper. Earlier releases (v1–v3.1.1) are described
in `RELEASE_NOTES.md`; their files are kept. Two folders of v3 have new names in v4 (`docs/PATH_MAP_v3.md`).

## Key results (v4.0.0)

300 K, 16³ transport mesh, q₁ = 12³ SCPH grid, primary family rule (threshold 0.5), group-averaged linewidths.

| quantity | model A | model B | source |
|---|---|---|---|
| κ<sub>P</sub> (W m⁻¹ K⁻¹) | 8.891 | 8.685 | Table 3 |
| κ<sub>C</sub> (W m⁻¹ K⁻¹) | 0.508 | 0.513 | Table 3 |
| share of κ<sub>P</sub> carried by rotational modes (%) | 13.59 | 13.49 | Table 5 |
| share of κ<sub>P</sub> carried by polar modes (%) | 0.45 | 0.51 | Table 5 |
| carriers (irreducible modes) carrying half of κ<sub>P</sub> | 203 | 211 | Table 5 |
| carrier-weighted linewidth with no soft partner (%) | 76.13 | 73.71 | Table 5 |
| … with one soft partner / both partners soft (%) | 21.36 / 2.51 | 23.37 / 2.91 | Table 5 |
| … partner pair neither / rotational (%) | 18.74 | 20.46 | Table 5 |
| per-mode ratio rotational / non-rotational partners, below 105 cm⁻¹ | 3.69 | 4.52 | Table 6 |
| same ratio below 80 cm⁻¹ | 5.68 | 6.32 | Table 7 |
| same ratio 80–95 / 95–105 cm⁻¹ | 1.80 / 1.53 | 1.99 / 1.90 | Table 7 |
| share of κ<sub>P</sub> carried by modes with η > 0.1 (%) | 0.001846 | 0.001982 | Table 4 |
| largest damping ratio η (R triplet) | 0.1284 | 0.1581 | Table 4 |
| κ<sub>P</sub> on 20³ (W m⁻¹ K⁻¹) | 9.209 | 9.012 | Table 5 |

- Full sum against the symmetry-reduced sum (16³): κ<sub>P</sub> changes by +0.003 % (A) and −0.018 % (B) (Table S31).
- 600 K: the rotational share of κ<sub>P</sub> is 12.29 % (A) and 12.18 % (B); the no-soft-partner linewidth share is
  85.64 % and 84.25 % (Sec. 5.3; SI S6). Model B at 600 K is outside its training distribution.

These are values of the stated model: three-phonon scattering, relaxation-time approximation, 2 × 2 × 2 supercell
force constants. Caveats are in `METHODS_AND_CAVEATS.md`.

## Repository structure (v4 additions)

| Path | Contents |
| --- | --- |
| `scripts/v4/` | analysis scripts of this paper: raw outputs → processed outputs → tables → number table (`scripts/v4/README.md`) |
| `scripts/v4/figures/` | figure scripts: Fig. 1 (`make_fig1_parity.py`), Fig. 3 and the soft-mode draft (`make_figure_drafts.py`), Figs. 4, 5, S2, S3 (`make_t4_figures.py`); graphical abstract of models A and B (`make_graphical_abstract_v3.py --without-archived`); [HELD] Fig. 2 (`make_fig2_softmode_grid.py`), Fig. S4 (`make_kappa_T_2015.py`), the graphical abstract with the archived fit |
| `src/fekappa/modes.py`, `src/fekappa/partition.py` | degenerate groups and group-averaged linewidths; readers of the linewidth-partition output |
| `data/processed/v4/a/` | model A and its refits: summaries of the carrier, band and damping-ratio analyses, benchmark tables |
| `data/processed/v4/b/` | model B, and tables with A and B side by side (see the data licence below) |
| `data/processed/v4/tables/` | the report tables behind the SI tables, one CSV per table with a provenance header (`INDEX.csv`) |
| `data/processed/v4/manuscript_numbers_v4_{a,b}.csv` | every number of the paper with its source file and run (parts a and b) |
| `data/provenance/v4/` | run map (every run by what it does, next to its job name) and scrubbed manifest (files, sha256, build ids) |
| `docs/v4/ANALYSIS_STEPS.md` | what each analysis step computes, with its script and runs |
| `tools/alamode_patches/t3_partition.diff` | linewidth-partition patch of `anphon` (build notes in `BUILD.md`) |

The tables of the earlier releases (`data/processed/sto_production/`, `data/processed/v3/`,
`data/processed/si_v4/`, …) are unchanged; some are inputs of v4 (`manuscript_numbers_v4_a.csv`, column `source`).

## Data

- **Raw outputs.** The run directories of v4 are in new archives of the Zenodo record (concept DOI
  10.5281/zenodo.20559513), named in `data/provenance/v4/run_map.csv` (column `archive`); the model A production
  runs are those of v3 (unchanged archives). Unpack them and link each run directory under its run-map name, then
  `python scripts/v4/regenerate_v4.py RUNS_ROOT data/processed/v4` regenerates parts a and b.
- **Per-mode outputs** (family labels on the full mesh, mode tables of the carrier analysis) are in the archives
  `v4_processed_a` and `v4_processed_b`; the repository holds the summaries.
- **Model B** uses the dataset of Lin et al., Materials Cloud Archive 2026.91 (DOI 10.24435/materialscloud:8b-q1),
  licensed CC BY-NC 4.0. Files derived from it are listed in `DATA_LICENSE.md`.
- **[HELD]** The files of the 2015 study, and the outputs derived from them, are not part of this release.

## Reproducing

```bash
pip install -e ".[dev]"
pytest
python scripts/v4/regenerate_v4.py RUNS_ROOT data/processed/v4 --a3-transport data/processed/v4/a/benchmark/A3_transport.csv \
       --source-bands MODEL_B_DATASET/data/SrTiO3/figures/bands
python scripts/v4/build_tables_v4.py data/processed/v4
python scripts/v4/build_manuscript_numbers_v4.py data/processed/v4 data/processed/v4/a/benchmark/A3_transport.csv --runs RUNS_ROOT
```

## License

- **Code**: MIT (`LICENSE`). The patch in `tools/alamode_patches/` modifies ALAMODE (MIT, © 2014 Terumasa Tadano) and
  is released under the same MIT License.
- **Data and documentation**: CC BY 4.0, except the files derived from the model B dataset, which are in their own
  folders and under CC BY-NC 4.0 (`DATA_LICENSE.md`).
