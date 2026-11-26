# Data and documentation license

Unless otherwise noted, the generated synthetic datasets, post-processed numerical results, figure data, configuration files, literature-extraction notes, and documentation in this repository are licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0).

License: https://creativecommons.org/licenses/by/4.0/

You are free to share and adapt these materials for any purpose, provided that appropriate credit is given to the author and the associated manuscript/repository is cited.

This license does not apply to third-party publications, figures, or source materials cited in the manuscript. Literature-derived numerical anchors included in this repository are provided as derived notes/extractions for reproducibility and should be cited through their original sources.

## Terms per part (v4)

| part | what | terms |
|---|---|---|
| a | model A and its refits: `data/processed/v4/a/`, `data/processed/v4/tables/a/`, `data/processed/v4/manuscript_numbers_v4_a.csv`; Zenodo archives `v4_runs_a`, `v4_processed_a` | CC BY 4.0 (as above) |
| b | material derived from the model B dataset (below) | CC BY-NC 4.0 |
| code | scripts, the `fekappa` package, tests, `tools/` including the patch | MIT (`LICENSE`) |

### Part b: material derived from the model B dataset

Source: C. Lin, J. Han, B. Xu and N. Marzari, Materials Cloud Archive 2026.91, DOI 10.24435/materialscloud:8b-q1,
licensed under the Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0); article
doi:10.1038/s41524-026-02163-1. Please cite the dataset and its article when using these files.

The following files are derived from that dataset and are licensed under CC BY-NC 4.0
(https://creativecommons.org/licenses/by-nc/4.0/). They are kept in their own folders and archives:
- `data/processed/v4/b/` and `data/processed/v4/tables/b/` (summaries, label-derived tables and report tables that
  contain model B values; `b/b4_model_comparison.csv` sets model A next to model B),
  `data/processed/v4/manuscript_numbers_v4_b.csv`;
- the Zenodo archives `v4_fits_b` (force data converted from the dataset, the fits H1 and H2, cross-validation,
  predicted forces on the test set, converted harmonic force constants and Born data), `v4_runs_b` (fitted and
  renormalized force constants and raw solver outputs of model B) and `v4_processed_b` (per-mode labels and mode
  tables of model B);
- the figures with model B values (Figs. 1, S1 and S3, and the model B series of Figs. 2, 3 and 5). The figure scripts
  are code (MIT).

What was changed: the force and displacement data of the dataset were converted to the input format of ALAMODE
(`DFSET_train`, `DFSET_test`); its DFPT harmonic force constants and Born data were converted to ALAMODE's formats;
anharmonic force constants of orders three to six were fitted to the converted training data (H1, with the harmonic
terms fixed; H2, with the harmonic terms fitted); the fitted force constants were used in self-consistent phonon and
relaxation-time calculations, whose outputs and summaries form the files above. The dataset's files themselves are
not redistributed.

### Material of the 2015 study

The production force constants, Born data and numerical results of the 2015 study (T. Tadano and S. Tsuneyuki, Phys.
Rev. B 92, 054301 (2015)) were shared privately by T. Tadano and were used in the calculations of the paper (the
benchmark and the fit test). They are not included in this release. The results obtained with them are given in the
paper and its Supplementary Information. Scripts that take those files as arguments are included. These scripts are
code (MIT); they contain no name and no value of those files.

### ALAMODE example files and the source patch

The model A inputs come from the SrTiO₃ example of ALAMODE (commit c5b015c), MIT License, Copyright (c) 2014
Terumasa Tadano; they are not redistributed in the archives (obtain them from ALAMODE at that commit).
`tools/alamode_patches/t3_partition.diff` modifies ALAMODE, carries that notice, and is released under the same MIT
License.
