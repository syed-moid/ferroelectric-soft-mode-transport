# Manuscript figures

Standalone scripts that re-plot the manuscript figures from committed data / kernel
outputs. Each script reproduces its original analysis-notebook plotting cell **verbatim**
(same colors, markers, sizes, line styles, fonts via `fekappa.plots.set_publication_style`,
axis scales/limits) and applies only surgical edits — titles removed, repo-only prose boxes
removed, legends relocated only where they collided, panel labels for multi-panel figures.
Figure captions live in the manuscript, not in the images.

## Regenerate

```bash
source ../../.venv/bin/activate      # repo venv with `pip install -e ".[dev]"`
python make_f1_dispersion.py         # ... make_f2 ... make_f5 (revised manuscript)
python make_fsupp_sto.py             # revised supplementary S1-S3
python make_fig01.py                 # ... make_fig06 (v5 set; some survive as appendix figures)
python make_supp.py                  # v5 supplementary
```

## v3.0.0 figures (revised manuscript and SI)

| Script | Figures | Data dependency |
| --- | --- | --- |
| `make_v8_figures.py` | main-text Figs. 1–4 (`out/v8_fig{1..4}`) | `data/processed/sto_production/`, `data/processed/revision/exp_datasets.csv`, `data/literature/STO/` |
| `make_si_v4_figures.py` | SI Figs. S1–S8 (`out/si_fig{1..8}`) | `data/processed/si_v4/`, `data/processed/revision/`, `data/processed/sto_production/` |

Both read tracked tables only: `python make_v8_figures.py && python make_si_v4_figures.py` from this directory.

## v2 figures (computed SrTiO3 core of the earlier revision)

All read committed `data/processed/` (+ `data/literature/`) only — never `data/raw/`.

| Script | Figure | Data dependency |
| --- | --- | --- |
| `make_f1_dispersion.py` | F1 — SCPH vs bare dispersion; stabilization | `sto_alamode/dispersion_path_*.csv`, `mode_kappa_T.csv`, `eta_validity_T.csv`, Tadano Fig. 5 |
| `make_f2_eta_map.py` | F2 — zone-resolved damping map η(Ω, T) | `sto_alamode/eta_map_T.csv`, `eta_validity_T.csv`, `eta_validity_sensitivity.csv` |
| `make_f3_kappa_decomposition.py` | F3 — κ(T) five-way stack + negative control | `sto_alamode/kappa_five_way_T.csv`, `bare_vs_scph_control.csv`, κ literature anchors |
| `make_f4_classification.py` | F4 — damping-axis classification (+ computed STO) | `fekappa.eta`, HAS Table I, `sto_alamode/eta_validity_T.csv` |
| `make_f5_soft_kernel.py` | F5 — soft kernel: computed vs literature-input | `sto_alamode/soft_kernel_T.csv`, `kappa_five_way_T.csv`, `section4_eta_kappa_triples.csv` |
| `make_fsupp_sto.py` | S1–S3 — mesh convergence; Λ² fit; τ validation (one file each: `fs1_mesh`, `fs2_lambda`, `fs3_tau`) | `sto_alamode/*`, `alamode_crosscheck/*` |

Manuscript v6.1 numbering (by first callout): Figure 4 = `make_f5_soft_kernel.py`,
Figure 5 = `make_f4_classification.py`. Supplement: S4 = `make_fig01.py`,
S5 = `make_fig02.py`, S6–S7 = `make_supp.py` (LST checks), S8 = `make_fig05.py`;
the old v5 transport-plane supplement figure is no longer referenced.

Each script prints a manifest (outputs + the exact list of what was removed/moved/kept) and
writes `out/<name>.pdf` (vector) + `out/<name>.png` (600 dpi). Rendering is deterministic
(Agg backend); the same inputs produce byte-identical output.

## Files

| Script | Figure | Source cell | Data dependency |
| --- | --- | --- | --- |
| `make_fig01.py` | Fig 1 — ω(T) triangulation (STO / KTaO₃), stacked 2×1 | NB 07 + NB 07b | `literature_io` Ω(T) loaders |
| `make_fig02.py` | Fig 2 — STO Γ(T): YS vs Vogt (HWHM basis) | NB 07 | YS + Vogt γ(T) loaders |
| `make_fig03.py` | Fig 3 — BTO hybrid negative result (Ω²/δ, η) | NB 09 | `data/raw/BTO/phonon_modes_*K.csv` + `archive/scph_calibration` |
| `make_fig04.py` | Fig 4 — damping-axis classification (main) | NB 08 | `fekappa.eta` + HAS Table I |
| `make_fig05.py` | Fig 5 — cubic BTO directional Λ² | NB 06 | HAS 1971 dispersion + Λ² loaders |
| `make_fig06.py` | Fig 6 — soft-mode kernel κ_soft(T) | NB 08 | `fekappa.kappa_kernel` + Vogt Ω(T) |
| `make_supp.py` | S1–S3 — STO/KTaO₃ LST, η–κ plane | NB 07 / 07b / 08 | LST loaders + `data/processed/section4_eta_kappa_triples.csv` |
| `pubstyle.py` | — | — | thin helpers: path shim, margin panel label, dual PDF+PNG saver |

## Tracking

`out/` and all `paper/**/*.png` / `paper/**/*.pdf` are gitignored — figures regenerate from
the scripts. Only the scripts (and the committed `data/processed/section4_eta_kappa_triples.csv`
that Fig S3 reads) are tracked. Source publication PDFs are never committed.
