# scripts/v4: from raw outputs to the tables and numbers of the v4 paper

Steps, scripts and runs: `docs/v4/ANALYSIS_STEPS.md`. Run names: `data/provenance/v4/run_map.csv`.
All inputs are command-line arguments or paths relative to the repository root. The files of the 2015 study are
arguments of the held scripts only; no script names them.

Order (repository root; `RUNS` = a directory with one link per run, named as in the run map):

```bash
python scripts/v4/regenerate_v4.py RUNS data/processed/v4 \
       --a3-transport data/processed/v4/a/benchmark/A3_transport.csv \
       --source-bands MODEL_B_DATASET/data/SrTiO3/figures/bands       # parts a and b, figures 3-5, S2, S3
python scripts/v4/figures/make_fig1_parity.py DFSET_TEST PRED_H1 PRED_H2 PRED_MODEL_A fig1.png
python scripts/v4/build_tables_v4.py data/processed/v4
python scripts/v4/figures/make_graphical_abstract_v3.py --without-archived --out OUTDIR   # graphical abstract, A and B
python scripts/v4/build_manuscript_numbers_v4.py data/processed/v4 data/processed/v4/a/benchmark/A3_transport.csv --runs RUNS
# held part (only with the files of the 2015 study):
python scripts/v4/regenerate_v4_held.py RUNS data/processed/v4 HELD --kappa-curve FILE --mode-file FILE \
       --arch-xml FILE --a-alm-in FILE --a-dfset FILE --b-dfset FILE --source-bands DIR
python scripts/v4/build_tables_v4.py data/processed/v4 --held HELD
python scripts/v4/build_manuscript_numbers_v4.py data/processed/v4 A3_TRANSPORT --runs RUNS --held HELD --held-runs RUNS
python scripts/v4/figures/make_graphical_abstract_v3.py --held HELD --out OUTDIR             # with the archived fit
```

The Markdown report writers (`t4_report.py`, `t5_report.py`, `t6_report.py`, `t7_report.py`, `t8_report.py`,
`p3_report.py`) write the tables that the SI tables were taken from; `build_tables_v4.py` turns each into a CSV.
Model inputs for the runs: `make_inputs_archived_fit.py` writes the inputs of the archived-fit runs from the
production inputs of model A.
