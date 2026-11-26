# Ferroelectric-Kappa-SCPH — Reproducibility Package (v4.0.0)

Release for the paper "Polar and rotational soft modes as heat carriers and scattering partners in cubic SrTiO₃".
Code, tables and new raw-output archives form one new version of the Zenodo record (concept DOI
10.5281/zenodo.20559513).

## What is new

- **Second model.** Model B: a fit (H1) of orders three to six to the 300 K force data of Lin et al. (Materials
  Cloud Archive 2026.91, CC BY-NC 4.0) with the DFPT harmonic force constants of that dataset fixed; a companion fit
  (H2) with the harmonic terms fitted, used for frequencies only.
- **Family labels.** Polar, rotational and neither, from the atomic motion of each degenerate group, with a rule fixed
  before any result was resolved by family; four definitions (primary 0.4 / 0.5 / 0.6, unweighted polar 0.5).
- **Linewidth partition.** A patch of `anphon` (`tools/alamode_patches/t3_partition.diff`) records, per carrier, the
  three-phonon linewidth split by the labels of the two partners, without changing the standard outputs.
- **Analyses.** Carriers of half of κ<sub>P</sub> and their scattering partners (16³ and 20³, 300 K and 600 K); the comparison
  of rotational and non-rotational partners within frequency bands (105 / 145 / 160 cm⁻¹, and 80 / 95 cm⁻¹ below
  105); full against symmetry-reduced tetrahedron sum (8³, 16³); the penalty endpoints and the reduced-basis refit of
  model A.
- **Files of the 2015 study.** The production force constants, Born data and numerical results of the 2015 study (T.
  Tadano and S. Tsuneyuki, Phys. Rev. B 92, 054301 (2015)) were shared privately by T. Tadano and were used in the
  calculations of the paper (the benchmark and the fit test). They are not included in this release. The results
  obtained with them are given in the paper and its Supplementary Information. Scripts that take those files as
  arguments are included. Among them are the scripts for the archived 2015 fit, in the configuration of the 2015 study
  and under the settings of model A.
- **Scripts and tables.** `scripts/v4/` regenerates the processed outputs from the raw outputs (`regenerate_v4.py`;
  `regenerate_v4_held.py`, which needs the files of the 2015 study), the tables (`build_tables_v4.py`) and the number table
  (`build_manuscript_numbers_v4.py` → `data/processed/v4/manuscript_numbers_v4_{a,b}.csv`). Figure scripts in
  `scripts/v4/figures/`.
- **Provenance.** `data/provenance/v4/run_map.csv` names every run by what it does next to its job name;
  `manifest_v4.json` lists every file with its sha256 and the build id.
- **Package.** `fekappa.modes` (degenerate groups, group-averaged linewidths, mesh and cell from `rta.in`) and
  `fekappa.partition` (readers of the partition output and of `.result` linewidths), with tests. Version 4.0.0.
- **Tests.** 313 passed / 3 skipped with the local raw ALAMODE outputs; 308 passed / 8 skipped in a clean checkout
  without raw outputs (the five further tests need those outputs).

## What is superseded from v3.1.1

- The key results and the status text of the v3 README (zone-resolved damping and the κ decomposition of the
  v3 paper) are replaced by those of this paper.
- The v3 folders of scripts and tables are renamed `scripts/v3/` and `data/processed/v3/` (old and new paths in
  `docs/PATH_MAP_v3.md`); path references inside the files are updated, contents otherwise unchanged. The v3.1.1 archives and `data/provenance/manifest_final.json` keep the old names.
- `data/processed/manuscript_numbers_v8.csv` is the number table of the v3 paper; the v4 paper uses
  `data/processed/v4/manuscript_numbers_v4_{a,b}.csv`.

## What is unchanged

- The model A production dataset of v3 (SCPH on q₁ = 12³, RTA on 8³–24³, eigenvectors) and its raw-output archives;
  v4 reuses these runs (`run_map.csv`, archive "v3 (reused)").
- `data/processed/sto_production/`, `data/processed/v3/` (renamed), `data/processed/si_v4/`, `data/literature/`,
  `queue/`, the four-phonon patches and `docs/recorded_rules.md`.

## Raw outputs (new Zenodo archives of this version)

Archives `v4_runs_a`, `v4_runs_b`, `v4_fits_b`, `v4_processed_a` and `v4_processed_b`, each holding the run
directories or files under their run-map names; `manifest_v4.json` gives the sha256 of every file. In 21
`renorm.xml` files of `v4_runs_a` and `v4_runs_b` the home path of the compute node in the header is replaced by
`<home>`; the manifest records both hashes (`sha256` as produced, `sha256_archived`).

# Ferroelectric-Kappa-SCPH — Reproducibility Package (v3.1.1)

Removes the journal name from the Zenodo packaging script (`tools/zenodo/make_release_v3.py`: the README template and
the metadata notes). Code, tables and data are those of v3.1.0. Package version 3.1.1.

# Ferroelectric-Kappa-SCPH — Reproducibility Package (v3.1.0)

Corrects a mis-attributed experimental value: the `tachibana_2008` entry (11 ± 1 W m⁻¹ K⁻¹ at 300 K) is removed,
because Appl. Phys. Lett. 93, 092902 (2008) does not report SrTiO₃. Figure 3 now shows κ<sub>P</sub> + κ<sub>C</sub> on the 16³ mesh.
Adds the SCPH force-constant-correction reproduction check and the single-rounding check for the four-phonon ratio
tables.

- Removed: `data/literature/STO/tachibana_2008/`, its loader in `fekappa.literature_io` and its tests; the Tachibana
  rows of `exp_datasets.csv` and of `manuscript_numbers_v8.csv`.
- Added: `scripts/v3/dfc2_reproduction_check.py` → `data/processed/v3/dfc2_reproduction.csv`;
  `scripts/v3/s7_ratios_single_rounding.py` → `data/processed/v3/s7_ratios_unrounded.csv`.
- The raw-output archives are those of v3.0.0, unchanged.
- Package version 3.1.0.

# Ferroelectric-Kappa-SCPH — Reproducibility Package (v3.0.0)

Release for the v3 manuscript on zone-resolved soft-mode damping and the particle/coherence decomposition of
the thermal conductivity of cubic SrTiO₃ (300–1000 K). Code, tables and the raw outputs (Zenodo archives) form one
new version of the Zenodo record (concept DOI 10.5281/zenodo.20559513), so the paper cites one version DOI.

## Summary

- **Production dataset.** Own SCPH on q₁ = 12³ with polarization mixing, on the force constants of the ALAMODE
  SrTiO₃ example (commit c5b015c). Wigner-RTA transport on 8³, 12³ and 16³ from 300 to 1000 K in 50 K steps, plus 15³,
  17³, 20³ and 24³ at selected temperatures. Every job is in the scrubbed manifest
  `data/provenance/manifest_final.json`: 170 jobs, each file with its sha256 and each job with the build id of its
  binary.
- **Convergence.**
  - Transport-mesh sweeps for κ<sub>P</sub> and κ<sub>C</sub>, and the mesh diagnostic r.
  - Two κ<sub>C</sub> predictions recorded before their runs: the 24³ prediction at 300 K failed, and the 20³ prediction at
    1000 K held (R2, R8).
  - The SCPH-grid check q₁ = 12³ against 16³ at 250, 300, 500 and 800 K (`n3_q1_16_check.csv`).
- **Benchmark.**
  - Matrix of SCPH grid × polarization mixing × transport mesh against the published calculation (SI S2).
  - R-point mode comparison, with all benchmark percentages computed once from unrounded values (SI S7b).
- **Mode labels.** Reference-character labels under two rules (character and frequency). Exactly degenerate
  subspaces are labelled as a whole (R10). Pair counts are given per irreducible record and star-weighted
  (`e1_pair_counts.csv`).
- **Four-phonon diagnostics (diagram c, selected modes).**
  - The identity-set exclusion patch against ALAMODE c5b015c and its 1e-6 Ry tolerance variant.
  - The trace build and three validation builds: omitted sum, strict-index exclusion and full-precision dump.
  - The worked single-quartet check (`tools/alamode_patches/TRACE.md`) and the build notes with full binary
    hashes (`BUILD.md`).
  - The full-precision printed inputs of SI Table S21 (`s21_tolerance_printed.txt`).
  - These are implementation checks of a stated finite-mesh prescription, not a converged four-phonon conductivity.
- **Recorded rules.** `docs/recorded_rules.md` gives each rule and prediction (R1–R11) with the UTC time it was
  written, the runs it governs, its outcome and the SI section that reports it.
- **Scripts and tables.** `scripts/v3/` regenerates every table from the raw outputs. Every table in
  `data/processed/v3/` and `data/processed/si_v4/` carries a provenance header (`fekappa.provenance`). The
  main-text numbers (`manuscript_numbers_v8.csv`) and Figs. 1–4 and S1–S8 regenerate from tracked tables alone.
- **Compute queue.** `queue/` holds the checkpointed runner (stdlib only), the job specs and its unit tests.
- **Package.**
  - `fekappa.evec_io` (anphon eigenvector reader) and `fekappa.provenance` are new.
  - Version 3.0.0; lock file refreshed.
  - Tests: **311 passed / 2 skipped** with the raw outputs present; **305 passed / 8 skipped** in a clean checkout
    (the raw-dataset tests skip).
- **Literature corrections.** Fauqué et al. 2022: R-point rows re-read and the paper title corrected
  (`data/literature/STO/fauque_2022/README.md`).
- **Housekeeping.** Notebook outputs stripped.

## Raw outputs (Zenodo archives of this version)

| archive | contents |
|---|---|
| `scph_q12` | production SCPH solution (q₁ = 12³, mixing), its dispersion and the temperature-dependent FC2 |
| `pm_q12_rta` | production RTA runs: κ<sub>P</sub>, κ<sub>C</sub> and the per-pair coherence elements at every T and mesh |
| `pm_q12_evec` | eigenvectors on each transport mesh and temperature (sector labels) |
| `scph_q16` | SCPH on q₁ = 16³ and the four RTA runs on it (250, 300, 500, 800 K) |
| `pm_q8` | SCPH on q₁ = 8³ with its RTA runs, and the q₁ = 8³ run with reduced mixing |
| `pm_q12_4ph_audit` | all per-mode four-phonon runs: unmodified pilot, rule, tolerance, trace and validation builds |
| `prequeue_runs` | runs made before the queue: SCPH/RTA benchmark matrix, settings sweep, bare-harmonic and `|ω|` runs, eigenvector sets, harmonic cross-check |

`SHA256SUMS` covers every file. The two timing probes of the four-phonon runs are not archived; they are listed
in the manifest as excluded.

---

# Ferroelectric-Kappa-SCPH — Reproducibility Package (v2.0.0)

Release of the manuscript as a computational study
of SrTiO₃ (SCPH-coupled anharmonic lattice dynamics; BaTiO₃/KTaO₃ as literature-benchmarked
context). Archived as a **new version of the existing Zenodo record**
(concept DOI lineage 10.5281/zenodo.20559513) so the version history stays linked.

## New in v2

- Computed SrTiO₃ dataset: SCPH-coupled RTA at 0–1000 K (21 temperatures, 8×8×8 mesh,
  525 modes each), with the Wigner coherence channel computed directly
  (`data/processed/sto_alamode/`, provenance headers on every table).
- Zone-resolved damping map η<sub>qs</sub>(T) with count- and κ-weighted validity fractions and a
  Vogt-rescaling sensitivity table (`eta_map_T.csv`, `eta_validity_T.csv`,
  `eta_validity_sensitivity.csv`).
- Mode-resolved κ<sub>P</sub> (closure vs the solver κ ≤ 0.3 %) and the five-way decomposition with
  partition-boundary sensitivity (`mode_kappa_T.csv`, `kappa_five_way_T.csv`).
- Bare-RTA vs SCPH-coupled negative-control table (`bare_vs_scph_control.csv`).
- Small-q Λ²(T) fits from per-temperature renormalized force constants
  (`lambda_sq_T.csv`), full-path dispersions (`dispersion_path_*.csv`), and the
  Tadano-variant match table (`tadano_fig5_match.csv`).
- Recomputed soft-mode kernel, computed and Vogt-anchored variants (`soft_kernel_T.csv`).
- Mesh-convergence table 8³ vs 12³ at 300 K (`mesh_convergence.csv`).
- New analysis modules (`alamode_io`, `damping_map`, `kappa_decomposition`, `lambda_fit`,
  `soft_kernel_fp`) with tests; suite now **297 passed / 2 skipped**.
- New figure scripts `make_f1_dispersion.py` … `make_f5_soft_kernel.py`,
  `make_fsupp_sto.py` (each figure regenerates from committed processed tables only).
- Literature layer: κ benchmark anchors added (Tadano 2015 digitized 300 K point;
  Tachibana 2008 experimental anchor); new standing caveat C8 (fixed-cell SCPH Ω₀
  systematic) in `METHODS_AND_CAVEATS.md`.
- Reference list of the manuscript audited against publisher records; several v5 entries
  corrected (Fauqué, Verdi, Ranalli, Ehsan, Meier, Masuki, Kim, Hlinka, Tachibana,
  Schmidt–Spaldin).

Raw solver outputs (`data/raw/alamode_sto/`) are local-only and regenerable from the
documented inputs and logs (`IMPORT_NOTE.md`, `PROVENANCE.md`).

---

# Ferroelectric-Kappa-SCPH — Reproducibility Package (v1.0.0)

Initial release of the code, convention-harmonized literature data, and figure pipeline
for the manuscript:

> **Soft-mode damping and reduced thermal-transport kernels in anharmonic displacive
> ferroelectrics** *(under review).*

## Overview

This work introduces a Green's-function-based soft-mode transport analysis for anharmonic
displacive ferroelectric perovskites, with specific emphasis on the regime where polar soft
modes are strongly renormalized and may cease to behave as well-defined quasiparticles. The
novelty is not a full first-principles prediction of total lattice thermal conductivity;
rather, it is a controlled, sector-resolved treatment of the long-wavelength soft-mode
contribution using self-consistent phonon renormalization, phonon self-energies, spectral
functions, damping ratios, and the Kubo–Wigner transport framework.

A central methodological contribution is the construction of a provenance-gated and
convention-harmonized literature dataset for BaTiO₃, SrTiO₃, and KTaO₃. The analysis shows
that linewidth definitions, reciprocal-space units, temperature-label conventions, and
digitization choices can qualitatively affect the inferred damping ratio. By enforcing a
consistent half-linewidth convention, the work establishes a robust damping-axis
classification: SrTiO₃ and KTaO₃ form an underdamped incipient-ferroelectric cluster,
whereas cubic BaTiO₃ lies in a distinct overdamped soft-mode regime.

The work also identifies and preserves a negative validation result for a previously
assembled BaTiO₃ hybrid dataset. Rather than treating this as a failure of the theoretical
framework, it shows that phenomenologically reconstructed frequencies and model-fitted
linewidths cannot validate a first-principles Green's-function/self-energy theory. This
clarifies the data requirements for future full thermal-conductivity calculations and
motivates the reduced soft-mode kernel developed here.

Overall, the study provides a reproducible bridge between soft-mode lattice dynamics and
thermal transport. It defines what can be predicted from currently available literature
data, what cannot be claimed without new anharmonic force-constant calculations, and how
future first-principles Wigner conductivity calculations for ferroelectric perovskites
should be benchmarked.

## Key highlights

- A finite-temperature Green's-function / SCPH / Wigner framework is formulated for
  soft-mode thermal transport in anharmonic displacive ferroelectric perovskites.
- A provenance-gated literature dataset is constructed for BaTiO₃, SrTiO₃, and KTaO₃, with
  explicit correction of linewidth conventions, reciprocal-space units, temperature-label
  issues, and digitization bias.
- A previously assembled BaTiO₃ hybrid dataset is rejected as an independent validation
  target because its frequencies and linewidths are not independent first-principles
  self-energy outputs.
- The damping ratio η = Γ/Ω separates the three materials into two regimes: SrTiO₃ and
  KTaO₃ are underdamped, while cubic BaTiO₃ is strongly overdamped.
- The reduced long-wavelength soft-mode transport kernel is of order 10⁻¹ W m⁻¹ K⁻¹ for
  SrTiO₃ and KTaO₃, confirming that the soft-mode sector is a partial transport weight
  rather than the total thermal conductivity.
- The work defines the additional first-principles inputs required for a full total-κ
  calculation: Brillouin-zone-wide third- and fourth-order force constants,
  SCPH-renormalized eigenvectors, linewidths, spectral functions, and Wigner velocity
  matrices.

## Results at a glance

| # | Result | Value (source) |
|---|--------|----------------|
| 1 | Damping cluster — two underdamped QPE | η<sub>STO</sub> ≈ 0.109, η<sub>KTaO3</sub> ≈ 0.126 ( \|Δ\| ≈ 0.016 ), both < 0.5 (NB 08 / Fig 4) |
| 2 | Overdamped ferroelectric | η<sub>BaTiO3</sub> ≈ 2.27 (HAS 1971 cubic); η-axis interior (0.5–2.0) empty |
| 3 | Soft-mode kernel ≪ total κ | κ<sub>soft</sub>: STO ≈ 0.137 [0.102, 0.167], KTaO₃ ≈ 0.118 [0.065, 0.157] W m⁻¹ K⁻¹ (NB 08 / Fig 6) |
| 4 | BaTiO₃ kernel order-of-magnitude only | soft-[100] vs stiff bracket ≈ 0.005–0.019 W m⁻¹ K⁻¹ (overdamped; quasiparticle picture invalid) |
| 5 | BaTiO₃ directional curvature | Λ²_soft[100] ≈ 972, Λ²_stiff[110]/[111] ≈ 4750 meV² Å² (anisotropy ≈ 5×) (NB 06 / Fig 5) |
| 6 | SrTiO₃ effective curvature | Λ² ≈ 889 meV² Å² (Verdi Γ→X-TO1 effective branch curvature, ±30%; KTaO₃ STO-proxy ±50%) |
| 7 | LST consistency | SrTiO₃ (Verdi) R² ≈ 0.998; KTaO₃ (Ranalli, internal) R² ≈ 0.993 (S1 / S2) |
| 8 | Rejected hybrid negative result | scalar-SCPH residual δ = +2806% at 400 K; hybrid η crosses 0.5 at 400 K (NB 09 / Fig 3) |

Convention basis: η on a half-linewidth (HWHM) basis; FWHM-reported widths are halved once,
at the consumption layer. Numbers are the committed benchmark values — see
`METHODS_AND_CAVEATS.md` for provenance and uncertainty.

## What's in this release

- **`src/fekappa/`** — loaders (`io`, `literature_io`), `scph`, `selfenergy`, `eta`,
  `kappa_kernel`, `kappa` (integrators in progress), `plots`.
- **`data/literature/`** — digitized, convention-harmonized published data with provenance;
  **`data/processed/section4_eta_kappa_triples.csv`** — the curated anchor table.
- **`notebooks/`** — 01–09 (+05b, 07b).
- **`paper/figures/`** — figure scripts (see `paper/figures/README.md`); rendered images are
  gitignored and regenerate from the scripts.
- **`tests/`** — unit tests for `io`, `scph`, `selfenergy`, `eta`, `kappa_kernel`,
  `literature_io`, and the literature-input schemas.
- **`METHODS_AND_CAVEATS.md`**, **`archive/scph_calibration/`** (preserved aborted
  scalar-SCPH calibration reproduced by NB 09).
- **Not included:** copyrighted source publication PDFs (excluded by design); rendered figure
  binaries (regenerate from scripts).

## Reproducing the results

CPU-only; each step runs in minutes. From the repo root with the environment installed
(`python3.11 -m venv .venv && source .venv/bin/activate && pip install -e ".[dev]"`):

```bash
pytest -q                                            # 261 passed / 3 skipped
jupyter nbconvert --to notebook --execute notebooks/08_three_material_cluster.ipynb   # cluster + kernel (Results 1–4, 6)
jupyter nbconvert --to notebook --execute notebooks/06_BTO_lambda_refit.ipynb          # BaTiO3 Λ² (Result 5)
jupyter nbconvert --to notebook --execute notebooks/09_BTO_hybrid_appendix.ipynb       # negative result (Result 8)
python paper/figures/make_fig04.py                   # regenerate any figure (Figs 1–6, S1–S3)
```

Each notebook prints a gate summary; each figure script prints a manifest of what it
removed/moved relative to the source cell.

## Requirements

- **Python 3.11+** (`.python-version` pins 3.11). Dependencies declared in `pyproject.toml`
  (`numpy`, `scipy`, `pandas`, `matplotlib`, `lmfit`; `pytest` + `ruff` for dev) — install
  with `pip install -e ".[dev]"`.
- **CPU only.** No GPU; modest memory; no large downloads (the literature anchors are small
  CSV/JSON files in-repo).

## Citation

If you use this package, please cite the repository and the manuscript (once published). A
release DOI is minted by Zenodo on archive — please cite that DOI. See `README.md` for the
BibTeX entry.
