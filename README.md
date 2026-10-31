# Zone-resolved soft-mode damping and the thermal-conductivity decomposition of SrTiO₃

Computational study of cubic SrTiO₃ from 300 to 1000 K: self-consistent phonon (SCPH) renormalization of the
force constants distributed with the ALAMODE SrTiO₃ example, three-phonon linewidths, a mode-resolved damping map
η = Γ/Ω across the Brillouin zone, and the thermal conductivity split into particle (κ<sub>P</sub>) and coherence (κ<sub>C</sub>)
contributions. Both contributions are grouped by reference character under two mode-labelling rules. BaTiO₃ and
KTaO₃ appear only as literature context from the earlier releases.

Version 3.1.1 holds the calculations, tables and scripts of the revised manuscript. It also holds the
four-phonon diagnostics with the ALAMODE patches they used, the job manifest of every calculation, and the
rules and predictions recorded before the runs. The raw solver outputs are archived on Zenodo, in the same
record as this release.

## Key results (v3.1.1, revised manuscript)

All values are for this model: fixed-cell cubic SCPH on the inherited force constants, three-phonon scattering, and
the relaxation-time approximation for κ<sub>P</sub>.

- **Damping.** Most of κ<sub>P</sub> is carried by weakly damped modes. The largest three-phonon damping ratio on the tested
  meshes is η₃ = 0.184, for the antiferrodistortive R-point mode near room temperature. Modes with η₃ > 0.1 carry at
  most 0.14 % of κ<sub>P</sub> at any temperature.
- **Particle conductivity.** κ<sub>P</sub>(300 K) = 8.89 W m⁻¹ K⁻¹ on the 16³ reference mesh.
- **Coherence conductivity.** κ<sub>C</sub> converges more slowly than κ<sub>P</sub>. At 300 K it spans 0.44–0.54 W m⁻¹ K⁻¹ across
  the meshes from 15³ to 24³. This is an observed range, not an error bar.
- **Sectors.** Pairs assigned to the hard-optic reference sector carry most of κ<sub>C</sub> under both labelling rules:
  65–73 % (character rule) and 74–88 % (frequency rule) on meshes from 15³ up, 300–1000 K. Under the character rule
  this comes mainly from the number of pairs.
- **Benchmark.** At the published settings κ<sub>P</sub>(300 K) is 9.24 W m⁻¹ K⁻¹. The coarse SCPH grid of the distributed
  example (q₁ = 2³) gives a κ<sub>P</sub> 16–18 % higher than q₁ = 8³; production uses q₁ = 12³ (`docs/recorded_rules.md (R1)`). The
  R-point mode remains 22–26 % softer than the published table at q₁ ≥ 8³, for a reason not resolved here.
- **Low-temperature limit.** No converged cubic SCPH solution was obtained below 250 K under the stated criteria.
  250 K fails the q₁ rule, so the principal range starts at 300 K.
- **Four-phonon diagnostics.** These are selected-mode diagram-(c) calculations with a stated regularization and
  exclusion; the implementation checks are in SI Sec. S7. They are not a converged four-phonon conductivity or a
  whole-zone bound.

The earlier analyses are kept and documented in `RELEASE_NOTES.md`:
- v2.0.0: the distributed ALAMODE example outputs, in `data/processed/sto_alamode/`;
- v1.0.0: the literature-only damping classification and kernel.

The revised manuscript supersedes their SrTiO₃ conclusions where the two differ.

## Status

v3.1.1 accompanies the revised manuscript. Standing caveats and checks are in
`METHODS_AND_CAVEATS.md`.

## Repository structure

| Path | Contents |
| --- | --- |
| `README.md`, `pyproject.toml`, `uv.lock`, `.python-version` | project metadata; Python 3.11 with a locked environment |
| `METHODS_AND_CAVEATS.md` | consolidated audit trail (checks, caveats, deferred items) |
| `RELEASE_NOTES.md` | version history |
| `LICENSE`, `DATA_LICENSE.md` | MIT (code) and CC BY 4.0 (data and docs) |
| `data/raw/` | BTO/STO IR + BTO hybrid CSVs; `alamode_sto/` raw solver outputs are local (gitignored) and archived on Zenodo |
| `data/literature/` | digitized published data with per-source provenance: `BTO/`, `STO/`, `KTO/` |
| `data/processed/` | gitignored except the curated interfaces below |
| `data/processed/manuscript_numbers.csv` | traceability index of the v2 manuscript: every number → source script + data file |
| `data/processed/section4_eta_kappa_triples.csv` | three-material $(\eta, \kappa)$ classification triples |
| `data/processed/alamode_crosscheck/` | band-resolved τ(300 K) validation tables (bare and SCPH) |
| `data/processed/sto_alamode/` | 17 tidy computed SrTiO₃ tables: $\Omega_{\mathbf{q}s}(T)/\Gamma_{\mathbf{q}s}(T)$ (SCPH and bare), $\eta$ map + validity + sensitivity, mode-resolved κ and RTA totals, five-way decomposition, six-block coherence elements, mesh convergence (incl. isotope), $\Lambda^2(T)$ fits, dispersion paths, soft kernel, bare-vs-SCPH control, harmonic reference, Tadano Fig. 5 match |
| `src/fekappa/io.py`, `literature_io.py` | raw + literature data loaders (with provenance) |
| `src/fekappa/alamode_io.py` | ALAMODE RTA/coherence parsers + tidy-table emitters |
| `src/fekappa/damping_map.py` | $\eta$ map, validity fractions, sensitivity tables |
| `src/fekappa/kappa_decomposition.py` | mode-resolved $\kappa_P$, five-way sector decomposition, coherence blocks |
| `src/fekappa/lambda_fit.py` | small-q $\Lambda^2(T)$ fits, dispersion tables |
| `src/fekappa/soft_kernel_fp.py` | recomputed soft-mode kernel (computed + anchored variants) |
| `src/fekappa/scph.py` | SCPH fits, soft-mode renormalization |
| `src/fekappa/selfenergy.py` | 3-/4-phonon damping (self-energy) decomposition |
| `src/fekappa/eta.py` | damping ratio $\eta = \Gamma_{\mathrm{HWHM}}/\Omega$ (convention-harmonized) |
| `src/fekappa/kappa_kernel.py` | reduced long-wavelength soft-mode transport kernel $K_{\mathrm{soft}}$ |
| `src/fekappa/kappa.py` | Boltzmann/Wigner/total integrators (in progress) |
| `src/fekappa/plots.py` | publication matplotlib style |
| `notebooks/` | 01–09 (+05b, 07b): inventory → benchmarks → kernel |
| `paper/figures/` | manuscript figure scripts (`out/` renders gitignored) |
| `scripts/revision/` | revision analysis scripts: raw outputs → tables → manuscript numbers (`scripts/revision/README.md`) |
| `data/processed/sto_production/` | production tables (q₁ = 12³, polarization mixing): κ per T and mesh, mode table, damping map, sectors, coherence blocks |
| `data/processed/revision/` | revision tables with provenance headers (benchmark, convergence, labels, pair factors, four-phonon diagnostics) |
| `data/processed/si_v4/` | the tables of the Supplementary Information |
| `data/processed/manuscript_numbers_v8.csv` | every number of the revised main text, with its source table |
| `data/provenance/` | scrubbed job manifest: 170 entries, 169 lattice-dynamics runs and one entry recording the SCPH-grid selection rule (inputs, outputs, sha256, build ids) |
| `docs/recorded_rules.md` | dated rules and predictions with their outcomes |
| `queue/` | the checkpointed job runner, job specs and tests |
| `tools/alamode_patches/` | ALAMODE c5b015c patches, build notes and binary hashes (`BUILD.md`), one-quartet trace (`TRACE.md`) |
| `tools/zenodo/` | builder of the raw-output archives |
| `archive/scph_calibration/` | preserved legacy BTO calibration code (negative-result illustration only) |
| `tests/` | regression suite (unit conversions, parsers, committed-table locks, production-table closure, four-phonon identity count, queue) |

## Installation

Targets Python 3.11+ with a `src/` layout:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

## Data

**Revision calculations (v3).** The raw outputs of every job are in the Zenodo record of this release, one archive
per variant. Unpack them under `data/raw/alamode_sto/revision/` to rerun `scripts/revision/`. The job manifest
`data/provenance/manifest_final.json` gives the sha256 of every input and output file and the build id of the binary
that produced it. The runs made before the job queue existed are in the archive `prequeue_runs`; set
`FEKAPPA_RUNS` to the unpacked directory. Zenodo: concept DOI 10.5281/zenodo.20559513 (all versions).

**Literature layer.** Inputs are digitized, convention-harmonized literature data under `data/literature/`
(BaTiO₃, SrTiO₃, KTaO₃), loaded through `fekappa.literature_io` with per-source provenance.
Soft-mode frequencies Ω(T), linewidths Γ(T), dispersions, and inverse-permittivity data are
drawn from independent experimental and first-principles sources (e.g. Vogt hyper-Raman,
Yamada–Shirane and Harada–Axe–Shirane neutron, Fauqué, Verdi, Tadano–Tsuneyuki,
Schmidt–Spaldin, Meier, Ranalli). Linewidth conventions (FWHM vs HWHM), reciprocal-space
units, and temperature scales are normalized at the consumption layer before computing
η = Γ/Ω; see `METHODS_AND_CAVEATS.md` for the harmonization decisions and standing caveats.

A legacy hybrid BaTiO₃ dataset in `data/raw/BTO/` is a Cochran-calibrated phenomenological
reconstruction, retained **only** as a negative-result illustration (it cannot validate a
first-principles self-energy framework); the preserved calibration code lives under
`archive/scph_calibration/`.

## Reproducing the figures

Revised manuscript and SI (from tracked tables only):

```bash
pip install -e ".[dev]"
python scripts/revision/build_manuscript_numbers_v8.py   # data/processed/manuscript_numbers_v8.csv
(cd paper/figures && python make_v8_figures.py && python make_si_v4_figures.py)   # Figs. 1–4, S1–S8 in out/
pytest
```

The tables themselves regenerate from the raw outputs with the scripts in `scripts/revision/` (order and inputs in
its README).

Earlier releases: analysis notebooks live in `notebooks/` (data inventory and pre-flight gate → per-material
soft-mode benchmarks → the damping-ratio and transport-kernel synthesis). The
manuscript-ready figures are re-plotted from committed data / kernel outputs by standalone
scripts in `paper/figures/` — see `paper/figures/README.md`. Rendered images are written to
`paper/figures/out/` and are gitignored (regenerate them by running the scripts); only the
scripts and the curated data interfaces are tracked.

## License

- **Code** (everything under `src/`, `paper/figures/`, `notebooks/`, `tests/`) — MIT, see
  [`LICENSE`](LICENSE).
- **Data, figure data, and documentation** (generated/derived datasets, post-processed
  results, literature-extraction notes, and docs) — Creative Commons Attribution 4.0
  International (CC BY 4.0), see [`DATA_LICENSE.md`](DATA_LICENSE.md).

Third-party publications, figures, and source materials cited in the manuscript are **not**
covered by these licenses; literature-derived numerical anchors in this repository are
provided as derived extractions for reproducibility and should be cited through their
original sources.
