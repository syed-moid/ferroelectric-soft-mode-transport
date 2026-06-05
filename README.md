# Soft-mode damping and reduced thermal-transport kernels in anharmonic displacive ferroelectrics

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

## Status

Under review — the manuscript has been submitted to a journal, and this repository is being
prepared for public release alongside the submission.

## Repository structure

```
ferroelectric-soft-mode-transport/
├── README.md
├── pyproject.toml
├── .python-version                  # 3.11
├── METHODS_AND_CAVEATS.md            # consolidated audit trail (checks/caveats/deferred)
├── data/
│   ├── raw/                         # BTO/STO IR + BTO hybrid soft-mode CSVs (see note)
│   ├── literature/                  # digitized published data: BTO/, STO/, KTO/
│   └── processed/                   # gitignored, except the curated Section-4 anchor CSV
├── src/fekappa/
│   ├── io.py / literature_io.py     # raw + literature data loaders (with provenance)
│   ├── scph.py                      # SCPH fits, soft-mode renormalization
│   ├── selfenergy.py                # 3-/4-phonon damping (self-energy) decomposition
│   ├── eta.py                       # damping ratio η = Γ_HWHM/Ω (convention-harmonized)
│   ├── kappa_kernel.py              # reduced long-wavelength soft-mode transport kernel K_soft
│   ├── kappa.py                     # Boltzmann/Wigner/total integrators (in progress)
│   └── plots.py                     # publication matplotlib style
├── notebooks/                       # 01–09 (+05b, 07b): inventory → benchmarks → kernel
├── paper/figures/                   # manuscript figure scripts (out/ renders gitignored)
├── manuscript/                      # LaTeX sources
└── tests/
```

## Installation

Targets Python 3.11+ with a `src/` layout:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

## Data

Inputs are digitized, convention-harmonized literature data under `data/literature/`
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

Analysis notebooks live in `notebooks/` (data inventory and pre-flight gate → per-material
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
