# Zone-resolved soft-mode damping and the thermal-conductivity decomposition of SrTiO₃

This is a computational study of SrTiO₃ based on self-consistent-phonon (SCPH)-coupled
anharmonic lattice dynamics, with BaTiO₃ and KTaO₃ retained as literature-benchmarked
context materials. The repository computes, validates, and releases: the SCPH-renormalized
frequencies $\Omega_{\mathbf{q}s}(T)$ and anharmonic linewidths $\Gamma_{\mathbf{q}s}(T)$
of all 525 mesh modes at 0–1000 K; the zone-resolved damping map
$\eta_{\mathbf{q}s} = \Gamma/\Omega$ with a quantitative validity boundary for the
phonon quasiparticle picture (count- and κ-weighted); the mode-resolved thermal
conductivity and its five-way decomposition (acoustic / soft-optic / hard-optic particle
channels + the directly computed Wigner coherence channel); a controlled bare-RTA vs
SCPH-coupled negative result; the small-q soft-branch stiffness Λ²(T); and the recomputed
soft-mode transport kernel with both computed and experiment-anchored input variants.

The earlier literature-only analysis (convention-harmonized damping classification and
reduced kernel for BaTiO₃/SrTiO₃/KTaO₃) is retained: its provenance-gated literature
dataset now serves as the benchmark layer for the computed SrTiO₃ results, and its
three-material damping-axis classification is anchored at SrTiO₃ by first-principles data.

## Key results

- **Quasiparticle validity, zone-wide:** no SrTiO₃ mode exceeds η = 0.5 at any
  temperature 0–1000 K (zone max η ≈ 0.10; median 0.003–0.021), by mode count or by
  κ-weight — robust against the disclosed fixed-cell SCPH soft-mode frequency systematic
  (explicit rescaling sensitivity test included).
- **κ(300 K) = 9.53 W m⁻¹ K⁻¹** on the 8×8×8 mesh (10.31 on 12×12×12), vs 8.7 published
  first-principles (12³) and ≈ 11 experimental; mode-resolved particle sum closes on the
  solver κ to ≤ 0.3 % at every temperature.
- **Wigner coherence channel computed directly:** $\kappa_C/\kappa$ grows from 0.8 % (100 K) to
  5.5 % (1000 K); mode-pair-resolved elements released.
- **Controlled negative result:** identical anharmonic force constants without SCPH
  coupling → 11/525 unstable modes, κ low by 2.6–4.7×, low-frequency band lifetimes up to
  ~9× short while stiff optic bands are unaffected.
- **Λ²(T) computed:** 1228 → 690 meV² Å² over 0–400 K (R² ≥ 0.999); the 300 K value
  (842) validates the previous effective literature curvature (889). Above ≈ 450 K the
  TO1–TO2 zone-center anticrossing ends the single-branch fit regime (documented).
- **Soft-mode kernel recomputed:** 0.32 [0.21, 0.43] W m⁻¹ K⁻¹ at 300 K vs 0.137
  [0.102, 0.167] literature-only; computed and Vogt-anchored variants nearly coincide,
  quantifying how the fixed-cell SCPH systematic propagates (largely cancels) in transport.
- **Damping-axis classification retained and anchored:** SrTiO₃/KTaO₃ underdamped cluster
  (η ≈ 0.109/0.126, literature) vs overdamped BaTiO₃ (η ≈ 2.27), with the computed SrTiO₃
  zone-center point (η₀ = 0.0225) shown alongside and the offset attributed (fixed-cell
  SCPH Ω₀ provenance).

Provenance note: the anharmonic force constants underlying the SrTiO₃ dataset are the
published, validated ALAMODE example set (disclosed in Methods and in
`data/raw/alamode_sto/IMPORT_NOTE.md`); every derived quantity is benchmarked against
independent published references, and regeneration from an in-house PBEsol pipeline is
flagged as future work. The known fixed-cell SCPH overestimate of the soft-mode frequency
(1.9–2.4× the hyper-Raman values below 300 K) is disclosed, benchmarked, and propagated
through sensitivity analyses — never adjusted away (`METHODS_AND_CAVEATS.md`, C8).

## Status

In revision — an earlier version was reviewed at a journal; this repository holds the
revised, computation-centered version being prepared for resubmission.

## Repository structure

| Path | Contents |
| --- | --- |
| `README.md`, `pyproject.toml`, `uv.lock`, `.python-version` | project metadata; Python 3.11 with a locked environment |
| `METHODS_AND_CAVEATS.md` | consolidated audit trail (checks, caveats, deferred items) |
| `RELEASE_NOTES.md` | version history |
| `LICENSE`, `DATA_LICENSE.md` | MIT (code) and CC BY 4.0 (data and docs) |
| `data/raw/` | BTO/STO IR + BTO hybrid CSVs; `alamode_sto/` inputs are local-only (gitignored) |
| `data/literature/` | digitized published data with per-source provenance: `BTO/`, `STO/`, `KTO/` |
| `data/processed/` | gitignored except the curated interfaces below |
| `data/processed/manuscript_numbers.csv` | traceability index: every manuscript number → source script + data file |
| `data/processed/section4_eta_kappa_triples.csv` | three-material $(\eta, \kappa)$ classification triples |
| `data/processed/alamode_crosscheck/` | band-resolved τ(300 K) validation tables (bare and SCPH) |
| `data/processed/sto_alamode/` | 17 tidy computed SrTiO₃ tables: $\Omega_{\mathbf{q}s}(T)$/$\Gamma_{\mathbf{q}s}(T)$ (SCPH and bare), $\eta$ map + validity + sensitivity, mode-resolved κ and RTA totals, five-way decomposition, six-block coherence elements, mesh convergence (incl. isotope), $\Lambda^2(T)$ fits, dispersion paths, soft kernel, bare-vs-SCPH control, harmonic reference, Tadano Fig. 5 match |
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
| `tests/` | regression suite (unit conversions, parsers, committed-table locks) |

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
