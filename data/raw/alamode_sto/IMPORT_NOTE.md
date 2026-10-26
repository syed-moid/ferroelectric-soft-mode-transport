# Import note — ALAMODE SrTiO3 SCPH-coupled RTA dataset

Copied 2026-08-18 from the companion lattice-dynamics repository (read-only
source; nothing is modified here):

- Origin directory:
  `data/raw/alamode_sto/` of the companion lattice-shear-viscosity repository
  - `STO_RTA_scph_{100,150,200,250,300,350,400}K.result` / `.kl`
    (SCPH-coupled anphon RTA, one run per temperature, 8x8x8 mesh)
  - `STO_RTA_production.result` / `.kl` (bare-RTA, 0-1000 K, 21 temperatures)
  - `RTA_production.in/.log`, `RTA_scph_*K.in/.log` (exact inputs + run logs)
  - `PROVENANCE.md` (full build recipe, calculation record, validation)
- Origin directory:
  `data/processed/reports/` of the companion lattice-shear-viscosity repository
  - `alamode_sto_tau300K_bands_{bare,scph}.csv` -> copied to
    `data/processed/alamode_crosscheck/` (committed; small band-resolved
    tau(300 K) validation tables against Tadano & Tsuneyuki 2015 Fig. 7)

This directory (`data/raw/alamode_sto/`) is gitignored except for this note. Only the tidy extracts in
`data/processed/` are committed.

## Revision calculations (v3.0.0): `revision/`

The calculations of the revised manuscript are under `revision/<variant>/<job id>/`, one directory per job of the
compute queue (`queue/run_queue.py`). The job manifest `data/provenance/manifest_final.json` records each job's
inputs, the sha256 of every file, its build id and its Zenodo archive (`archive_file`).

| variant directory | jobs | Zenodo archive |
|---|---|---|
| `pm_q12/` | production SCPH (q₁ = 12³, polarization mixing), RTA at every T and mesh, eigenvector runs | `scph_q12`, `pm_q12_rta`, `pm_q12_evec` |
| `pm_q16/` | SCPH on q₁ = 16³ and its four RTA runs | `scph_q16` |
| `pm_q8/`, `pm_q8_mix0p1/` | SCPH on q₁ = 8³, its RTA runs, and the reduced-mixing run | `pm_q8` |
| `pm_q12_4ph_audit/`, `pm_q12_4ph_audit_tol1e6/`, `pm_q12_4ph_trace/`, `pm_q12_4ph_v1omit/`, `pm_q12_4ph_v2strict/`, `pm_q12_4ph_v3dump/` | per-mode four-phonon runs of the rule, tolerance, trace and validation builds (`tools/alamode_patches/BUILD.md`) | `pm_q12_4ph_audit` |

To rerun `scripts/revision/`, unpack the archives of the Zenodo record (concept DOI 10.5281/zenodo.20559513) into
`data/raw/alamode_sto/revision/`. The runs made before the queue existed are in the archive `prequeue_runs`; point
`FEKAPPA_RUNS` at it.

## Force-constant provenance (state plainly wherever this data is used)

Quoted verbatim from the source `PROVENANCE.md` (kept alongside this note):

> These are ALAMODE's own tutorial force constants under ALAMODE's own
> DFT settings/functional — an external cross-check, not a first-principles
> result from this paper's own PBEsol pipeline.

Classification under this repository's provenance labels: `ext-FP-alamode`
(first-principles physics computed by us with published/tutorial force
constants; external FC provenance). Validation status: kappa(300 K) =
9.53 W/mK vs Tadano & Tsuneyuki 8.7 W/mK (~10%, on a coarser 8x8x8 vs 12x12x12
mesh); band-resolved tau(300 K) within the published Fig. 7 ranges for every
populated band; 0/525 unstable modes after SCPH coupling.

## SCPH supersedes bare at low frequency

The bare-RTA dataset `STO_RTA_production.result`/`.kl` (harmonic frequencies
from `STO_anharm.xml`, no SCPH renormalization) is SUPERSEDED by the
SCPH-coupled `STO_RTA_scph_*K.result` files for all paper-facing quantities:
its low-frequency manifold is quantitatively wrong (kappa(300 K) = 3.34 vs
9.53 W/mK; low-omega lifetimes 2-85x too short; 11/525 unstable modes). It is
retained solely as the controlled negative-control dataset (bare-RTA vs
SCPH-coupled comparison).

## Locally regenerated and extended runs (2026-08-18)

The ALAMODE 1.5.0 toolchain (anphon + dfc2) was rebuilt locally from
`github.com/ttadano/alamode` (default branch, commit c5b015c, 2026-06-17)
with the cmake recipe recorded in `PROVENANCE.md`, using the same example
inputs (`STO_anharm.xml`, `STO222.xml`, `STO_scph2-2.scph_dfc2`, `BORN`).
Validation: the regenerated 300 K SCPH-coupled run reproduces the imported
`STO_RTA_scph_300K.result` and `.kl` bit-for-bit.

New outputs added here under the same naming scheme (same tutorial force
constants -- the FC-provenance caveat above applies unchanged):

- `STO_RTA_scph_<T>K.{result,kl,kl_coherent,kc_elem}` for T = 0-1000 K in
  50 K steps (21 temperatures; the seven 100-400 K runs regenerated
  identically). `KAPPA_COHERENT = 2` added to the `&analysis` block:
  `.kl_coherent` is the coherent (Wigner off-diagonal) part of kappa,
  `.kc_elem` its mode-pair-resolved elements.
- `STO_RTA_scph_300K_mesh12.*`: single 12x12x12 SCPH-coupled convergence
  run at 300 K.
- `STO_band_scph_<T>K.bands`: renormalized-fc2 dispersion along
  Gamma -> (0.15, 0, 0) r.l.u. (61 points, NONANALYTIC = 3) per T, for the
  small-q Lambda^2(T) fit.

## Conventions pinned at source

- `#GAMMA_EACH` values are HWHM in cm-1 (convention verified by the Tadano
  cross-check recorded in `PROVENANCE.md`); tau = 1/(2*Gamma). Any unit
  conversion is applied exactly once, at the consumption layer.
- SCPH runs use `FCSXML = STO_anharm.xml` (anharmonic scattering terms) with
  `FC2XML = STO_renorm_<T>.xml` (dfc2-renormalized harmonic terms, one file
  per temperature); reference cell celldm(1) = 7.363 bohr (ALAMODE example
  cell, not this project's PBEsol cell).
