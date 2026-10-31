# Provenance record of the revision calculations

`manifest_final.json` is the final job manifest of the compute queue (`queue/run_queue.py`). It lists 170 entries, all
`done`: 169 lattice-dynamics runs and one `decision` entry (`d0_decision`) that applied the SCPH-grid selection
rule (R1 of `docs/recorded_rules.md`), with no solver run. It is written by `scripts/revision/export_manifest.py` from the live, untracked `queue/manifest.json`. The
raw outputs of every job are in the Zenodo record; the archive is named in each job's `archive_file`.

## Fields of a job record

| field | meaning |
|---|---|
| `id` | job name; `kind` + `params` define the run |
| `kind` | `scph` (SCPH temperature sweep), `rta` (Wigner-RTA transport at one T and mesh), `evec` (eigenvectors on a mesh), `mode4ph` (per-mode four-phonon self-energy, diagram (c)), `decision` (the q₁ grid rule) |
| `params` | `q1` (SCPH grid), `pm` (polarization mixing), `T`, `mesh` (transport or four-phonon mesh), `variant` (output directory), `extra` (other input keywords) |
| `deps` | jobs whose outputs this job reads (the SCPH solution for every transport run) |
| `status`, `attempts`, `history` | the run record: start, done, killed or failed events with UTC time and machine |
| `host` | machine that produced the archived outputs: `HOME` or `NODE` |
| `build` | `alamode@c5b015ca/anphon-sha256:<12 hex>` for the production build, `AUDIT:<tree>/anphon-sha256:<12 hex>` for a patched build (full hashes in `tools/alamode_patches/BUILD.md`), or `imported:` for three runs made before the queue existed |
| `start`, `end`, `duration_s`, `est_s` | wall-clock record |
| `outputs` | sha256 of every archived file, recorded when the job passed validation |
| `metrics` | values parsed from the outputs at validation (κ<sub>P</sub>, κ<sub>C</sub>, pair-sum closure, SCPH convergence, Γ₄ per mode, the `AUDIT` counters) |
| `archive` | output directory in the working repository (`data/raw/alamode_sto/revision/<variant>/<id>`) |
| `archive_file` | Zenodo archive with that directory |
| `excluded` | set for the two timing probes; they are not part of any reported result |
| `import` | for the imported runs: source directory, with home paths replaced by `<home>` |

## Scrub

- Host names are replaced by `HOME` and `NODE`.
- The host suffix of every build id is dropped; the anphon sha256 is kept.
- Local paths are replaced by `<runs>` (run directories), `<toolchain>`, `<repo>` or `<home>`, and account names by
  `<user>`.
- No archived output file contains a host name or a home path. The scrub therefore changes no recorded file hash.

## Machines and toolchain

| | HOME | NODE |
|---|---|---|
| hardware | Apple M3 Max, 36 GB | Apple M5 Pro, 48 GB |
| jobs (archived outputs) | 107 | 63 |

- ALAMODE: commit `c5b015caaf72c1a86a5250b6e34a86c0e09c7986`.
- Patches, binary hashes and the build recipe: `tools/alamode_patches/BUILD.md`.
- Toolchain: conda-forge, pinned in `tools/alamode_patches/env.explicit.txt` (sha256
  `7f762f10a185e08fa0a258db19af78c77cefca3b5b5abaaa4f0e45e58cbb89bf`).
- Cross-machine check (2026-10-12): NODE reran the HOME job `rta_pm_q12_T300_m12` in an isolated queue; all outputs
  were bit-identical.

## Queue runner version

The jobs from spec A1 on ran with `queue/run_queue.py` of sha256
`6ae11b6d0e67e4d586dfd6dc73009d954462d7b7004a39b503fc0cffba65fcea`. The tracked copy differs only in comments and
messages, in requiring `QUEUE_TOOLCHAIN`, and in naming the end-time option `--until`. The job inputs it writes are
unchanged, and `queue/test_queue.py` passes.
