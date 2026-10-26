# Compute queue

`run_queue.py` runs the ALAMODE (anphon) jobs of this work from job specs. For every job it writes the anphon input
from the spec, runs it in an emptied work directory, and validates the outputs. It then copies them to
`data/raw/alamode_sto/revision/<variant>/<job id>/` and records parameters, build id, timings and the sha256 of every
output in the job manifest. The scrubbed final manifest is `data/provenance/manifest_final.json`; its job ids and
build ids are the ones cited in SI Table S1. The runner uses the Python standard library only (3.9+).

## Commands

```
QUEUE_TOOLCHAIN=<dir> python3 queue/run_queue.py add queue/spec_<name>.json   # add the jobs of a spec
QUEUE_TOOLCHAIN=<dir> python3 queue/run_queue.py run [--until "YYYY-MM-DDTHH:MM:SS <IANA zone>"]
QUEUE_TOOLCHAIN=<dir> python3 queue/run_queue.py status                      # write and print queue/status.md
QUEUE_TOOLCHAIN=<dir> python3 queue/run_queue.py requeue <job id>            # failed -> pending
```

- `run` works through pending jobs in priority order once their dependencies are done.
- `--until` starts no job that would not finish before that time, and stops a running one at that time.
- `queue/STOP` stops after the current job; `queue/STOP_NOW` stops at once and re-queues the running job.

## Job specs

`spec_*.json` are the job definitions of every calculation in the manifest. Each is a JSON list of jobs with:
- `id`, `kind` (`scph`, `rta`, `evec`, `mode4ph`), `priority`, `deps` (job ids) and `est_s` (expected seconds);
- `params`:
  - `q1` (SCPH grid) and `pm` (polarization mixing);
  - `T` or `tmin`/`tmax`/`dt`, `mesh`, and `variant` (output directory);
  - for `mode4ph`: `epsilon` (cm⁻¹), `modes` (wave vector and branch, 1-based), `mpi_np`, and `anphon_tree` (a
    patched build in `tools/alamode_patches/`);
  - optional: `mixalpha`, `measure_mem`, `rss_kill_gb` (see the docstring of `run_queue.py`).

## Environment

| variable | meaning |
|---|---|
| `QUEUE_TOOLCHAIN` (required) | directory holding the ALAMODE builds (`alamode*/_build/`), the conda environment `env/` (`tools/alamode_patches/env.explicit.txt`), the example inputs `sto_example/` and the work area `runs/` |
| `QUEUE_DIR` | location of the manifest and runtime files (default: this directory) |
| `QUEUE_OMP_THREADS` | OpenMP threads per job (default 12) |

## Tests

`test_queue.py` covers validation, manifest atomicity, recovery, retries, STOP handling and the input writers. It
copies stored runs from `$QUEUE_TOOLCHAIN/runs/` before corrupting them. Run it as
`QUEUE_TOOLCHAIN=<dir> python3 queue/test_queue.py`.
