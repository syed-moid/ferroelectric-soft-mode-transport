# Provenance of the v4 runs

- `run_map.csv`: every run of v4 (parts a and b) by what it does (`run`: model, kind, grid, mesh, temperature,
  setting) next to its internal job name (`job`), with the build and the Zenodo archive that holds it. Runs reused
  from v3 point to the v3 manifest (`data/provenance/manifest_final.json`) for their build ids.
- `manifest_v4.json`: per run, every file of the run directory with its size and sha256 (symbolic links recorded as
  links, without target), and the build id: `alamode@c5b015ca/anphon-sha256:266b507a02c2` (base build on the compute
  node) or `PATCH:t3_partition/anphon-sha256:a0a5d30ca294` (linewidth-partition build; full hashes in
  `tools/alamode_patches/BUILD.md`).
- Scrub: no host name, user name or local path appears in these two files. `sha256` is the hash of each file as
  written by the solver. 21 `renorm.xml` files record a home path of the compute node in their header; in the Zenodo
  archives that path is replaced by `<home>`, and for these files the manifest also gives `sha256_archived` and
  `bytes_archived`, the hash and size of the file as archived.
- The runs that need the files of the 2015 study are not listed here; those files are not in this release
  (`README.md`, "Data").
