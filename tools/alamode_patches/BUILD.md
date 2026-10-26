# ALAMODE builds used in this work

Every `anphon` binary that ran a reported job, with the source patch, the binary sha256 and the build recipe. This
table matches SI Table S1 and `data/processed/si_v4/build_manifest.csv`. The per-job build ids are in
`data/provenance/manifest_final.json`, field `build`.

- **Base:** ALAMODE commit `c5b015caaf72c1a86a5250b6e34a86c0e09c7986`, unmodified, for the production build.
- **Patches:** each patch here applies to that commit on its own (`patch -p1` or `git apply` in a clean checkout). None
  is stacked on another.
  - Each patched tree differs from the base only in the file the patch touches: `anphon/selfenergy.cpp` (four-phonon
    builds) or `anphon/dynamical.cpp` (`absw`).
  - Applying each patch to a clean c5b015c checkout reproduces the source of the tree that was built, byte for byte.
    This was checked on 2026-10-24.
- **Toolchain:** conda-forge, osx-arm64, pinned in `env.explicit.txt` (sha256
  `7f762f10a185e08fa0a258db19af78c77cefca3b5b5abaaa4f0e45e58cbb89bf`). Create it with
  `conda create -p $TOOLCHAIN/env --file env.explicit.txt`.
- **Recipe:** `build_variant.sh`, the same CMake line for every tree. Run it as
  `TOOLCHAIN=<dir> sh build_variant.sh <tree> [<tree> …]`. The earlier per-tree scripts used for the production,
  audit and `absw` builds differ from it only in the tree name.
- **Binary hashes depend on the machine.** The absolute rpath of the toolchain is embedded, so the same source gives a
  different sha256 on HOME and on NODE. The build id of each job records the hash of the binary that ran it.
- **Machines:**
  - HOME: Apple M3 Max, 36 GB.
  - NODE: Apple M5 Pro, 48 GB.

## Builds that ran reported jobs

| variant | tree | patch (sha256 of the patch file) | machine | anphon sha256 | jobs |
|---|---|---|---|---|---|
| production | `alamode` | none | HOME | `23805d3b84c0ed625e04b52fbbe97930e47ecc55f61feac77861ea2e4e522569` | 103 |
| production | `alamode` | none | NODE | `266b507a02c2d42358326e9cb8f36c1dfee4f6fee3fe2fc0a7f7ce879e414881` | 27 |
| four-phonon rule: enlarged exclusion, identity tolerance 1e-9 Ry | `alamode_4ph_audit` | `4ph_audit.diff` (`653a330f4929cfea42a260a68dedbf5207132898724ce56fb647747519b6fa9d`) | NODE | `b6e7c732d04b9a1c50497e7cd4ff5a574170788406eda6eff1109cdfc09e5995` | 26 |
| tolerance sensitivity, 1e-6 Ry | `alamode_4ph_audit_tol1e6` | `4ph_audit_tol1e6.diff` (`a1d2ec7697b5cf70afc9cac4fec2ba89738a6c39f135aa71e6e901826e866e60`) | NODE | `66c4bd619882524aa5926840d2e11bb9aa8383adb8fc42afccbce22efc0024a0` | 2 |
| one-quartet trace (`TRACE.md`) | `alamode_4ph_trace` | `4ph_trace.diff` (`12bf571c7d22252f717592d3814d907759c2f9151575466201e5aa56a2014f15`) | NODE | `b4064feb18fb979b2058e172f940782bff131ca028ad92dc06e5d7a89409c544` | 1 |
| validation (1): omitted-sum accumulator | `alamode_4ph_v1omit` | `4ph_v1omit.diff` (`84e4dfae90afa66ee40684357dbec6b9a601a29a335266cba95dd8f691e66986`) | NODE | `ba9db2de378dc538acb01a632f7ec4ea2e1d940326e300159f59678abe3bc9b6` | 1 |
| validation (2): strict-index exclusion | `alamode_4ph_v2strict` | `4ph_v2strict.diff` (`4112045c1bf5ea35bdc60783cb2153837651d0b69e62bf9666b2841f641d8998`) | NODE | `9f899fc41967e2126b9a6c86889b87d2fdf5f10c7b23b8cb2f1c7c8a46230f9f` | 4 |
| validation (3): full-precision dump | `alamode_4ph_v3dump` | `4ph_v3dump.diff` (`3460354aae45210ec67db804643d40947863a0087430c2ecf69a6bc715a53fe9`) | NODE | `702793b11faa27eaff43d01a736633cc10c5580fde0bb166335933bfd88416f4` | 2 |

The job manifest also holds four entries without a binary. Three are imported runs from 2026-10-07 (the q1 = 8³ SCPH
run and two RTA runs on it), made with the HOME production toolchain before the queue existed. The fourth is the q₁ rule job
(`d0_decision`).

## Pre-queue build

| variant | tree | patch (sha256) | machine | anphon sha256 | used for |
|---|---|---|---|---|---|
| `|ω|` substitution for ω² < 0 | `alamode_absw` | `absw_omega_audit.diff` (`6d292815cd583ba538b1727a084a3fa25791123f08a305d7d8c0d145f9a9217f`) | HOME | `9244228d532bbccd46a8ac349be886cf27cdb4fc139d5f6dfc56ff88379b96d2` | the `|ω|` framing of the bare-harmonic diagnostic (SI Table S32); outputs in the Zenodo variant `prequeue_runs` |

A HOME build of `alamode_4ph_audit` exists, with sha256 `2807ac8d952847fb4f86827e07f1ca22e946feb35dca5d4f60923c41adf680c5`.
It ran no reported job.

## Other files here

- `apply_variants.py`: writes the trace build and the three validation builds into clean c5b015c copies. Their four
  `*.diff` files are `git diff` of the resulting trees. The `4ph_audit`, `4ph_audit_tol1e6` and `absw` patches were
  written by hand.
- `trace_one_quartet.py` and `TRACE.md`: the hand check of one printed quartet, covering sign, prefactor and units.
- `diff_excluded_sets.py`: compares the excluded quartets dumped by validation build (3) with the independent count
  from the eigenvalue file. It imports `fekappa.evec_io`.
