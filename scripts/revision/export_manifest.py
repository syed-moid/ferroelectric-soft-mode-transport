"""Write the scrubbed final job manifest, data/provenance/manifest_final.json, from the live queue manifest.

The live manifest (queue/manifest.json, not tracked) records host names and absolute paths. The export applies the
scrub of the Zenodo package (tools/zenodo/make_zenodo_package.py): host names become HOME / NODE, the host suffix of
each build id is dropped (the anphon sha256 stays), account names and home-directory paths become <user> / <home>,
and internal notes are dropped from `meta`. Input and output sha256 values are kept as recorded: no archived output
file contains a host name or a home path, so the scrub never changes a file hash. Each job gains `archive_file`,
the Zenodo archive that holds its outputs.

    python scripts/revision/export_manifest.py [--toolchain DIR]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools/zenodo"))
import make_zenodo_package as zp  # noqa: E402

OUT = ROOT / "data/provenance/manifest_final.json"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--toolchain", default=str(ROOT.parent / "_toolchain"))
    a = ap.parse_args()
    manifest = json.loads((ROOT / "queue/manifest.json").read_text())
    rules, info = zp.build_scrub(manifest, Path(a.toolchain))
    archive_of = {}
    for jid, j in manifest["jobs"].items():
        v = zp.variant_of(j)
        if v:
            archive_of[jid] = f"{v}.tar.gz"
    out = zp.scrub_manifest(manifest, rules, archive_of)
    txt = json.dumps(out, indent=1, sort_keys=True) + "\n"
    hits = [p.pattern for p, _ in rules if p.search(txt)]
    if hits:
        raise SystemExit(f"scrub patterns still present: {hits}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(txt)
    n = len(out["jobs"])
    print(f"wrote {OUT.relative_to(ROOT)}: {n} jobs, {len(archive_of)} with an archive; "
          f"hosts {sorted(set(info['hosts'].values()))}")


if __name__ == "__main__":
    main()
