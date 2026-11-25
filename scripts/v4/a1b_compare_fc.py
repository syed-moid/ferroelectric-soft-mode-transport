#!/usr/bin/env python3
"""Compare two ALAMODE force-constant XML files entry by entry.

Entries are matched by physical index: the order tag (FC2 ... FC6) and the
pair attributes (atom, Cartesian component and cell of every index).
Prints one CSV row per order to stdout.

With .fcs files, the irreducible tables are compared instead (8 printed
digits; needed for FC5 and FC6, which the XML does not hold).

Usage: a1b_compare_fc.py REFERENCE.{xml,fcs} CANDIDATE.{xml,fcs} LABEL
"""
import csv
import hashlib
import math
import sys
import xml.etree.ElementTree as ET

ORDERS = ("FC2", "FC3", "FC4", "FC5", "FC6")


def load(path):
    root = ET.parse(path).getroot()
    fcs = {o: {} for o in ORDERS}
    for el in root.iter():
        if el.tag in fcs:
            key = tuple(el.attrib[k] for k in sorted(el.attrib))
            if key in fcs[el.tag]:
                raise SystemExit(f"duplicate entry {el.tag} {key} in {path}")
            fcs[el.tag][key] = float(el.text)
    return fcs


def load_fcs(path):
    """Irreducible '*FCn' tables of an ALM .fcs file (8 printed digits), keyed by order and pair strings."""
    fcs, cur = {}, None
    with open(path) as f:
        for line in f:
            t = line.split()
            if not t:
                continue
            if t[0].startswith("**FC"):
                cur = None
            elif t[0].startswith("*FC"):
                cur = "FC" + t[0][3:]
                fcs[cur] = {}
            elif cur and t[0].isdigit() and len(t) >= 6:
                key = tuple(t[4:-1])
                fcs[cur][key] = float(t[2])
    return fcs


def sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def main():
    ref_path, cand_path, label = sys.argv[1:4]
    if ref_path.endswith(".fcs"):
        ref, cand = load_fcs(ref_path), load_fcs(cand_path)
    else:
        ref, cand = load(ref_path), load(cand_path)
    w = csv.writer(sys.stdout, lineterminator="\n")
    for o in ORDERS:
        if o not in ref and o not in cand:
            continue
        r, c = ref.get(o, {}), cand.get(o, {})
        keys = set(r) | set(c)
        only_r = sum(1 for k in keys if k not in c)
        only_c = sum(1 for k in keys if k not in r)
        diffs = [abs(r.get(k, 0.0) - c.get(k, 0.0)) for k in keys]
        if not keys:
            continue
        maxref = max((abs(v) for v in r.values()), default=0.0)
        maxdiff = max(diffs, default=0.0)
        rms = math.sqrt(sum(d * d for d in diffs) / len(diffs)) if diffs else 0.0
        rel = maxdiff / maxref if maxref else float("nan")
        w.writerow([label, o, len(r), len(c), only_r, only_c, f"{maxref:.6e}",
                    f"{maxdiff:.6e}", f"{rms:.6e}", f"{rel:.3e}",
                    "PASS" if (only_r == 0 and only_c == 0 and rel <= 1e-6) else "FAIL",
                    sha256(ref_path)[:16], sha256(cand_path)[:16]])


if __name__ == "__main__":
    main()
