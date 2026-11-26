#!/usr/bin/env python3
"""Tables of the paper (v4) as CSV, with provenance headers, from the regenerated outputs.

Every pipe table of the Markdown report tables (PROC/{a,b}/reports/*.md, written by regenerate_v4.py; HELD/reports
for the part that needs the files of the 2015 study) becomes one CSV: PROC/tables/{a,b}/<report>__<n>.csv (HELD/tables/ for that part). The
first lines of each CSV ('#') name the producing script, the report file, the heading above the table and the
part. The index PROC/tables/INDEX.csv lists every CSV with its heading and size.

Parts: a = model A only; b = model B, or A and B side by side; held = needs the files of the 2015 study; not in this release.
Usage (repository root): build_tables_v4.py PROC [--held HELD]
"""
import csv
import re
import sys
from pathlib import Path

SEP = re.compile(r"^\s*\|[-:| ]+\|\s*$")


def md_tables(path):
    out, head, h2, lines, i = [], "", "", Path(path).read_text().splitlines(), 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("#"):
            head = ln.lstrip("#").strip()
            if ln.startswith("## "):
                h2 = head
            elif ln.startswith("### ") and h2:
                head = f"{h2} — {head}"
        if ln.lstrip().startswith("|") and i + 1 < len(lines) and SEP.match(lines[i + 1]):
            cells = lambda s: [c.strip().replace("**", "").replace("`", "") for c in s.strip().strip("|").split("|")]  # noqa: E731
            header, rows, i = cells(ln), [], i + 2
            while i < len(lines) and lines[i].lstrip().startswith("|"):
                rows.append(cells(lines[i]))
                i += 1
            out.append({"heading": head, "header": header, "rows": rows})
            continue
        i += 1
    return out


def part_of(report, model_dir):
    name = report.stem
    if model_dir == "held":
        return "held"
    if name.endswith("_A") or model_dir == "a":
        return "a"
    return "b"


def main():
    argv = sys.argv[1:]
    held = None
    if "--held" in argv:
        i = argv.index("--held")
        held = Path(argv[i + 1])
        del argv[i:i + 2]
    proc = Path(argv[0])
    index = []
    srcs = [(proc / m / "reports", m) for m in ("a", "b")] + ([(held / "reports", "held")] if held else [])
    for rdir, m in srcs:
        if not rdir.exists():
            continue
        for rep in sorted(rdir.glob("*.md")):
            prt = part_of(rep, m)
            dest = (held / "tables" if prt == "held" else proc / "tables" / prt)
            dest.mkdir(parents=True, exist_ok=True)
            for n, t in enumerate(md_tables(rep), 1):
                f = dest / f"{rep.stem}__{n:02d}.csv"
                with open(f, "w", newline="") as fh:
                    fh.write(f"# script: scripts/v4/build_tables_v4.py from {rep.parent.parent.name}/reports/{rep.name}"
                             f" (table {n}); part {prt}\n# heading: {t['heading']}\n")
                    w = csv.writer(fh)
                    w.writerow(t["header"])
                    w.writerows(t["rows"])
                index.append({"file": str(f.relative_to(held if prt == "held" else proc)), "part": prt,
                              "report": rep.name, "table": n, "heading": t["heading"], "rows": len(t["rows"])})
    for prt, root in (("a/b", proc), ("held", held)):
        if root is None:
            continue
        rows = [r for r in index if (r["part"] == "held") == (prt == "held")]
        if not rows:
            continue
        with open(root / "tables" / "INDEX.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(prt, len(rows), "tables ->", root / "tables")


if __name__ == "__main__":
    main()
