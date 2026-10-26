"""Readers for anphon outputs: the `.mesh.evec` eigenvector file (PRINTEVEC = 1, KPMODE = 2) and the irreducible
wave-vector list of a `.result` file. Used by the revision analysis scripts (scripts/revision/)."""
import gzip
from pathlib import Path

import numpy as np

NS = 15


def read_evec(path):
    """(wave vectors [nk, 3], eigenvalues [nk, NS], eigenvectors [nk, mode, component]) of a `.mesh.evec` file
    (plain or gzipped). The eigenvalue column is the one anphon prints after "### mode"."""
    xk, om, ev, cur = [], [], [], None
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("## kpoint"):
                xk.append([float(v) for v in line.split(":")[1].split()])
            elif line.startswith("### mode"):
                om.append(float(line.split(":")[1])); cur = []; ev.append(cur)
            elif cur is not None and line.strip() and not line.startswith("#"):
                a, b = line.split(); cur.append(complex(float(a), float(b)))
    xk = np.array(xk)
    return xk, np.array(om).reshape(len(xk), NS), np.array(ev).reshape(len(xk), NS, NS)


def irred_q(result_path):
    lines = Path(result_path).read_text().splitlines()
    i = lines.index("#KPOINT"); n = int(lines[i + 2])
    return np.array([[float(v) for v in l.split(":")[1].split()[:3]] for l in lines[i + 3:i + 3 + n]])


def kmap(xk_mesh, q_irr):
    out = []
    for q in q_irr:
        d = np.abs(((xk_mesh - q) + 0.5) % 1.0 - 0.5).sum(axis=1)
        j = int(np.argmin(d)); assert d[j] < 1e-6; out.append(j)
    return np.array(out)


def irred_evecs(evec_path, result_path):
    """Eigenvectors at the irreducible q of `result_path`, shape [n_irred, mode, component]."""
    xk, om, ev = read_evec(evec_path)
    km = kmap(xk, irred_q(result_path))
    return ev[km], om[km]
