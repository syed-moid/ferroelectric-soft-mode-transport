"""fekappa.evec_io on a stored 4^3 eigenvector file (SCPH q1 = 12^3, polarization mixing, 300 K)."""
from pathlib import Path

import numpy as np

from fekappa.evec_io import NS, kmap, read_evec

FIXTURE = Path(__file__).parent / "data" / "evec_pm_q12_T300_m4.mesh.evec.gz"


def test_shapes():
    xk, om, ev = read_evec(FIXTURE)
    assert xk.shape == (64, 3)
    assert om.shape == (64, NS)
    assert ev.shape == (64, NS, NS)


def test_gamma_translations_and_stability():
    xk, om, ev = read_evec(FIXTURE)
    g = int(np.argmin(np.abs(xk).sum(axis=1)))
    assert np.allclose(xk[g], 0.0)
    assert np.all(np.abs(om[g, :3]) < 1e-20)  # the three translations
    assert np.all(om > -1e-20)  # stable renormalized spectrum: no negative eigenvalue on the mesh


def test_orthonormal_eigenvectors():
    _, _, ev = read_evec(FIXTURE)
    for k in range(ev.shape[0]):
        assert np.abs(ev[k].conj() @ ev[k].T - np.eye(NS)).max() < 1e-6


def test_kmap_finds_mesh_points():
    xk, _, _ = read_evec(FIXTURE)
    q = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [-0.25, 0.0, 0.25]])
    km = kmap(xk, q)
    d = ((xk[km] - q) + 0.5) % 1.0 - 0.5
    assert np.abs(d).max() < 1e-6
