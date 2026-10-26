"""Closure of the production tables (data/processed/sto_production/, own SCPH q1 = 12^3 with polarization mixing)."""
from pathlib import Path

import numpy as np
import pandas as pd

PROD = Path(__file__).resolve().parents[1] / "data/processed/sto_production"
rd = lambda f: pd.read_csv(PROD / f, comment="#")


def test_total_is_particle_plus_coherence():
    k = rd("kappa_mesh_T.csv")
    assert len(k) == 57
    np.testing.assert_allclose(k.kappa_total, k.kappa_P + k.kappa_C, rtol=1e-4)
    assert (k.closure_rel < 1e-3).all()  # pair-sum closure recorded at validation


def test_sector_shares_sum_to_mode_table_and_solver():
    s = rd("kappa_sectors_T.csv")
    np.testing.assert_allclose(s.kappa_P_A + s.kappa_P_soft + s.kappa_P_hard, s.kappa_P, rtol=1e-5)
    m = pd.read_csv(PROD / "mode_table_T.csv.gz")
    msum = m.groupby(["T_K", "mesh"]).kappa_mode_W_mK.sum()
    for r in s.itertuples():
        assert abs(msum[(r.T_K, r.mesh)] / r.kappa_P - 1) < 1e-5
    k = rd("kappa_mesh_T.csv").set_index(["T_K", "mesh"]).kappa_P
    e = rd("eta_validity_T.csv").set_index(["T_K", "mesh"]).closure_rel
    for r in s.itertuples():  # mode-table kappa_P against the solver value, within the recorded closure
        assert abs(r.kappa_P / k[(r.T_K, r.mesh)] - 1) <= e[(r.T_K, r.mesh)] + 1e-5


def test_coherence_blocks_sum_to_pair_sum():
    c = rd("coherence_blocks_T.csv")
    blocks = [x for x in c.columns if x.startswith("kC_")]
    assert len(blocks) == 6
    np.testing.assert_allclose(c[blocks].sum(axis=1), c.kappa_C_pairsum, rtol=1e-4)
    k = rd("kappa_mesh_T.csv").set_index(["T_K", "mesh"]).kappa_C
    for r in c.itertuples():
        assert abs(r.kappa_C_pairsum / k[(r.T_K, r.mesh)] - 1) < 1e-3


def test_character_projections_sum_to_one():
    m = pd.read_csv(PROD / "mode_table_T.csv.gz", usecols=["c_A", "c_soft", "c_hard", "sector_char"])
    m = m[m.sector_char.notna()]
    np.testing.assert_allclose(m.c_A + m.c_soft + m.c_hard, 1.0, atol=1e-5)
