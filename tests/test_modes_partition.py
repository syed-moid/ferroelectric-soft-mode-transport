"""Tests of fekappa.modes and fekappa.partition (helpers of the v4 analysis scripts)."""
import numpy as np
import pandas as pd
import pytest

from fekappa.modes import TOL_DEG_CM1, cell_volume_ang3, deg_groups, group_average_gamma, mesh_of
from fekappa.partition import RY_TO_CM, read_gamma, read_t3


def test_tolerance_is_anphon_tol_omega_in_cm1():
    assert TOL_DEG_CM1 == pytest.approx(1.0e-7 * 109737.31568160)


def test_deg_groups_joins_to_first_mode_of_group():
    tol = TOL_DEG_CM1
    om = [10.0, 10.0 + 0.5 * tol, 10.0 + 0.9 * tol, 10.0 + 1.1 * tol, 20.0]
    # joined while |w - w_first| < tol, measured from the first mode of the group, not the previous one
    assert deg_groups(om).tolist() == [0, 0, 0, 3, 4]


def test_group_average_gamma_averages_within_groups_only():
    df = pd.DataFrame({"q_index": [1, 1, 1, 2, 2], "branch": [1, 2, 3, 1, 2],
                       "omega_cm1": [5.0, 5.0, 9.0, 5.0, 6.0], "gamma_hwhm_cm1": [1.0, 3.0, 4.0, 2.0, 8.0]})
    g = group_average_gamma(df)
    assert g.tolist() == [2.0, 2.0, 4.0, 2.0, 8.0]
    assert df.gamma_hwhm_cm1.tolist() == [1.0, 3.0, 4.0, 2.0, 8.0]   # input unchanged


RTA_IN = """&general
 PREFIX = x
/
&cell
7.363
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
/
&kpoint
2
16 16 16
/
"""


def test_mesh_and_cell_volume_from_rta_in(tmp_path):
    (tmp_path / "rta.in").write_text(RTA_IN)
    assert mesh_of(tmp_path) == 16
    a_ang = 7.363 * 0.529177210903
    assert cell_volume_ang3(tmp_path) == pytest.approx(a_ang ** 3, rel=1e-6)


def test_read_t3_totals_and_bins(tmp_path):
    p = tmp_path / "t3.txt"
    p.write_text("T 1 4 x 300.0 x 2.0e-5\nB 1 4 300.0 0 2 1.0e-5 0.5e-5\nB 1 4 300.0 2 2 0.2e-5 0.3e-5\n")
    tot, bins = read_t3(p)
    assert tot == {(1, 4, 300.0): 2.0e-5}
    assert bins[(1, 4, 300.0)] == [(0, 2, 1.0e-5, 0.5e-5), (2, 2, 0.2e-5, 0.3e-5)]
    assert sum(d + a for _, _, d, a in bins[(1, 4, 300.0)]) == pytest.approx(tot[(1, 4, 300.0)])


def test_read_t3_rejects_duplicate_carrier(tmp_path):
    p = tmp_path / "t3.txt"
    p.write_text("T 1 4 x 300.0 x 2.0e-5\nT 1 4 x 300.0 x 2.0e-5\n")
    with pytest.raises(AssertionError):
        read_t3(p)


def test_read_gamma_from_result_blocks(tmp_path):
    p = tmp_path / "x.result"
    p.write_text("#TEMPERATURE\n300 400 100\n#END TEMPERATURE\n"
                 "#GAMMA_EACH\n2 3\n1\n0 0 0\n0.5\n0.7\n#END GAMMA_EACH\n")
    g = read_gamma(p)
    assert g == {(1, 2, 300.0): 0.5, (1, 2, 400.0): 0.7}


def test_ry_to_cm_matches_rydberg_in_wavenumbers():
    assert RY_TO_CM == pytest.approx(109737.3, rel=1e-6)
    assert np.isfinite(RY_TO_CM)
