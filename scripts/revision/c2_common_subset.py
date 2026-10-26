"""C2: bare-harmonic vs SCPH-coupled kappa_P on the full and the common stable mode set.

Bare = data/raw/alamode_sto/STO_RTA_production.result (21 T, bare fc2).
SCPH = data/raw/alamode_sto/STO_RTA_scph_<T>K.result (dfc2-renormalized fc2 per T).
Both 8x8x8, identical irreducible q list (35 q x 15 branches = 525 modes).

The unstable set U = (q_index, branch) with omega < 0 on the bare surface
(T-independent: bare frequencies do not depend on T). anphon 1.5.0 assigns
those modes lifetime 0 (conductivity.cpp:420-424) and removes them as
scattering partners (anharmonic_core.cpp:330-331, V3 = 0), so in the bare run
they contribute exactly zero to kappa_P. The common stable subset is the
complement of U; we evaluate the SCPH kappa_P restricted to it.

Output: data/processed/revision/c2_common_subset.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd

from fekappa.alamode_io import kappa_isotropic, read_kl
from fekappa.kappa_decomposition import mode_kappa_table

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/raw/alamode_sto"
OUT = ROOT / "data/processed/revision/c2_common_subset.csv"

bare = mode_kappa_table(RAW / "STO_RTA_production.result")
U = bare.loc[(bare["T_K"] == bare["T_K"].min()) & (bare["omega_cm1"] < 0), ["q_index", "branch", "omega_cm1", "multiplicity"]]
assert len(U) == 11, len(U)
ukeys = set(zip(U["q_index"], U["branch"]))
n_star_U = int(U["multiplicity"].sum())

bare_kl = kappa_isotropic(read_kl(RAW / "STO_RTA_production.kl"))
rows = []
for t in range(50, 1001, 50):
    b = bare[bare["T_K"] == float(t)]
    s = mode_kappa_table(RAW / f"STO_RTA_scph_{t}K.result")
    in_u_b = np.array([(q, br) in ukeys for q, br in zip(b["q_index"], b["branch"])])
    in_u_s = np.array([(q, br) in ukeys for q, br in zip(s["q_index"], s["branch"])])
    kb_full = b["kappa_mode_W_mK"].sum()
    kb_sub = b.loc[~in_u_b, "kappa_mode_W_mK"].sum()
    ks_full = s["kappa_mode_W_mK"].sum()
    ks_sub = s.loc[~in_u_s, "kappa_mode_W_mK"].sum()
    ks_U = s.loc[in_u_s, "kappa_mode_W_mK"].sum()
    rows.append(
        {
            "T_K": t,
            "bare_kappaP_solver": bare_kl.loc[bare_kl["T_K"] == float(t), "kappa_W_mK"].item(),
            "bare_kappaP_modesum_full": kb_full,
            "bare_kappaP_modesum_subset": kb_sub,
            "scph_kappaP_solver": kappa_isotropic(read_kl(RAW / f"STO_RTA_scph_{t}K.kl"))["kappa_W_mK"].item(),
            "scph_kappaP_modesum_full": ks_full,
            "scph_kappaP_modesum_subset": ks_sub,
            "scph_kappaP_on_U": ks_U,
            "scph_U_share_pct": 100 * ks_U / ks_full,
            "ratio_full": ks_full / kb_full,
            "ratio_subset": ks_sub / kb_sub,
            "ratio_shift_pct": 100 * (ks_sub / kb_sub - ks_full / kb_full) / (ks_full / kb_full),
        }
    )

df = pd.DataFrame(rows)
with open(OUT, "w") as fh:
    fh.write(
        "# C2 common-stable-subset audit; U = 11 bare-unstable irreducible modes "
        f"({n_star_U} of 7680 reducible mesh modes = 512 q x 15); subset = 514 irreducible modes\n"
        "# U list (q_index, branch, bare omega cm^-1): "
        + "; ".join(f"({q},{b},{w:.2f})" for q, b, w in U[["q_index", "branch", "omega_cm1"]].itertuples(index=False))
        + "\n"
    )
    df.to_csv(fh, index=False, float_format="%.6g")
pd.set_option("display.width", 250)
print(U.to_string(index=False))
print(df.round(4).to_string(index=False))


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"c2_common_subset.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['data/raw/alamode_sto (distributed example outputs; IMPORT_NOTE.md)'])
