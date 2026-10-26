"""C1: collect the 300 K benchmark matrix (PM x SCPH q1 mesh x transport mesh).

off_k2 rows = the submitted dataset (distributed tutorial dfc2; data/raw/alamode_sto, Aug build).
All other rows = our own MODE=SCPH runs ($FEKAPPA_RUNS/scph_<cfg>) -> dfc2 at 300 K
-> anphon RTA with the submitted RTA input (tetrahedron, quantum, no isotope, KAPPA_COHERENT=2).
Reference: Tadano & Tsuneyuki 2015 Fig. 7, 8.7 W/mK at 300 K (12^3 BTE, 8^3 SCPH q1 grid).
Output: data/processed/revision/c1_benchmark_matrix.csv
"""

from pathlib import Path

from fekappa.provenance import external_dir
import pandas as pd

from fekappa.alamode_io import kappa_isotropic, read_kl, read_kl_coherent

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/raw/alamode_sto"
RUNS = external_dir("FEKAPPA_RUNS", "the pre-queue run directory (Zenodo prequeue_runs)")
REF = 8.7
soft = pd.read_csv(ROOT / "data/processed/revision/c1_scph_softmodes.csv", comment="#")
soft = soft[soft.T_K == 300].pivot(index="config", columns="point", values="omega_meV")


def kap(prefix_path):
    p = kappa_isotropic(read_kl(Path(str(prefix_path) + ".kl")))["kappa_W_mK"].item()
    c = read_kl_coherent(Path(str(prefix_path) + ".kl_coherent"))[["k_xx", "k_yy", "k_zz"]].mean(axis=1).item()
    return p, c


rows = []
for cfg in ["off_k2", "on_k2", "off_k4", "on_k4", "off_k8", "on_k8"]:
    for m in (8, 12):
        if cfg == "off_k2":
            pref = RAW / ("STO_RTA_scph_300K" if m == 8 else "STO_RTA_scph_300K_mesh12")
            src = "submitted (distributed dfc2)"
        else:
            pref = RUNS / f"rta_{cfg}_m{m}" / f"rta_{cfg}_m{m}"
            src = "own MODE=SCPH"
        if not Path(str(pref) + ".kl").exists():
            continue
        p, c = kap(pref)
        rows.append({"config": cfg, "PM": cfg.startswith("on"), "kmesh_scph": int(cfg[-1]), "transport_mesh": m,
                     "source": src, "kappa_P": p, "kappa_C": c, "kappa_total": p + c,
                     "kappa_P_vs_tadano_pct": 100 * (p - REF) / REF,
                     "Omega_Gamma_TO1_meV": soft.loc[cfg, "Gamma"], "Omega_M_meV": soft.loc[cfg, "M"],
                     "Omega_R_meV": soft.loc[cfg, "R"]})
df = pd.DataFrame(rows)
with open(ROOT / "data/processed/revision/c1_benchmark_matrix.csv", "w") as fh:
    fh.write("# C1 300 K benchmark matrix; PM = SELF_OFFDIAG=1; kmesh_scph = KMESH_SCPH (KMESH_INTERPOLATE 2x2x2 fixed); "
             "reference kappa_P = 8.7 W/mK (Tadano & Tsuneyuki 2015, 12^3 BTE, 8^3 SCPH, PM implied)\n")
    df.to_csv(fh, index=False, float_format="%.5g")
pd.set_option("display.width", 220)
print(df.round(4).to_string(index=False))


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"c1_benchmark_matrix.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['data/processed/revision/c1_scph_softmodes.csv', 'data/raw/alamode_sto (distributed example outputs; IMPORT_NOTE.md)', '$FEKAPPA_RUNS/rta_{on,off}_k*_m* (zenodo:prequeue_runs.tar.gz)'])
