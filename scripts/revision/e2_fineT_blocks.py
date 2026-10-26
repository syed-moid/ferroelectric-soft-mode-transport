"""E2b: coherence pair-sector blocks on a 10 K grid, 250-350 K (8^3, tutorial SCPH settings).

250/300/350 K = submitted files (data/raw/alamode_sto); 260-340 K = own SCPH (TMIN 250 TMAX 350 DT 10,
tutorial settings; reproduces the distributed dfc2 at 250/300/350 K) -> dfc2 -> RTA, same input.
Also records the largest single pair contribution and its frequency splitting.
Output: data/processed/revision/e2_fineT_blocks.csv
"""
from pathlib import Path
from fekappa.provenance import external_dir
import pandas as pd
from fekappa.kappa_decomposition import coherence_blocks_row, read_kc_elem
from fekappa.alamode_io import read_kl_coherent
ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/raw/alamode_sto"; RUNS = external_dir("FEKAPPA_RUNS", "the pre-queue run directory (Zenodo prequeue_runs)")
rows = []
for t in range(250, 351, 10):
    p = RAW / f"STO_RTA_scph_{t}K" if t % 50 == 0 else RUNS / f"fineT_m8_{t}K" / f"fine_{t}K"
    row = {"T_K": t, **coherence_blocks_row(Path(str(p) + ".kc_elem"))}
    kc = read_kc_elem(Path(str(p) + ".kc_elem")); kc = kc[(kc.i == kc.j) & (kc.ibranch != kc.jbranch)]
    g = kc.groupby(["ik", "ibranch", "jbranch"]).agg(v=("kc_real", "sum"), w1=("omega1_cm1", "first"), w2=("omega2_cm1", "first"))
    g["v"] = g.v / 3 / 512 * 2  # unordered pair = (i,j)+(j,i)
    top = g.sort_values("v").iloc[-1]
    row.update({"top_pair_W_mK": top.v, "top_pair_ik": top.name[0], "top_pair_dw_cm1": abs(top.w1 - top.w2)})
    rows.append(row)
df = pd.DataFrame(rows)
df["soft_involved"] = df.kappa_C_A_soft + df.kappa_C_hard_soft + df.kappa_C_soft_soft
df.to_csv(ROOT / "data/processed/revision/e2_fineT_blocks.csv", index=False, float_format="%.5g")
pd.set_option("display.width", 250)
print(df.round(4).to_string(index=False))


# provenance header (fekappa.provenance)
from fekappa.provenance import stamp  # noqa: E402
for _f in (f"e2_fineT_blocks.csv",):
    stamp(Path(__file__).resolve().parents[2] / "data/processed/revision" / _f, __file__,
          ['data/raw/alamode_sto (distributed example outputs; IMPORT_NOTE.md)', '$FEKAPPA_RUNS/fineT_m8_*K (zenodo:prequeue_runs.tar.gz)'])
