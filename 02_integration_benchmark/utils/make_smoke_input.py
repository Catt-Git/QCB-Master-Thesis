"""02 utils: build the small input for the integration smoke test.

The smoke test (see the README, "Smoke test of the integration step") runs every
02_2 dispatcher on a tiny object to check the *code path* - that each method runs,
returns the right output type, and preserves cell order - before the real grid is
launched on 620k cells. This script produces that object from the 5,252-cell
`smoke_fixture.h5ad`:

  smoke_hvg.h5ad       the fixture restricted to the HVGs it contains, with
                       layers['counts'] preserved and written through
                       h5ad_compat so scgen-py (anndata 0.10) can read it too.
  smoke_hvg_list.csv   the HVG symbols actually present in the fixture, one per
                       line. The full 2,000-HVG list from 01_5 contains genes
                       absent from this fixture; passing all 2,000 as Seurat
                       anchor.features makes FindIntegrationMatrix index a gene
                       that is not there and fail with "subscript out of bounds".
                       On the real data all 2,000 are present, so this restriction
                       is a fixture artifact, not part of the benchmark.

Note: `smoke_fixture.h5ad` lists a None-keyed entry in `.layers` that aliases `.X`
in this anndata build; deleting it nulls `.X`. It is left untouched here.

Usage:
  export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
  python make_smoke_input.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import anndata as ad

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h5ad_compat import write_h5ad_compat  # noqa: E402

DATA_DIR = os.environ["DATA_DIR"]
FIXTURE = os.path.join(DATA_DIR, "smoke_fixture.h5ad")
HVG_CSV = os.path.join(DATA_DIR, "shiao_hvg_2k_unintegrated_list.csv")
OUT_H5AD = os.path.join(DATA_DIR, "smoke_hvg.h5ad")
OUT_LIST = os.path.join(DATA_DIR, "smoke_hvg_list.csv")

BATCH_KEY = "cohort"
LABEL_KEY = "cell_type"

print(f"[read] {FIXTURE}", flush=True)
adata = ad.read_h5ad(FIXTURE)
print(f"[read] {adata.n_obs:,} cells x {adata.n_vars:,} genes", flush=True)

assert "counts" in adata.layers, "fixture is missing layers['counts']"
for key in (BATCH_KEY, LABEL_KEY):
    assert key in adata.obs, f"fixture is missing obs[{key!r}]"

hvg = pd.read_csv(HVG_CSV, header=None)[0].astype(str).tolist()
present = [g for g in hvg if g in set(adata.var_names)]
mask = np.asarray(adata.var_names.isin(hvg))
print(f"[hvg] {len(present)} of {len(hvg)} HVGs present in the fixture", flush=True)

obs_before = adata.obs_names.to_numpy().copy()
sub = adata[:, mask].copy()
assert sub.n_vars == len(present), "HVG subset size mismatch"
assert np.array_equal(sub.obs_names.to_numpy(), obs_before), "cell order changed"

counts = sub.layers["counts"]
cvals = counts.data if hasattr(counts, "data") else np.asarray(counts).ravel()
assert np.allclose(cvals, np.floor(cvals)) and cvals.min() >= 0, (
    "layers['counts'] is no longer raw integer counts after the subset"
)

print(f"[write] {OUT_H5AD}", flush=True)
write_h5ad_compat(sub, OUT_H5AD, compression="gzip")

# One symbol per line, no header: the format hvg_csv_to_rds.R reads.
pd.Series(sub.var_names).to_csv(OUT_LIST, index=False, header=False)
print(f"[write] {OUT_LIST} ({len(present)} symbols)", flush=True)

print(
    f"[done] {sub.n_obs:,} x {sub.n_vars:,}, "
    f"{os.path.getsize(OUT_H5AD) / 1024 ** 2:.1f} MB", flush=True,
)
print("Next: h5ad_to_rds.R on smoke_hvg.h5ad, and hvg_csv_to_rds.R on "
      f"smoke_hvg_list.csv with --n-expected {len(present)}.", flush=True)
