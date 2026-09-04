"""
05_2 leiden clustering at the optimal resolution (mirrors 04_1/clustering_epi.py).
Produces the definitive malignant object.

Kept separate from reduce_data_tum.py because this is the stochastic, compute-heavy step
and may be rerun on its own.

## The one decision this phase could not inherit

04_1 picks the leiden resolution by maximising NMI against `cell_type`. On the malignant
subset that criterion has no target left: `cell_type` is the post-CNV label, which is the
constant `malignant` by construction, and NMI against a constant is zero at every
resolution - the sweep would return the first grid point and mean nothing.

The target used instead is `cell_type_01_4`, the pre-CNV CellTypist label the malignant
cells were carrying. It is worth being explicit about what that is and is not:

  * it is NOT a tumour taxonomy. It is a normal-breast model's answer to "which healthy
    state does this cell resemble most", which is exactly the contamination this phase
    exists to remove from the ANNOTATION;
  * it is nevertheless not noise, and not constant. Inside the malignant cells it splits
    18,943 `Lumsec-prol`, 12,028 `Lumsec-basal`, 3,975 `LummHR-SCGB` and 1,174
    `Lumsec-KIT`, and those groups differ transcriptionally for real reasons -
    differentiation state and proliferation - even though the names are borrowed from a
    healthy atlas;
  * it is used HERE and only here, to pick one number. The label is not propagated as a
    biological claim, and 05_6/05_7 do not standardise within it (see the note on GROUPBY
    in the phase README): inside one malignant compartment those groups are states, and
    state is the quantity this phase is trying to measure, so regressing it out would
    remove the signal rather than a confounder.

The honest alternative would be to choose the resolution by a criterion internal to the
data (modularity, stability under subsampling). That is a different method than the one 01,
03 and 04 all used, and switching it here would make this phase's clusters incomparable with
theirs for a reason unrelated to the biology. The borrowed label is the smaller compromise,
and it is written down rather than hidden.

Under CELL_SET=epi the target is the post-CNV `cell_type`, which is not constant there, and
this whole discussion does not apply - the script picks the right column either way.

The resolution grid is the full scib 0.1-2.0, as in 04_1.

Input : $DATA_DIR/05_tum/<prefix>_reduced.h5ad
        expects the neighbors graph already computed (from reduce_data_tum.py)
Output: $DATA_DIR/05_tum/<prefix>.h5ad
        adds per-resolution leiden columns 'optscib_tum_leiden_<res>', the selected
        'optscib_tum_leiden' and the resolution/NMI profile in
        .uns['optscib_tum_leiden_profile']
        $DATA_DIR/05_tum/<prefix>_leiden_resolution_profile.csv

Usage:
This is step 4 of 4 (`cluster`) of subsetting_all.sh. Re-running it over an existing output
needs --force.

export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
./subsetting_all.sh                    # the whole chain, resuming
./subsetting_all.sh cluster            # only this step, through the wrapper
./subsetting_all.sh --force cluster    # rerun it over an existing output

Standalone (DATA_DIR exported by hand):
python3 clustering_tum.py
"""

from __future__ import annotations
import numpy as np
import anndata2ri
from rpy2.robjects import conversion, default_converter

def _activate():
    conversion.set_conversion(default_converter + anndata2ri.converter)

def _deactivate():
    conversion.set_conversion(default_converter)

anndata2ri.activate = _activate
anndata2ri.deactivate = _deactivate

import scanpy as sc
import scib

import cell_set as C

np.random.seed(C.SEED)
sc.settings.verbosity = 1

CLUSTER_KEY = "optscib_tum_leiden"

C.banner("05_2 leiden at the optimal resolution")
IN_PATH = C.path("_reduced.h5ad")
OUT_PATH = C.path(".h5ad")
PROFILE_CSV_PATH = C.path("_leiden_resolution_profile.csv")

adata = sc.read_h5ad(IN_PATH)
print(f"Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes", flush=True)

assert "neighbors" in adata.uns, "neighbors graph missing; run reduce_data_tum.py first"

# Pick the NMI target: the post-CNV label when it carries information, the pre-CNV one when
# it does not (the malignant subset). See the module docstring for why this is a documented
# compromise and not a detail.
if adata.obs[C.LABEL_KEY].nunique() > 1:
    label_key = C.LABEL_KEY
    print(f"NMI target: {label_key!r} ({adata.obs[label_key].nunique()} levels)", flush=True)
else:
    label_key = C.NMI_LABEL_KEY
    assert label_key in adata.obs, (
        f"{C.LABEL_KEY!r} is constant on this subset and the fallback target "
        f"{label_key!r} is missing; subset_and_qc.ipynb must carry the pre-CNV label"
    )
    assert adata.obs[label_key].nunique() > 1, (
        f"both {C.LABEL_KEY!r} and {label_key!r} are constant: there is no target to "
        "maximise NMI against, and the resolution cannot be chosen this way"
    )
    print(f"{C.LABEL_KEY!r} is constant ('{adata.obs[C.LABEL_KEY].iloc[0]}') on this subset",
          flush=True)
    print(f"NMI target: {label_key!r} ({adata.obs[label_key].nunique()} levels) "
          "- the pre-CNV CellTypist label, used to pick the resolution only", flush=True)
print(adata.obs[label_key].value_counts().to_string(), flush=True)
print(flush=True)

resolutions = scib.clustering.get_resolutions(n=20)  # [0.1, 0.2, ..., 2.0]
print(f"Resolution grid: {resolutions[0]}..{resolutions[-1]} ({len(resolutions)} values)", flush=True)
res_max, score_max, score_all = scib.clustering.cluster_optimal_resolution(
    adata,
    label_key=label_key,
    cluster_key=CLUSTER_KEY,
    resolutions=resolutions,
    verbose=True,
    return_all=True,
    flavor="igraph",
    n_iterations=2,
)

assert CLUSTER_KEY in adata.obs, "cluster_optimal_resolution did not write the cluster column"
print(f"\nOptimal resolution: {res_max} (NMI {score_max:.4f}), "
      f"{adata.obs[CLUSTER_KEY].nunique()} clusters", flush=True)

# scib returns the profile but stores nothing; it is the evidence behind the chosen
# resolution and visualization_tum.ipynb plots it. Kept in the object and as a CSV.
adata.uns[f"{CLUSTER_KEY}_profile"] = {
    "resolution": score_all["resolution"].to_numpy(),
    "nmi": score_all["score"].to_numpy(),
}
adata.uns[f"{CLUSTER_KEY}_label_key"] = label_key   # which target produced this profile
score_all.to_csv(PROFILE_CSV_PATH, index=False)
print(f"Wrote the resolution profile to {PROFILE_CSV_PATH}", flush=True)

if res_max == resolutions[-1]:
    print(f"WARNING: the optimum is the largest resolution tested ({res_max}); "
          "extend the grid (get_resolutions(n=..., max=...)) before trusting it", flush=True)

# The NMI here is expected to be LOW compared with 04_1's, and a low value is not a failure:
# the target is a borrowed label with four effective levels, not a lineage annotation. What
# would be a failure is a flat profile, i.e. no resolution better than any other.
nmi = score_all["score"].to_numpy()
if np.nanmax(nmi) - np.nanmin(nmi) < 0.01:
    print("WARNING: the NMI profile is flat across the whole grid; the chosen resolution is "
          "effectively arbitrary. Report this rather than treating the optimum as meaningful.",
          flush=True)

adata.write_h5ad(OUT_PATH, compression="gzip")
print(f"Wrote {OUT_PATH}", flush=True)
