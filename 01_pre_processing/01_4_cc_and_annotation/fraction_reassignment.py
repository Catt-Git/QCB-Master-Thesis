"""
01_4 fraction_reassignment: recode `fraction` to a biological immune label from CellTypist

The `fraction` metadata originally stored the technical CD45 sort (CD45+/CD45-), which is
redundant with `dataset_origin` (immune/non_immune -- same sort). Here `fraction` is reassigned
to a BIOLOGICAL immune label taken from the CellTypist annotation: each cell becomes 'imm' if its
`cell_type` is an immune lineage, else 'non_imm'. The technical sort is preserved untouched in
`dataset_origin`, so the sort-vs-lineage mismatch stays inspectable downstream.

Motivation: a LARGE share of cells were annotated by CellTypist as the OPPOSITE lineage of their
CD45 sort - this is not a handful of outliers. In this dataset ~10.7% of the CD45+ (immune-sorted)
cells (~52k / ~490k) get a non-immune cell_type, and ~4.4% of the CD45- (non_immune-sorted) cells
get an immune cell_type. Keeping `fraction` equal to the sort would therefore mislabel tens of
thousands of cells; recoding it from the annotation makes `fraction` reflect biology, while the
raw sort remains available in `dataset_origin` for QC of that very mismatch.

The immune / non-immune lineage sets are specific to the CellTypist model used in
celltypist_annotation.py (Cells_Adult_Breast.pkl, Kumar et al. 2023 adult breast atlas, 58 labels).

Note: 'Lymph-*' are LYMPHATIC ENDOTHELIAL subtypes in this atlas's nomenclature (not lymphocytes),
so they are non-immune -- incl. 'Lymph-immune', an immune-interacting lymphatic endothelial subtype.

Input : $DATA_DIR/all_samples_combined_scrublet_norm_cc_annotated.h5ad
        Output of 01_4_celltypist_annotation. Must contain .obs['cell_type'] and .obs['fraction'].
Output: same file (in-place). Only .obs['fraction'] changes (CD45+/CD45- -> imm/non_imm);
        everything else, including .obs['dataset_origin'], is untouched.

Local usage:
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
python3 fraction_reassignment.py
"""

from __future__ import annotations
import os
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 1

DATA_DIR = os.environ["DATA_DIR"]
IN_PATH = os.path.join(DATA_DIR, "all_samples_combined_scrublet_norm_cc_annotated.h5ad")
OUT_PATH = IN_PATH  # in-place

# CellTypist labels (Cells_Adult_Breast.pkl, 58 total) grouped by lineage.
IMMUNE_CELL_TYPES = {
    # T / NK
    "CD4-Tem", "CD4-Th", "CD4-Th-like", "CD4-Treg", "CD4-activated", "CD4-naive",
    "CD8-Tem", "CD8-Trm", "CD8-activated", "GD", "NK", "NK-ILCs", "NKT", "T_prol",
    # Myeloid
    "Macro-IFN", "Macro-lipo", "Macro-m1", "Macro-m1-CCL", "Macro-m2", "Macro-m2-CXCL",
    "Mast", "Mono-classical", "Mono-non-classical", "Neutrophil",
    "cDC1", "cDC2", "mDC", "pDC", "mye-prol",
    # B / Plasma
    "b_naive", "bmem_switched", "bmem_unswitched", "plasma_IgA", "plasma_IgG",
}
NON_IMMUNE_CELL_TYPES = {
    # Epithelial (luminal hormone-responsive / luminal secretory / basal)
    "LummHR-SCGB", "LummHR-active", "LummHR-major",
    "Lumsec-HLA", "Lumsec-KIT", "Lumsec-basal", "Lumsec-lac", "Lumsec-major",
    "Lumsec-myo", "Lumsec-prol", "basal",
    # Stromal / vascular
    "Fibro-SFRP4", "Fibro-major", "Fibro-matrix", "Fibro-prematrix",
    "Lymph-immune", "Lymph-major", "Lymph-valve1", "Lymph-valve2",
    "Vas-arterial", "Vas-capillary", "Vas-venous", "pericytes", "vsmc",
}

print("Loading data...", flush=True)
adata = sc.read_h5ad(IN_PATH)
print(adata, flush=True)

assert "cell_type" in adata.obs, "Expected CellTypist labels in .obs['cell_type']"
assert "fraction" in adata.obs, "Expected .obs['fraction']"

# Guard: every observed cell_type must be classified as immune or non-immune, otherwise the
# reassignment would silently dump unknown labels into 'non_imm'. Fail loudly instead.
observed = set(adata.obs["cell_type"].astype(str).unique())
unclassified = observed - IMMUNE_CELL_TYPES - NON_IMMUNE_CELL_TYPES
assert not unclassified, (
    f"{len(unclassified)} cell_type value(s) not in either immune/non-immune list; "
    f"update the lineage sets before recoding: {sorted(unclassified)}"
)

is_immune_celltype = adata.obs["cell_type"].isin(IMMUNE_CELL_TYPES)

old_fraction = adata.obs["fraction"].copy()
adata.obs["fraction"] = pd.Categorical(
    np.where(is_immune_celltype, "imm", "non_imm"),
    categories=["imm", "non_imm"],
)

print("\nNew `fraction` counts:", flush=True)
print(adata.obs["fraction"].value_counts(), flush=True)
print("\nOld CD45 sort (fraction) vs new CellTypist-based fraction:", flush=True)
print(pd.crosstab(old_fraction, adata.obs["fraction"]), flush=True)
print("\ndataset_origin (unchanged) vs new fraction:", flush=True)
print(pd.crosstab(adata.obs["dataset_origin"], adata.obs["fraction"]), flush=True)

# Post-conditions: fraction fully recoded, technical sort preserved.
assert set(adata.obs["fraction"].unique()) <= {"imm", "non_imm"}
assert adata.obs["fraction"].notna().all(), "recoded fraction has NaNs"
assert "dataset_origin" in adata.obs, "dataset_origin must be preserved"

# Write atomically: dump to a temp file next to the target, then os.replace. This avoids
# corrupting a large, hard-to-regenerate source file if the write is interrupted mid-way.
print("\nSaving (in-place, atomic)...", flush=True)
tmp_path = OUT_PATH + ".tmp"
adata.write_h5ad(tmp_path, compression="gzip")
os.replace(tmp_path, OUT_PATH)
print(f"Wrote {OUT_PATH}", flush=True)
