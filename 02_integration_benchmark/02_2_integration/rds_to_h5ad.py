#!/usr/bin/env python
"""
02_2 integration: convert an integrated Seurat .rds back to .h5ad for the metrics.

The R methods (fastmnn, seurat_cca, seurat_rpca) write a Seurat object; the
metrics in 02_4_metrics read only .h5ad. This bridges the two, extracting the right slot
for each output type:

    full  (seurat_cca, seurat_rpca) -> corrected 'integrated' assay -> adata.X
    full  (fastmnn RNA assay)       -> corrected RNA assay          -> adata.X
    embed (fastmnn reduction)       -> 'fastmnn' embedding          -> obsm['X_emb']

Why in python?
The obvious tool, zellkonverter::writeH5AD, has no native R writer: it runs
through basilisk (`basiliskRun(..., testload="anndata")`), which provisions a
private CPython from source on first use. Instead this drives the same
rpy2 + anndata2ri bridge validated elsewhere in the phase: R builds a
SingleCellExperiment, anndata2ri converts it to an AnnData in Python, and anndata
writes it. No basilisk, and it works identically on the cluster.

Cell order is preserved through the conversion (anndata2ri keeps colData row
order), and checked here against the object read back.

Usage:
  export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
  python rds_to_h5ad.py -i $DATA_DIR/02_integration/seurat_cca_unscaled.rds \
                        -o $DATA_DIR/02_integration/seurat_cca_unscaled.h5ad \
                        --types full
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

# scib_compat calls require_r() (fail early if R is unreachable) and installs the
# anndata2ri converter shim used below.
UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import scib_compat  # noqa: F401,E402

import anndata2ri  # noqa: E402
import rpy2.robjects as ro  # noqa: E402
from rpy2.robjects import conversion, default_converter  # noqa: E402


# R side: read the Seurat object and assemble a SingleCellExperiment whose main
# assay is the corrected full matrix and whose 'X_emb' reducedDim, if present, is
# the corrected embedding. Genes are rows (SCE convention); anndata2ri transposes
# to the cells-by-genes AnnData layout.
_BUILD_SCE = """
function(rds_path) {
  suppressPackageStartupMessages({
    library(Seurat); library(SeuratObject); library(SingleCellExperiment)
  })
  so <- readRDS(rds_path)

  # SummarizedExperiment (pulled in by SingleCellExperiment) masks Assays() and
  # would return an S4 object instead of the assay names; qualify the Seurat
  # accessors so the %in% checks below get plain character vectors.
  full_assay <- if ("integrated" %in% SeuratObject::Assays(so)) "integrated" else SeuratObject::DefaultAssay(so)
  X <- SeuratObject::GetAssayData(so, assay = full_assay, layer = "data")

  sce <- SingleCellExperiment(
    assays  = list(X = X),
    colData = so@meta.data
  )

  # Only the fastMNN reduction is a batch-corrected embedding worth exporting.
  # (X_pca / X_umap carried in from 02_1 describe the unintegrated space.)
  if ("fastmnn" %in% SeuratObject::Reductions(so)) {
    emb <- SeuratObject::Embeddings(so, "fastmnn")
    rownames(emb) <- colnames(so)
    reducedDim(sce, "X_emb") <- emb
  }
  sce
}
"""

def parse_args():
    p = argparse.ArgumentParser(description="Convert an integrated Seurat .rds to .h5ad")
    p.add_argument("-i", "--input", required=True, help="integrated Seurat .rds")
    p.add_argument("-o", "--output", required=True, help="output .h5ad")
    p.add_argument("--types", default="full",
                   help="expected output type(s), comma-separated: full, embed "
                        "[default %(default)s]")
    p.add_argument("--emb-out", default=None,
                   help="path to write obsm['X_emb'] as .npy (embed outputs)")
    return p.parse_args()


def main():
    args = parse_args()
    types = {t.strip() for t in args.types.split(",") if t.strip()}

    if not os.path.exists(args.input):
        raise SystemExit(f"input not found: {args.input}")

    build_sce = ro.r(_BUILD_SCE)

    print(f"[read] {args.input}", flush=True)
    with conversion.localconverter(default_converter + anndata2ri.converter):
        adata = build_sce(args.input)

    print(f"[convert] {adata.n_obs:,} cells x {adata.n_vars:,} genes, "
          f"obsm={list(adata.obsm)}", flush=True)

    # Cross-check the object against the declared output types.
    if "embed" in types:
        assert "X_emb" in adata.obsm, (
            "type 'embed' expected but the object has no obsm['X_emb'] "
            "(no fastMNN reduction found)"
        )
    if "full" in types:
        assert adata.X is not None and adata.n_vars > 0, (
            "type 'full' expected but the object has no expression matrix"
        )
    assert np.isfinite(
        adata.X.data if hasattr(adata.X, "data") else np.asarray(adata.X)
    ).all(), "the corrected matrix contains non-finite values"

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    print(f"[write] {args.output}", flush=True)
    adata.write_h5ad(args.output, compression="gzip")

    if args.emb_out and "X_emb" in adata.obsm:
        os.makedirs(os.path.dirname(os.path.abspath(args.emb_out)), exist_ok=True)
        np.save(args.emb_out, np.asarray(adata.obsm["X_emb"], dtype=np.float32))
        print(f"[write] {args.emb_out} ({adata.obsm['X_emb'].shape})", flush=True)

    size_gb = os.path.getsize(args.output) / 1024 ** 3
    print(f"[done] {adata.n_obs:,} x {adata.n_vars:,}, {size_gb:.2f} GB", flush=True)


if __name__ == "__main__":
    main()
