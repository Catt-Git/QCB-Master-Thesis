# 02_2 integration: method implementations for the R side.
#
# One function per method, sourced by Luecken, M.D. et al. Benchmarking atlas-level data 
# integration in single-cell genomics. Nat Methods (2022). 
# https://doi.org/10.1038/s41592-021-01336-8.
#
# Covers the three R methods:
#
#   fastmnn      (batchelor)  -> corrected RNA assay (full) + fastmnn reduction (embed)
#   seurat_cca   (Seurat)     -> corrected 'integrated' assay (full)
#   seurat_rpca  (Seurat)     -> corrected 'integrated' assay (full)
#
# Notes and differences from the original Luecken et al. 2022 implementations:
#
# Seurat v5 (5.5.1):
#   - `do.cpp = TRUE` was removed from IntegrateData in v5; passing it errors.
#   - runSeuratRPCA's ScaleData()/RunPCA() must receive the HVG *vector*, not the
#     integer 2000 the old default passed.
#
# Legacy anchor API kept: FindIntegrationAnchors + IntegrateData
# return a corrected expression matrix (a `full` output). Seurat v5's
# IntegrateLayers returns an embedding instead, which would turn Seurat into an
# `embed` method and change which metrics are computable.
#
# Reference-based: 34 batches with full pairwise anchoring is not
# feasible at 620k cells, so CCA and RPCA anchor against a few large patients
# (Patient53, Patient16, Patient43 by default). This is a declared deviation from
# Luecken et al. 2022, recorded in Materials & Methods.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

# Map reference patient names to their positions in a SplitObject list, failing
# loudly if any requested reference batch is absent.
.reference_indices <- function(batch_list, reference_batches) {
  idx <- match(reference_batches, names(batch_list))
  if (any(is.na(idx))) {
    missing <- reference_batches[is.na(idx)]
    stop(sprintf("reference batch(es) not found in the data: %s\navailable: %s",
                 paste(missing, collapse = ", "),
                 paste(names(batch_list), collapse = ", ")), call. = FALSE)
  }
  idx
}

# Seurat CCA (reference-based). Returns the integrated Seurat object; the
# corrected matrix is in assay 'integrated'.

runSeuratCCA <- function(sobj, batch, hvg, reference_batches,
                         dims = 1:30, k.weight = 100) {
  batch_list <- SplitObject(sobj, split.by = batch)
  ref_idx <- .reference_indices(batch_list, reference_batches)
  message(sprintf("[cca] %d batches, reference: %s",
                  length(batch_list), paste(reference_batches, collapse = ", ")))

  anchors <- FindIntegrationAnchors(
    object.list     = batch_list,
    anchor.features = hvg,
    reference       = ref_idx,
    reduction       = "cca",
    scale           = TRUE,
    l2.norm         = TRUE,
    dims            = dims,
    k.anchor        = 5,
    k.filter        = 200,
    k.score         = 30,
    max.features    = 200,
    eps             = 0
  )
  IntegrateData(
    anchorset        = anchors,
    new.assay.name   = "integrated",
    dims             = dims,
    k.weight         = k.weight,
    sd.weight        = 1,
    preserve.order   = FALSE,
    eps              = 0,
    verbose          = TRUE
  )
}

# Seurat RPCA (reference-based). Each batch is scaled and PCA'd on the HVGs
# before anchoring with reduction = 'rpca'.

runSeuratRPCA <- function(sobj, batch, hvg, reference_batches,
                          dims = 1:30, k.weight = 100) {
  batch_list <- SplitObject(sobj, split.by = batch)
  ref_idx <- .reference_indices(batch_list, reference_batches)
  message(sprintf("[rpca] %d batches, reference: %s",
                  length(batch_list), paste(reference_batches, collapse = ", ")))

  # features = hvg must be the vector of HVG symbols, not a count.
  batch_list <- lapply(batch_list, function(x) {
    x <- ScaleData(x, features = hvg, verbose = FALSE)
    x <- RunPCA(x, features = hvg, verbose = FALSE)
    x
  })

  anchors <- FindIntegrationAnchors(
    object.list     = batch_list,
    anchor.features = hvg,
    reference       = ref_idx,
    reduction       = "rpca",
    scale           = TRUE,
    l2.norm         = TRUE,
    dims            = dims,
    k.anchor        = 5,
    k.filter        = 200,
    k.score         = 30,
    max.features    = 200,
    eps             = 0
  )
  IntegrateData(
    anchorset        = anchors,
    new.assay.name   = "integrated",
    dims             = dims,
    k.weight         = k.weight,
    sd.weight        = 1,
    preserve.order   = FALSE,
    eps              = 0,
    verbose          = TRUE
  )
}

# fastMNN (batchelor). Returns the Seurat object with the corrected expression in
# the RNA assay (full) and the corrected embedding in the 'fastmnn' reduction
# (embed). Operates on the 'data' layer, which is the log-normalised (or scaled)
# expression prepared in 02_1.

runFastMNN <- function(sobj, batch) {
  suppressPackageStartupMessages({
    library(batchelor)
    library(Matrix)
    library(SingleCellExperiment)
  })

  expr <- GetAssayData(sobj, assay = DefaultAssay(sobj), layer = "data")
  sce <- fastMNN(expr, batch = sobj@meta.data[[batch]])

# Store the reconstructed matrix only in the RNA data slot, leaving counts empty to avoid 
# duplicating a dense matrix and nearly doubling the .rds size. This is also compatible 
# with Seurat 5.5.1.
  corrected <- as(assay(sce, "reconstructed"), "CsparseMatrix")
  sobj[["RNA"]] <- CreateAssayObject(data = corrected)
  DefaultAssay(sobj) <- "RNA"

  # Corrected low-dimensional embedding -> a reduction (embed output).
  emb <- reducedDim(sce, "corrected")
  rownames(emb) <- colnames(sobj)
  sobj[["fastmnn"]] <- CreateDimReducObject(
    embeddings = emb, key = "fastmnn_", assay = "RNA"
  )
  sobj
}

# Registry consumed by run_integration.R.
INTEGRATION_METHODS <- list(
  fastmnn     = function(sobj, batch, hvg, reference_batches)
    runFastMNN(sobj, batch),
  seurat_cca  = function(sobj, batch, hvg, reference_batches)
    runSeuratCCA(sobj, batch, hvg, reference_batches),
  seurat_rpca = function(sobj, batch, hvg, reference_batches)
    runSeuratRPCA(sobj, batch, hvg, reference_batches)
)
