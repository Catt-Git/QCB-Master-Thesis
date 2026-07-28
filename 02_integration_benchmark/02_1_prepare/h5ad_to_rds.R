#!/usr/bin/env Rscript

# 02_1_prepare: convert an .h5ad benchmark input to a Seurat v3 .rds
#
# The three R methods (fastMNN, Seurat CCA, Seurat RPCA) read a Seurat object,
# not an AnnData. This produces it, once per input variant, so the integration
# scripts in 02_2 just readRDS() and go.
#
# Why zellkonverter and why reader="R"?
# zellkonverter is the maintained Bioconductor bridge for .h5ad <-> SCE, and
# Seurat converts natively from a SingleCellExperiment. Its native reader
# (reader="R") reads the file directly in R through rhdf5/HDF5Array - no Python,
# so it behaves the same here and on the cluster. 
#
# Why a v3 Seurat object?
# The integration path kept in 02_2 is the legacy FindIntegrationAnchors +
# IntegrateData, which is reliable on a v3 Assay and not on Seurat v5's Assay5.
# options(Seurat.object.assay.version = "v3") forces the classic class.
#
# How the matrix slots are filled:
#   counts slot  <- layers['counts']   (raw integer counts, sparse, non-negative)
#   data slot    <- .X                 (log-normalised, or z-scored for the scaled
#                                        variant: whatever the methods should read)
#
# Seurat v3 requires a non-negative counts slot, so the raw counts always go
# there - including for the scaled input, where .X is negative and could not serve
# as counts. The raw counts stay a sparse dgCMatrix (~0.5 GB), so carrying them
# does not blow up the scaled object, whose cost is the dense .X.
#
# Storage is preserved as the reader hands it over: the unscaled .X is sparse
# (log-normalised, ~9% dense) and stays a dgCMatrix; the scaled .X is genuinely
# dense (per-batch z-scoring turns every zero into a non-zero) and is kept dense.
# Coercing the scaled matrix to sparse would store ~1.24 billion entries plus
# their indices - larger than the dense form it came from.
#
# Two minor differences from the Python reader:
#   - It ignores X_name: .X always lands in an assay literally called "X", which
#     is renamed here.
#   - It does not read .uns (a harmless "H5Identifier not valid" warning). A Seurat
#     object would not carry it anyway; the cell_type palette is a Python-side
#     concern for the figures.
#
# Usage (called once per variant):
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   Rscript h5ad_to_rds.R -i $DATA_DIR/shiao_hvg_2k.h5ad        -o $DATA_DIR/shiao_hvg_2k.rds
#   Rscript h5ad_to_rds.R -i $DATA_DIR/shiao_hvg_2k_scaled.h5ad -o $DATA_DIR/shiao_hvg_2k_scaled.rds

suppressPackageStartupMessages({
  library(optparse)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "input .h5ad"),
  make_option(c("-o", "--output"), type = "character",
              help = "output .rds"),
  make_option("--x-name", type = "character", default = "logcounts",
              help = "name to give the .X assay after reading [default %default]"),
  make_option("--counts-layer", type = "character", default = "counts",
              help = "layer holding raw counts for the Seurat counts slot [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$output)) {
  stop("both -i/--input and -o/--output are required", call. = FALSE)
}
if (!file.exists(opt$input)) {
  stop(sprintf("input not found: %s", opt$input), call. = FALSE)
}

x_name <- opt[["x-name"]]
counts_layer <- opt[["counts-layer"]]

# Read

message(sprintf("[read] %s", opt$input))
t0 <- Sys.time()
sce <- readH5AD(opt$input, reader = "R")
message(sprintf("[read] %.1f s, %d genes x %d cells, assays: %s",
                as.numeric(difftime(Sys.time(), t0, units = "secs")),
                nrow(sce), ncol(sce),
                paste(assayNames(sce), collapse = ", ")))

# The R reader names the .X assay "X" regardless of X_name; rename it.
if (!"X" %in% assayNames(sce)) {
  stop(sprintf("expected an assay named 'X' from the R reader, found: %s",
               paste(assayNames(sce), collapse = ", ")), call. = FALSE)
}
names(assays(sce))[names(assays(sce)) == "X"] <- x_name

has_counts <- counts_layer %in% assayNames(sce)
if (!has_counts) {
  message(sprintf("[warn] no '%s' assay; using %s as the counts slot too",
                  counts_layer, x_name))
}

# Build V3 Seurat object

options(Seurat.object.assay.version = "v3")

# Keep the storage the reader produced: sparse assays become dgCMatrix (what
# Seurat wants), dense assays stay dense. See the note at the top on why the
# scaled .X must not be coerced to sparse.
as_seurat_matrix <- function(m) {
  if (is(m, "sparseMatrix")) as(m, "CsparseMatrix") else as.matrix(m)
}

data_mat <- as_seurat_matrix(assay(sce, x_name))
counts_mat <- if (has_counts) as_seurat_matrix(assay(sce, counts_layer)) else data_mat
message(sprintf("[build] data slot is %s, counts slot is %s",
                if (is(data_mat, "sparseMatrix")) "sparse" else "dense",
                if (is(counts_mat, "sparseMatrix")) "sparse" else "dense"))

message("[build] assembling v3 Seurat object")
so <- CreateSeuratObject(
  counts    = counts_mat,
  meta.data = as.data.frame(colData(sce)),
  assay     = "RNA"
)
so <- SetAssayData(so, layer = "data", new.data = data_mat)

# Carry over the reduced dimensions anndata stored (X_pca, and X_umap on the
# unscaled object), so downstream code can reach them without recomputation.
for (dr in reducedDimNames(sce)) {
  emb <- reducedDim(sce, dr)
  rownames(emb) <- colnames(so)
  key <- paste0(gsub("[^A-Za-z0-9]", "", tolower(dr)), "_")
  so[[dr]] <- CreateDimReducObject(embeddings = emb, key = key, assay = "RNA")
}

# Check the Seurat object is consistent with the SCE it came from
stopifnot(
  "assay is not v3"          = inherits(so[["RNA"]], "Assay"),
  "cell order changed"       = identical(colnames(so), colnames(sce)),
  "gene order changed"       = identical(rownames(so), rownames(sce)),
  "data slot lost its shape" = identical(dim(GetAssayData(so, layer = "data")),
                                         dim(assay(sce, x_name)))
)
if (has_counts) {
  cnts <- GetAssayData(so, layer = "counts")
  stopifnot("counts slot holds non-integers" = all(cnts@x == floor(cnts@x)),
            "counts slot holds negatives"     = min(cnts@x) >= 0)
}

# Write the Seurat object to .rds
# compress=FALSE: the scaled object's data slot is dense; gzip would spend many
# minutes on ~10 GB for little gain, and disk is the cheaper resource here.
message(sprintf("[write] %s", opt$output))
t1 <- Sys.time()
saveRDS(so, opt$output, compress = FALSE)
message(sprintf("[write] %.1f s, %.2f GB on disk",
                as.numeric(difftime(Sys.time(), t1, units = "secs")),
                file.size(opt$output) / 1024^3))
message("[done]")
