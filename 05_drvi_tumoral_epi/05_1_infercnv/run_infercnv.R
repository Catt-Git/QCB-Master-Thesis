#!/usr/bin/env Rscript
#
# 05_1 step 2: inferCNV on one cohort, and only the small summaries kept afterwards.
#
# Reads what prepare_infercnv_input.py wrote for one cohort, runs infercnv::run() with the
# two immune reference groups, and reduces the result to a per-cell table. The heavy part
# of the output is deleted when the run succeeds (--keep-work to keep it), because this
# machine has ~29 GB free and inferCNV's working directory is 1-3 GB PER COHORT: 34 of them
# would not fit, and nothing downstream reads those files.
#
# What is kept, per cohort:
#   $DATA_DIR/05_tum/summary/<cohort>_cnv.csv   one row per cell (see below)
#   <FIG_DIR>/infercnv_<cohort>.png                 the heatmap, the thing to actually look at
#
# The per-cell columns:
#   group        ref_tcell / ref_myeloid / stromal / epi, as prepared
#   cnv_score    mean squared residual across genes. The residual is taken against the mean
#                profile OF THE REFERENCE CELLS OF THIS RUN rather than against the value 1
#                that inferCNV centres on, so the score does not depend on how a given
#                version chose to scale `expr.data`. Reference cells score near zero by
#                construction; that is the point of also carrying `stromal`, which does not.
#   cnv_corr     Pearson correlation of the cell's residual profile with the mean residual
#                profile of the top-CNV epithelial cells (the top TOP_FRAC by cnv_score).
#                This is the second axis of the standard two-axis call (Puram et al. 2017,
#                Neftel et al. 2019): magnitude alone cannot separate a real aneuploid clone
#                from a cell that is noisy everywhere, whereas a cell that is both large and
#                correlated with the cohort's CNV profile is carrying that clone's karyotype.
#   chr1..chr22  mean residual per chromosome, for the figures and for reading which arms
#                drive a call. chrX/chrY/chrM are excluded by inferCNV itself (chr_exclude).
#
# The call itself is NOT made here: thresholds on these two axes are a decision taken by
# looking at the distributions, so they live in call_malignant.ipynb.
#
# Usage (in the infercnv-r env, NOT in benchmark-py-r):
#     conda activate infercnv-r
#     export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#     Rscript run_infercnv.R --cohort Patient64
#     Rscript run_infercnv.R --cohort Patient64 --threads 8 --keep-work
#     Rscript run_infercnv.R --cohort Patient64 --hmm        # + the i6 HMM (slow, needs JAGS)
#
# infercnv_all.sh loops this over every prepared cohort; this script does one.

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(infercnv)
})

option_list <- list(
  make_option("--cohort", type = "character", help = "cohort to run (required)"),
  make_option("--data-dir", type = "character", default = Sys.getenv("DATA_DIR"),
              dest = "data_dir", help = "datasets directory [env DATA_DIR]"),
  make_option("--fig-dir", type = "character", default = Sys.getenv("FIG_DIR"),
              dest = "fig_dir", help = "figure directory [env FIG_DIR]"),
  make_option("--threads", type = "integer", default = 8L,
              help = "threads for infercnv::run() [%default]"),
  make_option("--hmm", action = "store_true", default = FALSE,
              help = "also run the i6 HMM (hours per cohort; not needed for the binary call)"),
  make_option("--keep-work", action = "store_true", default = FALSE, dest = "keep_work",
              help = "keep inferCNV's working directory (1-3 GB per cohort)"),
  make_option("--force", action = "store_true", default = FALSE,
              help = "re-run a cohort whose summary already exists")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$cohort)) stop("--cohort is required")
if (!nzchar(opt$data_dir)) stop("set DATA_DIR (or pass --data-dir)")

# --------------------------------------------------------------------------------------
# Parameters. TOP_FRAC is the only one that is a judgement call rather than an inferCNV
# default; 0.05 is what the two papers above use for the same purpose, and cnv_corr is
# reported rather than thresholded here, so the notebook can see how it behaves.
# --------------------------------------------------------------------------------------
REF_GROUPS <- c("ref_tcell", "ref_myeloid")
CUTOFF <- 0.1          # inferCNV's recommended value for 10x data
TOP_FRAC <- 0.05       # epithelial cells defining the cohort's CNV reference profile
MIN_TOP_CELLS <- 20    # ...but never fewer than this many
GENE_ORDER_NAME <- "gene_order_hg38_gencode_v27.txt"
set.seed(0)

cohort   <- opt$cohort
data_dir <- normalizePath(path.expand(opt$data_dir), mustWork = TRUE)
cnv_dir  <- file.path(data_dir, "05_tum")
in_dir   <- file.path(cnv_dir, "input", cohort)
work_dir <- file.path(cnv_dir, "work", cohort)
sum_dir  <- file.path(cnv_dir, "summary")
sum_path <- file.path(sum_dir, paste0(cohort, "_cnv.csv"))

fig_dir <- if (nzchar(opt$fig_dir)) opt$fig_dir else
  file.path(dirname(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))),
            "figures", "05_1_infercnv")
# infercnv::run() calls dir.create(out_dir) non-recursively, so work/ has to exist first.
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(sum_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

if (file.exists(sum_path) && !opt$force) {
  cat(sprintf("[have] %s: %s already exists, skipping (--force to re-run)\n", cohort, sum_path))
  quit(status = 0)
}
if (!dir.exists(in_dir)) stop(sprintf("no prepared input at %s; run prepare_infercnv_input.py", in_dir))

cat(sprintf("cohort   : %s\ninput    : %s\nwork dir : %s\nfig dir  : %s\nthreads  : %d | HMM: %s\n\n",
            cohort, in_dir, work_dir, fig_dir, opt$threads, opt$hmm))

# --------------------------------------------------------------------------------------
# Load. The counts arrive as Matrix Market (genes x cells) with the names in side files,
# which is what keeps a 12k-cell cohort at tens of MB instead of the ~1 GB a dense
# tab-delimited matrix would take - the format inferCNV documents but this step cannot
# afford 34 times over.
# --------------------------------------------------------------------------------------
counts <- readMM(file.path(in_dir, "counts.mtx"))
counts <- as(as(counts, "CsparseMatrix"), "dgCMatrix")
rownames(counts) <- readLines(file.path(in_dir, "genes.tsv"))
colnames(counts) <- readLines(file.path(in_dir, "barcodes.tsv"))
cat(sprintf("counts: %d genes x %d cells\n", nrow(counts), ncol(counts)))

annot <- read.delim(file.path(in_dir, "annotations.tsv"), header = FALSE,
                    col.names = c("cell", "group"), stringsAsFactors = FALSE)
print(table(annot$group))
stopifnot(all(REF_GROUPS %in% annot$group))

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file  = file.path(in_dir, "annotations.tsv"),
  gene_order_file   = file.path(cnv_dir, GENE_ORDER_NAME),
  ref_group_names   = REF_GROUPS
)
rm(counts); invisible(gc())

# --------------------------------------------------------------------------------------
# Run. The flags that are not defaults, and why:
#   denoise=TRUE          the residuals this script summarises are the denoised ones
#   cluster_by_groups     cluster observations within epi / stromal instead of one global
#                         dendrogram, so the heatmap shows the two blocks separately
#   no_prelim_plot        the preliminary heatmap is a second full-size png of a matrix we
#                         do not use
#   save_rds / save_final_rds / resume_mode / write_expr_matrix = FALSE
#                         all four only exist to write GB-scale intermediates to disk; the
#                         object is returned in memory, which is where the summary is
#                         computed from
#   analysis_mode         left at "samples". 'subclusters' costs hours per cohort and buys
#                         resolution on clonal structure, which is a question 05_1 does not
#                         ask - it asks malignant vs not
# --------------------------------------------------------------------------------------
t0 <- Sys.time()
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff            = CUTOFF,
  out_dir           = work_dir,
  cluster_by_groups = TRUE,
  denoise           = TRUE,
  HMM               = opt$hmm,
  HMM_type          = "i6",
  analysis_mode     = "samples",
  num_threads       = opt$threads,
  no_prelim_plot    = TRUE,
  output_format     = "png",
  write_expr_matrix = FALSE,
  save_rds          = FALSE,
  save_final_rds    = FALSE,
  resume_mode       = FALSE
)
cat(sprintf("\ninfercnv::run() took %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# --------------------------------------------------------------------------------------
# Reduce to the per-cell table.
# --------------------------------------------------------------------------------------
expr <- infercnv_obj@expr.data
ref_idx <- sort(unlist(infercnv_obj@reference_grouped_cell_indices, use.names = FALSE))
stopifnot(length(ref_idx) > 0)

# Residual against the reference mean per gene, so the centring inferCNV applied does not
# enter the score. Reference cells then sit at ~0 by construction.
ref_mean <- rowMeans(expr[, ref_idx, drop = FALSE])
resid <- expr - ref_mean
rm(expr); invisible(gc())

cnv_score <- colMeans(resid^2)

group <- setNames(annot$group, annot$cell)[colnames(resid)]
stopifnot(!anyNA(group))

# cnv_corr: correlate every cell against the mean profile of the most-aneuploid epithelium.
epi_cells <- which(group == "epi")
n_top <- max(MIN_TOP_CELLS, ceiling(TOP_FRAC * length(epi_cells)))
n_top <- min(n_top, length(epi_cells))
top_cells <- epi_cells[order(cnv_score[epi_cells], decreasing = TRUE)[seq_len(n_top)]]
top_profile <- rowMeans(resid[, top_cells, drop = FALSE])
cnv_corr <- as.vector(cor(resid, top_profile))
cat(sprintf("CNV reference profile built on the top %d of %d epithelial cells\n",
            n_top, length(epi_cells)))

# Mean residual per chromosome: small, and the only way to read WHICH arms drive a call.
chrs <- infercnv_obj@gene_order$chr[match(rownames(resid), rownames(infercnv_obj@gene_order))]
chr_levels <- paste0("chr", 1:22)
chr_means <- sapply(chr_levels, function(ch) {
  rows <- which(as.character(chrs) == ch)
  if (length(rows) == 0) rep(NA_real_, ncol(resid)) else colMeans(resid[rows, , drop = FALSE])
})

out <- data.frame(
  cell      = colnames(resid),
  cohort    = cohort,
  group     = unname(group),
  cnv_score = unname(cnv_score),
  cnv_corr  = cnv_corr,
  chr_means,
  check.names = FALSE, stringsAsFactors = FALSE
)
write.csv(out, sum_path, row.names = FALSE)
cat(sprintf("Wrote %s (%d cells)\n", sum_path, nrow(out)))

cat("\nmedian cnv_score by group:\n")
print(round(tapply(out$cnv_score, out$group, median), 5))

# --------------------------------------------------------------------------------------
# Keep the heatmap, drop the rest.
# --------------------------------------------------------------------------------------
heatmap_src <- file.path(work_dir, "infercnv.png")
if (file.exists(heatmap_src)) {
  file.copy(heatmap_src, file.path(fig_dir, sprintf("infercnv_%s.png", cohort)), overwrite = TRUE)
  cat(sprintf("Wrote %s\n", file.path(fig_dir, sprintf("infercnv_%s.png", cohort))))
} else {
  warning(sprintf("no heatmap at %s", heatmap_src))
}

if (opt$keep_work) {
  cat(sprintf("[keep] working directory left at %s\n", work_dir))
} else {
  unlink(work_dir, recursive = TRUE)
  cat(sprintf("[clean] removed %s\n", work_dir))
}
