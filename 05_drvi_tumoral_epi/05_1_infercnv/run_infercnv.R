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
#   group        ref_tcell / ref_myeloid / stromal / immune_heldout / epi, as prepared
#   subcluster   the per-run leiden partition, "<group>_s<i>" (see analysis_mode below).
#                This is the unit the malignant verdict is taken on: a median over the
#                cells of a subcluster is stable where a per-cell threshold is not, and
#                the partition is computed INSIDE this run, so it is private to this
#                patient and a subcluster is a candidate clone rather than a
#                cross-patient object.
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

# Not cosmetic, and not optional with analysis_mode="subclusters". inferCNV stitches the
# per-subcluster trees through as.phylo()/Newick, and a branch length written as "1e-05"
# is not parseable on the way back in - the run dies after the smoothing has already been
# paid for. Forcing fixed notation avoids it. It changes no result, only how R prints.
options(scipen = 100)

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
  make_option("--analysis-mode", type = "character", default = "subclusters",
              dest = "analysis_mode",
              help = "'subclusters' (default) or 'samples' to reproduce the pre-subcluster runs"),
  make_option("--leiden-resolution", type = "character", default = "0.005",
              dest = "leiden_resolution",
              help = "leiden resolution for the subclustering, or 'auto' [%default]"),
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

cat(sprintf(paste0("cohort   : %s\ninput    : %s\nwork dir : %s\nfig dir  : %s\n",
                   "threads  : %d | HMM: %s | mode: %s | leiden res: %s\n\n"),
            cohort, in_dir, work_dir, fig_dir, opt$threads, opt$hmm,
            opt$analysis_mode, opt$leiden_resolution))

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
#   analysis_mode         "subclusters" (inferCNV's own default; this script used to
#                         override it to "samples"). It partitions each observation group
#                         with leiden INSIDE this run, and that partition is what the
#                         malignant verdict is aggregated over downstream. Two reasons.
#                         First, robustness: the per-cell call has no gap to find - 40% of
#                         the epithelium sits within 1.5x of its own cut - and a median
#                         over the cells of a subcluster is stable where a per-cell
#                         threshold is not. This is also how Shiao et al. take the verdict.
#                         Second, and this is what the previous aggregation got wrong: it
#                         used 04's epithelial leiden, which is computed on the INTEGRATED
#                         object and therefore mixes patients - its cluster 0 held 16,478
#                         cells across 29 cohorts. A CNV profile is private to a patient,
#                         so a cross-patient cluster cannot be a clone and could not carry
#                         the claim being made about it. A per-run partition can.
#   tumor_subcluster_partition_method
#                         "leiden", inferCNV's own default and its documented preference.
#   leiden_resolution     "auto" = (11.98/n_cells)^(1/1.165), so it FALLS as a cohort gets
#                         bigger and the granularity is not comparable across cohorts.
#                         Exposed as --leiden-resolution because that is a property worth
#                         seeing rather than inheriting.
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
  analysis_mode     = opt$analysis_mode,
  tumor_subcluster_partition_method = "leiden",
  leiden_resolution = if (opt$leiden_resolution == "auto") "auto" else
                        as.numeric(opt$leiden_resolution),
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

# --------------------------------------------------------------------------------------
# The subcluster assignment. It lives ONLY in the returned object - nothing writes it to a
# file this script keeps - so it has to be pulled out here, before the working directory is
# removed. inferCNV stores it as a flat named list, "<group>_s<i>" -> integer column
# indices; the character branch is defensive, older layouts stored cell names.
# --------------------------------------------------------------------------------------
subcluster <- rep(NA_character_, ncol(resid))
sc <- infercnv_obj@tumor_subclusters$subclusters
if (is.null(sc)) {
  cat("no subclusters in the object (analysis_mode='samples'?)\n")
} else {
  for (nm in names(sc)) {
    v <- sc[[nm]]
    if (is.list(v)) {                      # nested one level, in some versions
      for (nm2 in names(v)) {
        ix <- v[[nm2]]
        if (is.character(ix)) ix <- match(ix, colnames(resid))
        subcluster[ix] <- nm2
      }
    } else {
      if (is.character(v)) v <- match(v, colnames(resid))
      subcluster[v] <- nm
    }
  }
  cat(sprintf("subclusters: %d over %d cells (%d cells unassigned)\n",
              length(unique(na.omit(subcluster))), ncol(resid), sum(is.na(subcluster))))
  print(table(group, ifelse(is.na(subcluster), "<none>", "assigned")))
}

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
  cell       = colnames(resid),
  cohort     = cohort,
  group      = unname(group),
  subcluster = subcluster,
  cnv_score  = unname(cnv_score),
  cnv_corr   = cnv_corr,
  chr_means,
  check.names = FALSE, stringsAsFactors = FALSE
)
write.csv(out, sum_path, row.names = FALSE)
cat(sprintf("Wrote %s (%d cells)\n", sum_path, nrow(out)))

cat("\nmedian cnv_score by group:\n")
print(round(tapply(out$cnv_score, out$group, median), 5))

# The epithelial subclusters, which are what the verdict will be taken on. A cohort whose
# epithelial subclusters all sit at the level of its stromal ones is a cohort with no
# detectable tumour, and that is a result rather than a run to be retried.
epi_sub <- out[out$group == "epi" & !is.na(out$subcluster), ]
if (nrow(epi_sub)) {
  cat("\nepithelial subclusters (n cells, median cnv_score):\n")
  print(round(do.call(rbind, lapply(split(epi_sub, epi_sub$subcluster), function(d)
    data.frame(n = nrow(d), median_cnv_score = median(d$cnv_score)))), 5))
  cat(sprintf("stromal median for reference: %.5f\n",
              median(out$cnv_score[out$group == "stromal"])))
}

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
