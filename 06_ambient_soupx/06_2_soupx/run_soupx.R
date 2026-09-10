#!/usr/bin/env Rscript
#
# 06_2 step 2: SoupX on one channel, and only the removals kept afterwards.
#
# Reads what prepare_soupx_input.py wrote for one sample, estimates the contamination
# fraction of that channel, and subtracts the ambient counts.
#
# What is kept, per sample:
#   $DATA_DIR/06_amb/adjusted/<sample>/removed.mtx.gz  genes x cells, integer, >= 0
#   $DATA_DIR/06_amb/adjusted/<sample>/rho.csv         one row: the estimate and how it went
#   $DATA_DIR/06_amb/adjusted/<sample>/cells.csv       per cell: UMIs before and after
#   <FIG_DIR>/soupx_<sample>.png                       autoEstCont's own plot, with --plot
#
# **The removals, not the adjusted matrix.** `adjustCounts()` returns a full corrected
# matrix, and writing 148 of those back out would be ~16 GB of Matrix Market on a laptop
# with ~48 GB free (see the arithmetic in prepare_soupx_input.py). What actually changed is
# a much smaller object: SoupX only ever SUBTRACTS, so `toc - adjusted` is non-negative,
# sparse, and enough to reconstruct the result exactly. `assemble_soupx.py` applies it back
# to the counts layer of `shiao.h5ad`. The non-negativity is asserted here rather than
# assumed - it is the property that makes the delta a lossless representation.
#
# **How rho is estimated.** `autoEstCont()` finds genes that are highly specific to some
# clusters (by tf-idf against the soup profile), then asks how much of that gene shows up in
# the clusters that should not express it at all. Its default `tfidfMin = 1.0` can leave a
# channel with no usable marker gene - a small or homogeneous sample, of which this dataset
# has plenty - and the function then stops with an error. Rather than losing the sample, the
# threshold is walked down through TFIDF_STEPS and the value that worked is recorded in
# rho.csv, so a sample estimated at a looser threshold can be told apart in the notebook. If
# none of them work, the channel falls back to FALLBACK_RHO and says so in the `method`
# column; nothing is silently dropped and nothing is silently made up.
#
# `forceAccept = TRUE` is passed on purpose: without it, an estimate landing on the edge of
# SoupX's `contaminationRange` (0.01-0.8) is refused. A channel that really is 80% soup is a
# result this phase exists to find, and the place to look at it is 06_3's notebook, not a
# stack trace here.
#
# Usage (in the soupx-r env, NOT in benchmark-py-r):
#     conda activate soupx-r
#     export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#     Rscript run_soupx.R --sample P01_A_P
#     Rscript run_soupx.R --sample P01_A_P --plot            # save autoEstCont's diagnostic
#     Rscript run_soupx.R --sample P01_A_P --rho 0.10        # skip the estimate, fix rho
#     Rscript run_soupx.R --sample P01_A_P --force           # re-run an existing sample
#
# soupx_all.sh loops this over every prepared sample; this script does one.

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(SoupX)
})

option_list <- list(
  make_option("--sample", type = "character", help = "sample / channel to run (required)"),
  make_option("--data-dir", type = "character", default = Sys.getenv("DATA_DIR"),
              dest = "data_dir", help = "datasets directory [env DATA_DIR]"),
  make_option("--fig-dir", type = "character", default = Sys.getenv("FIG_DIR"),
              dest = "fig_dir", help = "figure directory [env FIG_DIR]"),
  make_option("--method", type = "character", default = "subtraction",
              help = "adjustCounts method: subtraction / soupOnly / multinomial [%default]"),
  make_option("--rho", type = "double", default = NA_real_,
              help = "skip autoEstCont and use this contamination fraction"),
  make_option("--plot", action = "store_true", default = FALSE,
              help = "save autoEstCont's diagnostic plot for this sample"),
  make_option("--force", action = "store_true", default = FALSE,
              help = "re-run a sample whose removals already exist")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$sample)) stop("--sample is required")
if (!nzchar(opt$data_dir)) stop("set DATA_DIR (or pass --data-dir)")

# --------------------------------------------------------------------------------------
# Parameters. Only the first two are judgement calls; the rest are SoupX defaults, named
# here so a run is reproducible from this file alone.
# --------------------------------------------------------------------------------------
TFIDF_STEPS <- c(1.0, 0.8, 0.5, 0.3)  # autoEstCont's marker-specificity threshold, walked down
FALLBACK_RHO <- 0.05                  # SoupX's own priorRho; used only if every step failed
SOUP_QUANTILE <- 0.90                 # autoEstCont default
ROUND_TO_INT <- TRUE                  # 05_2 asserts integer counts; see the phase README
set.seed(0)

sample_id <- opt$sample
data_dir  <- normalizePath(path.expand(opt$data_dir), mustWork = TRUE)
amb_dir   <- file.path(data_dir, "06_amb")
in_dir    <- file.path(amb_dir, "input", sample_id)
out_dir   <- file.path(amb_dir, "adjusted", sample_id)
rho_path  <- file.path(out_dir, "rho.csv")

fig_dir <- if (nzchar(opt$fig_dir)) opt$fig_dir else
  file.path(dirname(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))),
            "figures", "06_2_soupx")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (opt$plot) dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

if (file.exists(rho_path) && !opt$force) {
  cat(sprintf("[have] %s: %s already exists, skipping (--force to re-run)\n", sample_id, rho_path))
  quit(status = 0)
}
if (!dir.exists(in_dir)) stop(sprintf("no prepared input at %s; run prepare_soupx_input.py", in_dir))

cat(sprintf("sample   : %s\ninput    : %s\noutput   : %s\nmethod   : %s\n\n",
            sample_id, in_dir, out_dir, opt$method))

# --------------------------------------------------------------------------------------
# Load. genes.tsv is shared by every sample and lives one level up: it is the row space of
# both counts.mtx.gz and soup.tsv, so all three are indexed the same way by construction.
# --------------------------------------------------------------------------------------
genes <- readLines(file.path(amb_dir, "input", "genes.tsv"))
toc <- readMM(gzfile(file.path(in_dir, "counts.mtx.gz")))
toc <- as(as(toc, "CsparseMatrix"), "dgCMatrix")
cells <- readLines(file.path(in_dir, "barcodes.tsv"))
rownames(toc) <- genes
colnames(toc) <- cells
cat(sprintf("counts: %d genes x %d cells, %d non-zero\n", nrow(toc), ncol(toc), nnzero(toc)))

clusters <- readLines(file.path(in_dir, "clusters.tsv"))
stopifnot(length(clusters) == ncol(toc))
names(clusters) <- cells
cat(sprintf("clusters: %d distinct\n", length(unique(clusters))))

soup <- read.delim(file.path(in_dir, "soup.tsv"), stringsAsFactors = FALSE)
stopifnot(identical(soup$gene, genes))
soup_profile <- data.frame(est = soup$est, counts = soup$counts, row.names = soup$gene)
cat(sprintf("soup: %.0f UMIs over %d genes, top: %s\n\n",
            sum(soup$counts), nrow(soup_profile),
            paste(soup$gene[order(-soup$counts)][1:3], collapse = ", ")))

# --------------------------------------------------------------------------------------
# The channel. `SoupChannel()` normally derives the soup profile from the table of droplets
# (tod), which is the unfiltered Cell Ranger matrix - tens of GB, on the cluster, and the
# reason 06_1 exists. Passing `calcSoupProfile = FALSE` and then `setSoupProfile()` is
# SoupX's own documented path for exactly this situation; `tod = toc` is a placeholder that
# is never read once the profile is set.
# --------------------------------------------------------------------------------------
sc <- SoupChannel(tod = toc, toc = toc, calcSoupProfile = FALSE)
sc <- setSoupProfile(sc, soup_profile)
sc <- setClusters(sc, clusters)

# --------------------------------------------------------------------------------------
# rho. Either given on the command line, or estimated with the threshold walk described in
# the header.
# --------------------------------------------------------------------------------------
rho <- NA_real_
method_used <- NA_character_
tfidf_used <- NA_real_

if (!is.na(opt$rho)) {
  sc <- setContaminationFraction(sc, opt$rho, forceAccept = TRUE)
  rho <- opt$rho
  method_used <- "fixed"
  cat(sprintf("rho fixed at %.4f (--rho), autoEstCont not run\n\n", rho))
} else {
  if (opt$plot) png(file.path(fig_dir, sprintf("soupx_%s.png", sample_id)),
                    width = 1400, height = 1000, res = 150)
  for (tfidf in TFIDF_STEPS) {
    res <- try(autoEstCont(sc, tfidfMin = tfidf, soupQuantile = SOUP_QUANTILE,
                           doPlot = opt$plot, forceAccept = TRUE, verbose = TRUE),
               silent = TRUE)
    if (!inherits(res, "try-error")) {
      sc <- res
      rho <- sc$metaData$rho[1]
      method_used <- "autoEstCont"
      tfidf_used <- tfidf
      cat(sprintf("\nrho = %.4f (autoEstCont, tfidfMin = %.1f)\n\n", rho, tfidf))
      break
    }
    cat(sprintf("[retry] tfidfMin = %.1f failed: %s", tfidf,
                conditionMessage(attr(res, "condition"))))
  }
  if (opt$plot) dev.off()

  if (is.na(rho)) {
    # Every threshold failed. The channel is kept, corrected at SoupX's own prior, and
    # flagged: 06_3's notebook groups the per-cell statistics by `method` precisely so a
    # fallback sample can be looked at separately, or excluded there.
    sc <- setContaminationFraction(sc, FALLBACK_RHO, forceAccept = TRUE)
    rho <- FALLBACK_RHO
    method_used <- "fallback_prior"
    cat(sprintf("\n[WARN] autoEstCont failed at every tfidfMin; falling back to rho = %.4f\n\n",
                FALLBACK_RHO))
  }
}

# --------------------------------------------------------------------------------------
# Adjust, and reduce to the delta. `roundToInt = TRUE` because phase 05 reads this object as
# raw counts and asserts they are integers (05_2/subset_and_qc.ipynb), and because scran,
# DRVI and inferCNV all model counts rather than fractions.
# --------------------------------------------------------------------------------------
adj <- adjustCounts(sc, method = opt$method, roundToInt = ROUND_TO_INT)
removed <- drop0(toc - adj)

stopifnot(nrow(removed) == nrow(toc), ncol(removed) == ncol(toc))
if (min(removed@x, 0) < 0) {
  stop(sprintf("%s: adjustCounts ADDED counts somewhere (min delta %g); the removals are ",
               sample_id, min(removed@x)),
       "not a lossless delta and assemble_soupx.py would be wrong")
}

umis_before <- Matrix::colSums(toc)
umis_after  <- Matrix::colSums(adj)
total_before <- sum(umis_before)
total_after  <- sum(umis_after)
cat(sprintf("removed %.0f of %.0f UMIs (%.2f%%), %d non-zero deltas\n",
            total_before - total_after, total_before,
            100 * (total_before - total_after) / total_before, nnzero(removed)))

# writeMM() takes a file NAME and not a connection, so the gzip is done afterwards by
# streaming the plain file through gzfile(). Keeps the handoff to Python compressed like
# the input, without depending on a `gzip` binary being on PATH.
plain <- tempfile(tmpdir = out_dir, fileext = ".mtx")
invisible(writeMM(removed, plain))   # writeMM prints its NULL return under Rscript
con_in <- file(plain, "rb"); con_out <- gzfile(file.path(out_dir, "removed.mtx.gz"), "wb")
repeat {
  buf <- readBin(con_in, "raw", 1e7L)
  if (length(buf) == 0L) break
  writeBin(buf, con_out)
}
close(con_in); close(con_out); unlink(plain)

write.csv(data.frame(
  cell = cells, umis_before = umis_before, umis_after = umis_after,
  row.names = NULL
), file.path(out_dir, "cells.csv"), row.names = FALSE)

write.csv(data.frame(
  sample = sample_id, rho = rho, method = method_used, tfidf_min = tfidf_used,
  adjust_method = opt$method, n_cells = ncol(toc), n_genes = nrow(toc),
  umis_before = total_before, umis_after = total_after,
  frac_removed = (total_before - total_after) / total_before,
  n_deltas = nnzero(removed), soupx_version = as.character(packageVersion("SoupX"))
), rho_path, row.names = FALSE)

cat(sprintf("\n[ok] %s -> %s\n", sample_id, out_dir))
