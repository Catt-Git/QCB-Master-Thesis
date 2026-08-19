#!/usr/bin/env Rscript
#
# 02_4 metrics - final summary table
#
# methods ranked by
#
#     overall = 0.6 * bio conservation + 0.4 * batch correction
#
# on min-max scaled metrics, as in Luecken et al. 2022.
#
# By default it reproduces the **published scIB figure** with the vendored
# plotting code in `../utils/` (`plotSingleTaskRNA.R` + `knit_table.R` + the
# `img/` output-type icons, the "funky heatmap": circles per metric, grouped into
# Batch correction / Bio conservation, aggregate score bars, methods ranked). The
# merged CSV is already in the exact shape that code expects: the first column is
# the path `<task>/metrics/<scaling>/<hvg>/<method>_<type>`, which it splits on
# `/` to recover task, scaling, feature space and method.
#
# Usage:
#   Rscript make_summary_table.R -i <merged.csv>
#   Rscript make_summary_table.R -i <merged.csv> -o <outdir> --viz-dir <path>
#   Rscript make_summary_table.R -i <merged.csv> --exclude scgen,scanvi --tag unsup
#
# The last form is the unsupervised-only table: scGen and scANVI are the two
# semi-supervised methods, they see the cell-type labels the bio-conservation
# metrics are computed against, so they are the two starred rows at the top of the
# full figure. `--tag` keeps that run in its own files next to the full one; the
# two tables are NOT comparable row by row, because the scores are means of
# min-max scaled metrics and dropping methods rescales the remaining ones.
#
# TO RUN THIS LOCALLY:
#   #   D=/users/genomics/albertoc/Tesi/hopes_and_dreams/datasets
#   rsync -av albertoc@shiva.prib.upf.edu:$D/02_metrics/ \
#             ~/Desktop/QCB-Master-Thesis/datasets/02_metrics/
#
#   cd ~/Desktop/QCB-Master-Thesis/02_integration_benchmark/02_4_metrics
#   conda activate benchmark-py-r
#   L=~/Desktop/QCB-Master-Thesis/datasets
#   python merge_metrics.py -o $L/02_metrics_merged.csv -r $L/02_metrics \
#          --glob "$L/02_metrics/shiao/metrics/*/hvg/*.csv"
#   Rscript make_summary_table.R -i $L/02_metrics_merged.csv
#
# No `-o` on that last line on purpose: the default is the phase's
# figures/02_4_metrics/, and `-o figures` would write into whatever directory you
# happen to be standing in (which is this one).

suppressPackageStartupMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-i", "--input"), type = "character", help = "merged metrics CSV"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "output directory [default: the phase's figures/02_4_metrics]"),
  make_option("--viz-dir", type = "character", default = NULL,
              dest = "viz_dir", help = "directory with plotSingleTaskRNA.R + knit_table.R + img/"),
  make_option("--weight-batch", type = "double", default = 0.4, dest = "weight_batch",
              help = "weight of the batch-correction score [default %default]"),
  make_option("--exclude", type = "character", default = NULL,
              help = paste("comma-separated method names to drop before scoring,",
                           "e.g. 'scgen,scanvi'. Matched on the <method> token of",
                           "<method>_<type>, so it drops every output type and both",
                           "scalings of that method")),
  make_option("--tag", type = "character", default = NULL,
              help = paste("suffix added to the task label, and so to every output",
                           "file name (<task>_<tag>_summary_*), so a filtered run",
                           "does not overwrite the full one"))
)))
stopifnot(!is.null(opt$input))

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(script_dir) == 0 || script_dir == "") script_dir <- "."
viz_dir <- if (!is.null(opt$viz_dir)) opt$viz_dir else file.path(script_dir, "..", "utils")

# Figures live in ONE directory per phase, `02_integration_benchmark/figures/`,
# with a subdirectory per step - the layout of 01 and of 02_3 (figures/<run_id>).
# So the default is anchored to the script, not to the working directory: a bare
# `-o figures` used to make a second, competing figures/ inside 02_4_metrics/
# whenever the command was run from here, which is where it is run from.
out_dir <- if (!is.null(opt$output)) opt$output else
  file.path(script_dir, "..", "figures", "02_4_metrics")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
opt$output <- out_dir

# Excluding methods is not a display filter: the scores are means of *min-max
# scaled* metrics, so the scaling range - and therefore every score in the table -
# depends on which methods are in it. The rows are dropped before the scorer sees
# them, so the filtered table is rescaled over the methods that remain, and its
# numbers are not comparable with the full table's.
excluded <- if (is.null(opt$exclude)) character(0) else
  trimws(strsplit(opt$exclude, ",")[[1]])
excluded <- excluded[nzchar(excluded)]

# The vendored plotter names its files after the task (the first component of the
# run label) and has no other handle on the output name, so the tag rides along on
# the task label rather than being passed to it.
tag_suffix <- if (is.null(opt$tag) || !nzchar(opt$tag)) "" else paste0("_", opt$tag)


# Metric groups, the split declared in the README, for the fallback scorer.
BATCH_METRICS <- c("PCR_batch", "ASW_label/batch", "iLISI", "graph_conn")
BIO_METRICS   <- c("NMI_cluster/label", "ARI_cluster/label", "ASW_label",
                   "isolated_label_F1", "isolated_label_silhouette", "cLISI",
                   "cell_cycle_conservation", "hvg_overlap")


# the published figure, via the vendored scIB code 
run_scib_plot <- function() {
  plotter <- file.path(viz_dir, "plotSingleTaskRNA.R")
  knit    <- file.path(viz_dir, "knit_table.R")
  if (!file.exists(plotter) || !file.exists(knit)) return(FALSE)

  # plotSingleTaskRNA.R needs the required packages at source time, and does
  # source("knit_table.R") relative to the working directory.
  needed <- c("tibble", "RColorBrewer", "dynutils", "stringr", "Hmisc", "plyr",
              "ggplot2", "cowplot", "dplyr", "ggimage", "scales", "ggforce")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    message("[summary] vendored scIB plotting needs these R packages: ",
            paste(missing, collapse = ", "),
            "\n           install with: Rscript -e \"install.packages(c(",
            paste(sprintf("'%s'", missing), collapse = ", "), "))\"")
    return(FALSE)
  }

  # The scIB parser recovers the method by splitting the path's last component
  # `<method>_<type>` on "_" and taking the first token, so it assumes the method
  # name has no underscore. Two of ours do (seurat_cca, seurat_rpca), which would
  # make it read "cca"/"rpca" as the output type. Normalise only the label handed
  # to the plotter (a display fix; the CSV tree keeps the readable run names).
  fixed <- read.csv(opt$input, check.names = FALSE)

  raw_method <- vapply(as.character(fixed[[1]]), function(lbl) {
    parts <- strsplit(lbl, "/")[[1]]
    toks <- strsplit(parts[length(parts)], "_")[[1]]
    paste(toks[-length(toks)], collapse = "_")
  }, character(1), USE.NAMES = FALSE)

  if (length(excluded) > 0) {
    drop <- raw_method %in% excluded
    unknown <- setdiff(excluded, unique(raw_method))
    if (length(unknown) > 0)
      message("[summary] --exclude names no run: ", paste(unknown, collapse = ", "))
    message("[summary] excluding ", sum(drop), " run(s): ",
            paste(sort(unique(raw_method[drop])), collapse = ", "))
    fixed <- fixed[!drop, , drop = FALSE]
    if (nrow(fixed) == 0) stop("--exclude removed every run")
  }

  fixed[[1]] <- vapply(as.character(fixed[[1]]), function(lbl) {
    parts <- strsplit(lbl, "/")[[1]]
    toks <- strsplit(parts[length(parts)], "_")[[1]]
    type <- toks[length(toks)]
    method <- paste(toks[-length(toks)], collapse = "_")
    method <- switch(method,
                     "seurat_cca" = "seurat",       # -> "Seurat v3 CCA"
                     "seurat_rpca" = "seuratrpca",  # -> "Seurat v3 RPCA"
                     gsub("_", "", method))
    parts[length(parts)] <- paste0(method, "_", type)
    parts[1] <- paste0(parts[1], tag_suffix)
    paste(parts, collapse = "/")
  }, character(1))

  in_abs <- file.path(normalizePath(opt$output),
                      paste0(".plot_input", tag_suffix, ".csv"))
  write.csv(fixed, in_abs, row.names = FALSE)

  out_abs <- normalizePath(opt$output)
  old_wd  <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  # plotSingleTaskRNA.R sources "knit_table.R" and knit_table.R reads the
  # output-type icons from "./img/", both relative to the working directory.
  setwd(viz_dir)
  source("plotSingleTaskRNA.R")
  plotSingleTaskRNA(csv_metrics_path = in_abs, outdir = out_abs,
                    weight_batch = opt$weight_batch)
  message("[summary] wrote the published-style figure(s) to ", out_abs)
  TRUE
}

# fallback: self-contained min-max scorer + plain heatmap 
run_fallback <- function() {
  suppressPackageStartupMessages(library(ggplot2))
  message("[summary] using the built-in scorer (plain heatmap, not the published look).")

  tab <- read.csv(opt$input, row.names = 1, check.names = FALSE)
  run_label <- sub(".*/", "", rownames(tab))

  if (length(excluded) > 0) {
    raw_method <- vapply(strsplit(run_label, "_"), function(toks)
      paste(toks[-length(toks)], collapse = "_"), character(1))
    keep <- !(raw_method %in% excluded)
    message("[summary] excluding ", sum(!keep), " run(s): ",
            paste(sort(unique(raw_method[!keep])), collapse = ", "))
    tab <- tab[keep, , drop = FALSE]
    run_label <- run_label[keep]
    if (nrow(tab) == 0) stop("--exclude removed every run")
  }

  present <- function(cols) cols[cols %in% colnames(tab)]
  batch_cols <- present(BATCH_METRICS)
  bio_cols   <- present(BIO_METRICS)

  minmax <- function(x) {
    if (all(is.na(x))) return(rep(NA_real_, length(x)))
    rng <- range(x, na.rm = TRUE)
    if (!is.finite(rng[1]) || diff(rng) == 0) return(rep(NA_real_, length(x)))
    (x - rng[1]) / diff(rng)
  }
  scaled <- as.data.frame(lapply(tab, minmax), check.names = FALSE)
  rownames(scaled) <- rownames(tab)

  row_mean <- function(df, cols) {
    if (length(cols) == 0) return(rep(NA_real_, nrow(df)))
    rowMeans(df[, cols, drop = FALSE], na.rm = TRUE)
  }
  bio_score   <- row_mean(scaled, bio_cols)
  batch_score <- row_mean(scaled, batch_cols)
  overall <- (1 - opt$weight_batch) * bio_score + opt$weight_batch * batch_score

  scores <- data.frame(run = run_label,
                       bio_conservation = round(bio_score, 4),
                       batch_correction = round(batch_score, 4),
                       overall = round(overall, 4), check.names = FALSE)
  scores <- scores[order(-scores$overall), ]
  write.csv(scores, file.path(opt$output, paste0("summary", tag_suffix, "_scores.csv")),
            row.names = FALSE)

  long <- data.frame(
    run = factor(rep(run_label, ncol(scaled)), levels = rev(scores$run)),
    metric = factor(rep(colnames(scaled), each = nrow(scaled)), levels = c(bio_cols, batch_cols)),
    value = as.vector(as.matrix(scaled)))
  long <- long[long$metric %in% c(bio_cols, batch_cols), ]

  p <- ggplot(long, aes(metric, run, fill = value)) +
    geom_tile(colour = "grey90") +
    scale_fill_viridis_c(na.value = "grey95", limits = c(0, 1)) +
    labs(x = NULL, y = NULL, fill = "min-max\nscaled",
         title = "Integration benchmark (overall = 0.6 bio + 0.4 batch)") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(opt$output, paste0("summary", tag_suffix, "_heatmap.png")), p,
         width = 9, height = 0.4 * nrow(scores) + 2, dpi = 150)
  message("[summary] wrote summary", tag_suffix, "_scores.csv + summary",
          tag_suffix, "_heatmap.png to ", opt$output)
}


if (!run_scib_plot()) run_fallback()
