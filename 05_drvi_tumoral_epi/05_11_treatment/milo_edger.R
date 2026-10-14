#!/usr/bin/env Rscript
# 05_11, the GLM half of Milo: one negative-binomial test per neighbourhood, with edgeR.
#
# Called by milo_treatment_tum.py, which builds the neighbourhoods and the counts and does
# the spatial FDR; this script only fits. It is the same model miloR's `testNhoods` fits
# (edgeR quasi-likelihood, TMM, robust dispersion), kept in R because that is where edgeR is
# and re-implementing the empirical-Bayes QL in Python would be a second, untested edgeR.
#
#   Rscript milo_edger.R <workdir>
#
# Reads   <workdir>/nhood_counts.csv   neighbourhoods x biopsies, integer cell counts
#         <workdir>/samples.csv        sample, cohort, treatment, response, n_cells
# Writes  <workdir>/edger_<model>_<contrast>.csv, one per estimable contrast, rows in the
#         input order, columns logFC, logCPM, F, PValue.
#
# Two models, both with the patient as a blocking factor so every contrast is WITHIN-patient:
#
#   all          ~ cohort + treatment                   the main analysis
#   by_response  ~ cohort + group, where group is BASE for every baseline biopsy and
#                <response>_<timepoint> otherwise       exploratory: 4-8 patients a group
#
# In `by_response` the baseline is one level for everybody because the cohort term already
# absorbs each patient's baseline; what is left is the change from it, per response group.
#
# The library size is the number of malignant cells in the biopsy, so a log-fold change reads
# as "this neighbourhood's share of the biopsy's malignant cells", and TMM then corrects for
# the composition of the other neighbourhoods, as in miloR.

suppressPackageStartupMessages({
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: Rscript milo_edger.R <workdir>")
wd <- args[[1]]

counts <- as.matrix(read.csv(file.path(wd, "nhood_counts.csv"), row.names = 1,
                             check.names = FALSE))
samples <- read.csv(file.path(wd, "samples.csv"), stringsAsFactors = FALSE)
rownames(samples) <- samples$sample
samples <- samples[colnames(counts), , drop = FALSE]
stopifnot(identical(rownames(samples), colnames(counts)))

TREATMENTS <- c("BASE", "PD1", "RTPD1")
CONTRASTS <- list(c("PD1", "BASE"), c("RTPD1", "BASE"), c("RTPD1", "PD1"))

samples$cohort <- factor(samples$cohort)
samples$treatment <- factor(samples$treatment,
                            levels = intersect(TREATMENTS, unique(samples$treatment)))
samples$group <- ifelse(samples$treatment == "BASE", "BASE",
                        paste(samples$response, samples$treatment, sep = "_"))
samples$group <- relevel(factor(samples$group), ref = "BASE")

fit_model <- function(design) {
  # A coefficient no biopsy informs (e.g. a response group with no RTPD1 biopsy over the
  # floor) is dropped rather than left to make the design singular; its contrasts are skipped.
  ne <- limma::nonEstimable(design)
  if (!is.null(ne)) {
    message("not estimable, dropped: ", paste(ne, collapse = ", "))
    design <- design[, !colnames(design) %in% ne, drop = FALSE]
  }
  y <- DGEList(counts = counts, lib.size = samples$n_cells)
  y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, design)
  list(fit = glmQLFit(y, design, robust = TRUE), design = design)
}

# The coefficient of one level, or NA when the level is the reference (0) or was dropped.
coef_vector <- function(design, prefix, level, ref) {
  v <- setNames(numeric(ncol(design)), colnames(design))
  if (level == ref) return(v)
  name <- paste0(prefix, level)
  if (!name %in% colnames(design)) return(NULL)
  v[name] <- 1
  v
}

write_contrast <- function(m, model, prefix, later, earlier, ref, tag) {
  a <- coef_vector(m$design, prefix, later, ref)
  b <- coef_vector(m$design, prefix, earlier, ref)
  if (is.null(a) || is.null(b)) {
    message(sprintf("[skip] %s %s: a level is not in the design", model, tag))
    return(invisible(NULL))
  }
  res <- glmQLFTest(m$fit, contrast = a - b)
  tab <- topTags(res, n = Inf, sort.by = "none")$table
  out <- file.path(wd, sprintf("edger_%s_%s.csv", model, tag))
  write.csv(tab[, c("logFC", "logCPM", "F", "PValue")], out)
  message(sprintf("[write] %s  (%d nhoods, %d with p < 0.05)", basename(out), nrow(tab),
                  sum(tab$PValue < 0.05)))
}

# all: ~ cohort + treatment
m <- fit_model(model.matrix(~ cohort + treatment, data = samples))
for (cc in CONTRASTS) {
  write_contrast(m, "all", "treatment", cc[1], cc[2], "BASE",
                 sprintf("%s_vs_%s", cc[1], cc[2]))
}

# by_response: ~ cohort + group
m <- fit_model(model.matrix(~ cohort + group, data = samples))
for (r in sort(unique(samples$response))) {
  for (cc in CONTRASTS) {
    lv <- function(t) if (t == "BASE") "BASE" else paste(r, t, sep = "_")
    write_contrast(m, "by_response", "group", lv(cc[1]), lv(cc[2]), "BASE",
                   sprintf("%s_%s_vs_%s", r, cc[1], cc[2]))
  }
}
