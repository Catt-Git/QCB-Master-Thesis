#!/usr/bin/env python3
"""05_10: the standard pipeline against the factor pipeline, on the same cells and the same
gene sets. How many programmes does each one RECOVER?

05_9 asked whether a programme is PRESENT in a coordinate system and answered yes for
Harmony - on most readouts more cleanly than for DRVI. That answer is real and it is the
reason this step exists: presence is not the question the phase turned on. The question is
whether a pipeline, run on these cells with no prior naming of an axis, ARRIVES at the
programme. So this step does not compare two spaces. It compares two ways of turning cells
into a gene list, and then runs the identical ORA on both.

    standard   space -> k-NN graph -> Leiden -> one-vs-rest DE -> top N genes -> ORA
    factor     DRVI decoder -> top N genes per dimension-direction -> ORA

Everything downstream of "top N genes" is one code path. The gene sets, the background, the
list depth, the hypergeometric test and the Benjamini-Hochberg treatment are the same object
in both arms, so the only thing that differs is WHERE THE GENE LIST CAME FROM. That is the
whole design, and it is what lets the difference be attributed.

THE ARMS, AND WHY THERE ARE THREE OF THEM

  harmony_leiden   05_3b's corrected components, clustered      the pipeline anyone would run
  drvi_leiden      05_3's live DRs, clustered the same way      the control that separates
                                                                the SPACE from the UNIT
  drvi_axes        05_3's decoder, read per dimension-direction the phase's own route (05_7)
  pca_leiden       the uncorrected PCA, clustered               optional, what integration bought

`drvi_leiden` is the arm that makes this an experiment rather than a demonstration. If it
lands with `harmony_leiden` and both land far below `drvi_axes`, the loss is a property of
CLUSTERING and not of Harmony - which is the honest reading and the stronger claim, because
it does not require Harmony to be the worse method at anything it promises. If instead
`drvi_leiden` sits above `harmony_leiden`, the space mattered too, and the gap between them
is how much.

THE TWO OUTCOMES, AND WHY THE OBVIOUS ONE IS USELESS

The first thing anyone measures is "did some gene list of this arm enrich for this
signature at FDR < 0.05". On these cells that number SATURATES: it is 17/17 for the
metaprogram collection in both arms, because a 200-gene marker list of a 2,000-gene
background overlaps something in a 22-set catalogue almost always. That saturation is not a
nuisance to be tuned away, it is the first result of the step, and it is reported as
`recall_any` precisely so it cannot be quoted as agreement between the two pipelines.

The outcome that discriminates is the one that matches what naming a dimension MEANS in
05_8:

    a unit of analysis - a cluster or a dimension-direction - is IDENTIFIED BY the
    signature it enriches for most strongly, and a programme counts as NAMED by an arm when
    some unit of that arm has it as its best hit.

`recall_named` is the fraction of the testable catalogue that ends up as some unit's
identity, and it is the headline. The gap between the two numbers is the whole phenomenon:
an arm whose `recall_any` is 1.00 and whose `recall_named` is 0.59 did not find seventeen
programmes, it found ten coarse identities each of which is significant for a dozen
overlapping lists. `mean_hits_per_list` is the same statement from the other side, and
`gene_list_hits` is the table it is counted in.

Neither number is a ranking of integration methods, and neither is charged to Harmony: the
`drvi_leiden` arm runs the identical clustering on DRVI's own space, so whatever the
clustering costs is visible without Harmony in the picture at all.

THE TARGET LIST, AND WHY IT IS NOT "THE PROGRAMMES DRVI NAMED"

Recall needs a denominator that neither pipeline chose. Using 05_8's 19 named programmes
would be circular: they are the output of one of the arms. The denominator here is instead a
property of the CATALOGUE and the OBJECT - every signature of the collection with at least
`C.MIN_SIGNATURE_GENES` genes inside the ORA background, i.e. every set that either arm
could in principle enrich for. It is computed before either arm runs, it is identical for
all of them, and it is printed.

That denominator is generous to the standard pipeline in one direction and harsh in another,
and both are stated rather than corrected: generous because a signature no cell in this
object expresses is in the denominator and neither arm can find it, harsh because a set that
is genuinely absent counts against both equally. `--presence-filter` narrows it to the
signatures 05_9 found present in the Harmony space, when that table is on disk; the headline
number is the unfiltered one.

THE BUDGET OBJECTION, WHICH IS REAL

The factor arm brings 2 x n_latent gene lists to the test and Leiden at resolution 0.2 brings
eight. More lists is more chances to hit, and also a bigger BH denominator. Neither effect is
argued away here; the resolution scan is what answers it empirically. It runs up to
`--resolutions 3.0`, which puts the cluster count in the same range as the direction count,
and every row of the output table carries `n_gene_lists` and `n_tests`, so recall can be read
against the budget that bought it. A standard pipeline that is still below the factor arm at
matched list count has not been starved of tests.

THE NULL, WHICH IS NOT OPTIONAL

A gene list of 200 genes drawn from a 2,000-gene background overlaps a 50-gene signature by
chance often enough that "some cluster was significant for something" is not evidence. The
floor is a SIZE-MATCHED RANDOM PARTITION: the cluster labels of the reference resolution
shuffled among the cells, so the number of groups and every group size are preserved and
nothing else is. It goes through the identical DE and ORA. A random partition recovering
programmes at the same rate as a real one would mean the DE is reading list length rather
than biology, and the step says so.

WHAT THIS STEP DOES NOT CLAIM

It does not rank Harmony and DRVI as integration methods - that is phase 02, on what they
both promise - and it does not say the standard pipeline is wrong. A cluster is the right
unit for a question about discrete populations, and this object has none: it is one
compartment of one lineage, and its structure is continuous by construction. What the step
measures is what that costs, in programmes, on these cells.

Usage:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    cd 05_drvi_tumoral_epi/05_10_pipeline_recall

    N_LATENT=64 python3 pipeline_recall_tum.py --collection gavish
    N_LATENT=64 python3 pipeline_recall_tum.py --collection scie --arms harmony_leiden drvi_axes
    N_LATENT=64 python3 pipeline_recall_tum.py --resolutions 0.4 1.0 2.0 --no-null

`N_LATENT=64` is not optional in practice: 05_3b only ever ran at 64, so at any other value
the Harmony embedding this step needs is not on disk. `CELL_SET`, `HVG_SET` and
`PRUNE_VANISHED` behave exactly as in 05_4 - 05_9.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import warnings

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp
from statsmodels.stats.multitest import multipletests

UTILS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils")
sys.path.insert(0, UTILS_DIR)
import signature_common as C  # noqa: E402
import sig_collections as SC  # noqa: E402

CS = C.CS

STEP = "05_10_pipeline_recall"
FDR = 0.05                       # as 05_7, and applied the same way
K_NEIGHBOURS = 15                # as 05_3, 05_3b and 05_9: the graph is the phase's graph
N_TOP_GENES = 200                # as 05_7: the list depth is shared by both arms
DEFAULT_RESOLUTIONS = (0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0)
REF_RESOLUTION = 1.0
# The scan the DRVI states are read at (`drvi_states`). Separate from DEFAULT_RESOLUTIONS on
# purpose: that grid is shared by every arm of the recall curve and must not move with it.
STATE_RESOLUTIONS = (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)

# The cluster -> direction pairs Alberto kept from the SMI matches at 0.5 (29/09/2026): the
# clean one-to-one ties. The same selection as STATES_SELECTION in 05_8/dr_naming_tum.py,
# which draws them on the naming barplot; `--states-de` characterises these clusters.
STATE_PICKS = {
    "3": ["DR 38+"],
    "7": ["DR 49+", "DR 17+"],
    "10": ["DR 17-"],
    "11": ["DR 55-"],
    "12": ["DR 15+"],
    "15": ["DR 7+"],
    "17": ["DR 45+"],
}
STATE_DE_MIN_LFC = 0.5      # log2 fold change a characterisation marker has to clear
STATE_DE_MIN_PCT = 0.10     # and the share of the cluster's cells that must express it
STATE_DE_MIN_CELLS = 20     # patient x (cluster | rest) cells for the per-patient check
STATE_DE_N_CHECK = 50       # top markers per cluster checked patient by patient
STATE_DR_AUROC = 0.80       # a direction counts as ENRICHED in a cluster from here
STATE_DR_TOP = 3            # directions per cluster carried into the summary

# The max(DR) assignment (`--states-argmax`): every cell goes to the live DR DIRECTION with its
# largest value, on the RAW coordinate - the live dimensions share the prior's scale (sd
# 0.48-0.84), while rescaling each direction to its own q99 hands every direction ~1% of the
# cells, near-dead ones (DR 55-, q99 0.16) included. Two gates turn a winner into an assignment:
ARGMAX_MIN_ACTIVITY = 1.0   # the top value, one prior sd; below it the cell is `inactive`
ARGMAX_MAX_RATIO = 0.8      # second value / top value; above it the cell is `ambiguous`
ARGMAX_MIN_CELLS = 50       # a group smaller than this is kept in the table but not tested

# THE READING of each cluster at resolution 0.5, curated from the --states-de / --states-dr
# outputs of 29/09/2026 (composition, QC, markers, enriched DRs and their names). The factor
# is the ONE thing that best explains the cluster; the numbers on the figure are recomputed
# every run, the words below are not - re-read them if the partition changes.
#   factor: programme | single_patient | depth | contamination | low_quality | none
CLUSTER_READING = {
    "0": ("low_quality", "immediate-early / dissociation stress (FOS, ZFP36); no DR >= 0.8"),
    "1": ("programme", "cell cycle G2/M (UBE2C, CDC20, BIRC5)"),
    "2": ("single_patient", "interferon / MHC (HLA-B, HLA-DRA, LY6E)"),
    "3": ("contamination", "ambient: haemoglobin + fibroblast (HBB, DCN, CFD), 424 genes"),
    "4": ("programme", "secretory, MP23 SECRETED 2 (WFDC2, SLPI) - partly one patient"),
    "5": ("none", "ribosomal-rich, low genes; no programme, no DR >= 0.8"),
    "6": ("single_patient", "Route B only: respiration (MP21)"),
    "7": ("single_patient", "baseline sample; VIM, BMP7, S100B, high genes"),
    "8": ("single_patient", "LIM_STEM axis (DR 32-) in one patient"),
    "9": ("depth", "depth axis DR 2-, 716 genes"),
    "10": ("single_patient", "secretory (PRR27, CSN3, ODAM)"),
    "11": ("too_small", "24 cells"),
    "12": ("contamination", "immune: PTPRC, CCL5, SRGN (doublets / ambient)"),
    "13": ("depth", "depth axis DR 1+, fibroblast ambient (SPARC, DCN)"),
    "14": ("programme", "luminal secretory (MUCL1, APOD, AZGP1) - weak within patient"),
    "15": ("low_quality", "nuclear lncRNA (MALAT1, NEAT1), 735 genes"),
    "16": ("too_small", "25 cells"),
    "17": ("single_patient", "interferon (LY6E, BST2), baseline"),
    "18": ("too_small", "13 cells"),
    "19": ("too_small", "26 cells"),
    "20": ("too_small", "23 cells"),
    "21": ("too_small", "19 cells"),
}
# THE READING of each plain max(DR) group (`--states-argmax`, both gates off), curated from the
# --states-de / --states-dr outputs of 29/09/2026 on the same criteria as CLUSTER_READING, plus
# `few_patients` for a group whose top patient holds 50-80% of it. Median genes per cell over
# the object is 1,383; "genes" below is the group's median. Re-read if the assignment changes.
ARGMAX_PLAIN_READING = {
    "DR 3+": ("programme", "cell cycle G2/M (UBE2C, CDC20, PLK1)"),
    "DR 5+": ("programme", "cell cycle G1/S (TYMS, PCNA, MCM7)"),
    "DR 25-": ("programme", "G2/M, HMG-rich (BIRC5, CDKN3, HMGB3)"),
    "DR 37-": ("programme", "histones / S-G2 (HIST1H4C, DEK, SMC4)"),
    "DR 52-": ("programme", "replication histones (HIST1H1C, HIST1H2BC)"),
    "DR 13-": ("programme", "hypoxia / glycolysis (BNIP3, SLC2A1, LDHA)"),
    "DR 18-": ("programme", "MHC-II (CD74, HLA-DRA, HLA-DPA1)"),
    "DR 23-": ("programme", "IFN-gamma response (GBP1, IDO1, IRF1)"),
    "DR 34-": ("programme", "type I interferon (ISG15, MX1, IFIT1)"),
    "DR 31-": ("programme", "MHC-I (HLA-A/B/C, B2M)"),
    "DR 32+": ("programme", "immunoproteasome (PSMB8, PSMB9, PSME2)"),
    "DR 43-": ("programme", "metallothioneins (MT1X, MT2A, MT1G)"),
    "DR 26-": ("programme", "contractile / EMT_2 (TAGLN, MYL9, ACTG2)"),
    "DR 20+": ("programme", "S100 / EMT_3 (S100A6, S100A10, LGALS1)"),
    "DR 14-": ("programme", "mesenchymal-like (FABP7, VIM, SCRG1)"),
    "DR 42-": ("programme", "YAP targets (CTGF, CYR61, THBS1)"),
    "DR 32-": ("programme", "LIM_STEM (CTNNB1, DSP, CKAP4) - 2 patients"),
    "DR 11-": ("programme", "secretory (FDCSP, WFDC2, SLPI)"),
    "DR 22+": ("programme", "secretory (MGP, CLU, LY6E)"),
    "DR 9-": ("programme", "secretory (WFDC2, SLPI, SMIM22)"),
    "DR 30+": ("programme", "LTF / inflammatory secretory (LTF, CFB, SAA1)"),
    "DR 27+": ("programme", "luminal secretory (AZGP1, PIP, SLPI) - 4 patients"),
    "DR 46-": ("programme", "luminal secretory (MUCL1, AZGP1, APOD) - 4 patients"),
    "DR 19+": ("programme", "luminal progenitor (ELF5, GABRP, S100B)"),
    "DR 28+": ("programme", "epithelial / claudins (CLDN4, TACSTD2, CD24)"),
    "DR 36-": ("programme", "S100A8/A9 senescence (S100A8, S100A9, S100A7)"),
    "DR 40-": ("programme", "basal keratins (KRT5, KRT14, KRT16)"),
    "DR 41-": ("programme", "kallikreins / squamous (KLK5-7, LCN2)"),
    "DR 29-": ("programme", "NF-kB / inflammatory (NFKBIA, CXCL2, TNFAIP3)"),
    "DR 7-": ("programme", "amino-acid stress / ATF4 (ASNS, PSAT1, TRIB3)"),
    "DR 50+": ("programme", "ER stress / UPR (DDIT3, TRIB3, HERPUD1)"),
    "DR 16-": ("programme", "MYC / biosynthesis (TOMM40, GTF3A), 3,029 genes"),
    "DR 9+": ("programme", "metabolic (ATP5F1B, LDHB, PKM)"),
    "DR 10-": ("programme", "ER / collagen folding (SERPINH1, CALR, PDIA3)"),
    "DR 4+": ("single_patient", "ribosomal / respiration (RPL39, COX4I1), 485 genes"),
    "DR 17-": ("single_patient", "secretory (PRR27, CSN3, ODAM)"),
    "DR 49-": ("single_patient", "secretory (CSN3, PRR27, PPP1R1B)"),
    "DR 21-": ("single_patient", "neuroendocrine-like (PCSK1N, PEG10, CRABP1)"),
    "DR 4-": ("single_patient", "MHC-I / ICAM1 (LGMN, SERPING1, ICAM1)"),
    "DR 54+": ("single_patient", "S100B, VIM, SNCA, 4,277 genes"),
    "DR 45+": ("single_patient", "CRABP2, RBP1, PDLIM1, 3,288 genes"),
    "DR 21+": ("few_patients", "PHLDA2, LAMB3, SFN"),
    "DR 35-": ("few_patients", "IFN-gamma chemokines (CXCL9-11, IDO1)"),
    "DR 19-": ("few_patients", "EMT_2 / LIM_STEM (FN1, CD44, CTSC)"),
    "DR 48-": ("few_patients", "secretory (SLPI, WFDC2, FOLR1)"),
    "DR 29+": ("few_patients", "respiration (NDUFA3, CRABP1, TOMM7)"),
    "DR 39+": ("few_patients", "SIVA1, GGCT, CXCL17"),
    "DR 47-": ("few_patients", "SPP1, CD9, PPIB"),
    "DR 51+": ("few_patients", "NPB, PYCARD, BCL2A1"),
    "DR 39-": ("few_patients", "FSCN1, EMP3, PLOD3, 4,111 genes"),
    "DR 16+": ("few_patients", "MHC-I / TF (TF, HLA-C, B2M)"),
    "DR 18+": ("few_patients", "S100 / annexins (S100A10, ANXA2)"),
    "DR 1+": ("depth", "depth axis, fibroblast ambient (DCN, SPARC, COL1A2), 604 genes"),
    "DR 2-": ("contamination", "myeloid (C1QA-C, APOE, TYROBP) on depth axis DR 2-"),
    "DR 38+": ("contamination", "ambient haemoglobin (HBB, HBA1, HBA2), 435 genes"),
    "DR 33+": ("contamination", "plasma-cell ambient (IGKC, JCHAIN, IGHG1), 631 genes"),
    "DR 15+": ("contamination", "immune (PTPRC, CCL5, SRGN)"),
    "DR 8+": ("low_quality", "immediate-early / dissociation (FOS, JUN, EGR1)"),
    "DR 44-": ("low_quality", "heat shock / dissociation (HSPA1A, HSPA1B, DNAJB1)"),
    "DR 53-": ("low_quality", "immediate-early (NR4A1-3, DUSP1, ATF3)"),
    "DR 7+": ("low_quality", "nuclear lncRNA (NEAT1, MALAT1), 769 genes"),
    "DR 6+": ("low_quality", "nuclear lncRNA (MALAT1, NEAT1, XIST), 901 genes"),
    "DR 6-": ("low_quality", "housekeeping only (NDUFA4, H3F3A), 845 genes"),
    "DR 12-": ("none", "mixed (CD81, CTSD, MUC1, JUND)"),
    "DR 24-": ("none", "ribosomal + luminal (RPS4X, XBP1, TFF3)"),
    "DR 14+": ("none", "HSPB1, NUPR1, CRYAB - 3 patients"),
    "DR 10+": ("none", "NUPR1, HSPB1, SOD2, 956 genes"),
    "DR 12+": ("none", "housekeeping (FKBP1A, SDCBP)"),
    "DR 8-": ("none", "housekeeping (SERF2, ALDOA)"),
    "DR 20-": ("none", "3 markers only (EPCAM, HSP90B1)"),
    "DR 22-": ("none", "2 markers only (EEF1A1, CFL1)"),
}

# Slots 1-5 of the validated categorical palette, in fixed order; greys for the two
# categories that are the absence of an explanation.
FACTOR_STYLE = {
    "programme": ("#2a78d6", "biological programme, across patients"),
    "single_patient": ("#eb6834", "single patient (>= 80% of cells)"),
    "few_patients": ("#f5b48f", "few patients (top patient 50-80% of cells)"),
    "depth": ("#1baf7a", "sequencing depth"),
    "contamination": ("#eda100", "contamination (ambient / immune)"),
    "low_quality": ("#e87ba4", "low quality / dissociation stress"),
    "none": ("#9a9a9a", "no clear programme"),
    "too_small": ("#cfcfcf", "too small to read (< 50 cells)"),
}

# What each arm's space IS, in the cache key: a change here drops that arm's cache and no
# other. `drvi_leiden` moved from 64 dimensions to the live ones only.
SPACE_VERSION = {"drvi_leiden": "live_dims_only"}

HEADER_NOTE = ("recall of a signature collection by two pipelines on ONE set of cells: "
               "Leiden clusters vs DRVI decoder directions, identical ORA downstream")

# One arm, one colour, fixed so dropping an arm never recolours the others.
ARM_COLOUR = {
    "harmony_leiden": "#2f6f9f",
    "drvi_leiden": "#7fb3d5",
    "pca_leiden": "#7a7a7a",
    "drvi_axes": "#c25e00",
    "drvi_axes_matched": "#e8993f",
    "null_partition": "#b0b0b0",
}
ARM_LABEL = {
    "harmony_leiden": "Harmony + Leiden + DE",
    "drvi_leiden": "DRVI + Leiden + DE",
    "pca_leiden": "uncorrected PCA + Leiden + DE",
    "drvi_axes": "DRVI decoder directions",
    "drvi_axes_matched": "DRVI decoder, list budget matched",
    "null_partition": "size-matched random partition",
}
LEIDEN_ARMS = ("harmony_leiden", "drvi_leiden", "pca_leiden")
ARM_SPACE = {"harmony_leiden": "harmony", "drvi_leiden": "drvi", "pca_leiden": "pca"}
INK = "#222222"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    SC.add_argument(p)
    p.add_argument("--arms", nargs="+",
                   default=["harmony_leiden", "drvi_leiden", "drvi_axes"],
                   choices=sorted(ARM_LABEL),
                   help="which pipelines to run (default: the three that make the comparison "
                        "an experiment; add pca_leiden for the uncorrected space)")
    p.add_argument("--resolutions", nargs="+", type=float, default=list(DEFAULT_RESOLUTIONS),
                   help="the Leiden resolution scan (default %(default)s). The top of the "
                        "range is there to match the factor arm's list budget")
    p.add_argument("--ref-resolution", type=float, default=REF_RESOLUTION,
                   help="the resolution the per-programme and per-cluster tables are written "
                        "at (default %(default)s); must be in --resolutions")
    p.add_argument("--state-resolutions", nargs="+", type=float,
                   default=list(STATE_RESOLUTIONS),
                   help="Leiden resolutions on the live DRVI dimensions whose clusters get a "
                        "DR heatmap and a size table (default %(default)s); the resolution is "
                        "chosen on the cluster count. Only with drvi_leiden among --arms")
    p.add_argument("--smi-resolution", type=float, default=0.5,
                   help="the state resolution whose clusters get the SMI screening against "
                        "the DR directions (default %(default)s); must be in "
                        "--state-resolutions")
    p.add_argument("--smi-threshold", type=float, default=0.4,
                   help="SMI at or above which a direction-cluster pair is drawn (default "
                        "%(default)s, the 05_3 notebook's value for a leiden target)")
    p.add_argument("--states-only", action="store_true",
                   help="run only the DRVI states (partitions, heatmaps, UMAPs) and stop: "
                        "no DE, no ORA, no recall tables or figures")
    p.add_argument("--states-de", action="store_true",
                   help="characterise the STATE_PICKS clusters at --smi-resolution and stop: "
                        "one-vs-rest DE on all genes, per-patient consistency of the markers, "
                        "composition, marker overlap with the decoder's top genes, ORA")
    p.add_argument("--states-dr", action="store_true",
                   help="every cluster of --smi-resolution: which live DR directions are "
                        "enriched in it (AUROC, cluster vs rest, and within patient), what each "
                        "direction is (name, Route B hits, Hallmark, confounders, top genes), "
                        "and a one-row summary per cluster. Reads the --states-de DE table")
    p.add_argument("--states-argmax", action="store_true",
                   help="assign every cell to the live DR direction with its largest value "
                        "(top 3 kept, activity and dominance gates), then run the --states-de "
                        "and --states-dr analyses on those groups instead of the Leiden "
                        "clusters, and stop")
    p.add_argument("--argmax-min-activity", type=float, default=ARGMAX_MIN_ACTIVITY,
                   help="top value a cell needs to be assigned (default %(default)s)")
    p.add_argument("--argmax-max-ratio", type=float, default=ARGMAX_MAX_RATIO,
                   help="largest second/top value ratio a cell may have to be assigned "
                        "(default %(default)s; 1.0 turns the dominance gate off)")
    p.add_argument("--argmax-min-cells", type=int, default=ARGMAX_MIN_CELLS,
                   help="assigned cells a group needs to be characterised (default "
                        "%(default)s)")
    p.add_argument("--k", type=int, default=K_NEIGHBOURS,
                   help=f"neighbours in the graph (default {K_NEIGHBOURS}, as 05_3/05_3b)")
    p.add_argument("--n-top-genes", type=int, default=N_TOP_GENES,
                   help=f"genes per gene list, BOTH arms (default {N_TOP_GENES}, as 05_7)")
    p.add_argument("--de-fdr", type=float, default=0.05,
                   help="BH-adjusted p a one-vs-rest DE gene has to clear to enter a "
                        "cluster's gene list (default %(default)s). Without it a random "
                        "partition scores like a real one - see de_top_genes")
    p.add_argument("--de-min-lfc", type=float, default=0.0,
                   help="log fold change floor on the same genes (default %(default)s, i.e. "
                        "up-regulated is enough; Seurat's own default is 0.25)")
    p.add_argument("--no-null", action="store_true",
                   help="skip the size-matched random partition, i.e. leave the floor "
                        "unmeasured")
    p.add_argument("--presence-filter", action="store_true",
                   help="narrow the denominator to the signatures 05_9 found present in the "
                        "Harmony space, when that table is on disk")
    p.add_argument("--overwrite", action="store_true",
                   help="recompute the partitions and DE lists instead of reusing the "
                        "cache in $DATA_DIR/05_tum/pipeline_lists_<arm>_<run>.json")
    return p.parse_args()


# --------------------------------------------------------------------------- #
# The spaces
# --------------------------------------------------------------------------- #

def load_coords(name: str, harmony_run: str, drvi_run: str):
    """(cells x dimensions indexed BY CELL NAME, a one-line description) for one space.

    Realignment is by cell name everywhere in this phase and by position nowhere, which is
    the guard 05_3b's own DRVI comparison and 05_9 both use. The uncorrected PCA is the one
    array with no names of its own; it borrows the Harmony embedding's index and is DROPPED
    rather than lined up on a guess if the lengths disagree.
    """
    tum = CS.tum_dir()
    if name in ("harmony", "pca"):
        path = tum / f"embed_{harmony_run}.h5ad"
        if not path.exists():
            sys.exit(f"missing {path}: run 05_3b_harmony_run/run_harmony_tum.py first")
        e = sc.read_h5ad(path)
        cells = e.obs_names.astype(str)
        if name == "harmony":
            X = pd.DataFrame(np.asarray(e.X), index=cells)
            return X, f"05_3b: Harmony, {X.shape[1]} corrected components"
        npy = tum / f"pca_{harmony_run}.npy"
        if not npy.exists():
            print(f"[skip] {npy.name} is not on disk; the uncorrected arm is dropped",
                  flush=True)
            return None, ""
        raw = np.load(npy)
        if raw.shape[0] != len(cells):
            print(f"[skip] {npy.name} has {raw.shape[0]} rows against {len(cells)} cells in "
                  f"the Harmony embedding; the uncorrected arm is dropped", flush=True)
            return None, ""
        X = pd.DataFrame(raw, index=cells)
        return X, f"the same {X.shape[1]} components BEFORE the correction"

    path = tum / f"embed_{drvi_run}.h5ad"
    if not path.exists():
        sys.exit(f"missing {path}: run 05_3_drvi_run first")
    e = sc.read_h5ad(path)
    # The DRs are used the way PCs would be, and a vanished dimension is not one: DRVI's KL
    # pressure collapsed it to ~0.01 sd against a median of 0.6 for the live ones. It would
    # add next to nothing to a Euclidean distance, and it is dropped so the graph is built on
    # exactly the dimensions the heatmaps show. Read from `var['vanished']`, never by hand.
    live = ~e.var["vanished"].astype(bool).to_numpy() if "vanished" in e.var \
        else np.ones(e.n_vars, dtype=bool)
    X = pd.DataFrame(np.asarray(e.X)[:, live], index=e.obs_names.astype(str),
                     columns=e.var["title"].to_numpy()[live] if "title" in e.var else None)
    return X, (f"05_3: DRVI, {X.shape[1]} of {e.n_vars} latent dimensions "
               f"({e.n_vars - X.shape[1]} vanished dropped)")


# --------------------------------------------------------------------------- #
# The two ways of producing a gene list
# --------------------------------------------------------------------------- #

def leiden_labels(adata, coords: pd.DataFrame, resolutions, k: int, seed: int, arm: str):
    """{resolution: labels}, one Leiden partition per resolution on one space.

    The graph is built the way 05_3 and 05_3b build theirs - k = 15, the whole space as the
    representation, no rescaling of the coordinates - so "neighbours" means here what it
    means in those two steps' UMAPs and in 05_9's Moran's I.
    """
    rep = f"X_{arm}"
    adata.obsm[rep] = coords.reindex(adata.obs_names.astype(str)).values.astype(np.float32)
    t0 = time.time()
    sc.pp.neighbors(adata, n_neighbors=k, use_rep=rep, random_state=seed)
    print(f"    graph: k = {k} on {coords.shape[1]} dimensions "
          f"({time.time() - t0:.0f}s)", flush=True)

    out = {}
    for res in resolutions:
        key = f"leiden_{arm}_{res}"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.tl.leiden(adata, resolution=res, key_added=key, flavor="igraph",
                         n_iterations=2, directed=False, random_state=seed)
        lab = adata.obs[key].astype(str)
        out[res] = lab
        print(f"    resolution {res:<5} -> {lab.nunique():3d} clusters "
              f"(smallest {lab.value_counts().min():,}, largest {lab.value_counts().max():,})",
              flush=True)
    return out


def de_top_genes(adata, labels: pd.Series, n_top: int, de_fdr: float, min_lfc: float,
                 key: str = "_grp"):
    """({cluster: [genes]}, per-cluster DE stats) from one-vs-rest Wilcoxon.

    THE SIGNIFICANCE FILTER IS NOT OPTIONAL, and the first run of this step is why. Handing
    the ORA the top `n_top` genes of every cluster REGARDLESS of whether the cluster has any
    marker at all gives a size-matched RANDOM partition the same recall as a real one - it
    scored 10/10 on `scie` against 10/10 for both real arms - because 200 arbitrary genes of
    a 2,000-gene background overlap a 100-gene signature about as often as 200 real ones do.
    That is not a property of clustering, it is the instrument failing to discriminate, and
    it is exactly what the null partition is in this step to catch.

    So a gene enters a cluster's list only if the one-vs-rest test calls it: BH-adjusted
    p < `de_fdr` and log fold change > `min_lfc`, the two filters any real pipeline applies
    before it opens Enrichr. `n_top` then truncates what is left, so a cluster with a
    thousand markers is still read at the same depth as a decoder direction. A cluster with
    NO marker contributes an empty list and tests nothing, which is the correct behaviour
    and is counted: `n_de_genes` is reported per cluster and the arm is not compensated for
    it.

    UP-REGULATED ONLY, and truncated by test statistic rather than by fold change. A
    cluster's ORA question is "what is this cluster's programme", so a gene it is DEPLETED
    of is not part of the answer - the same asymmetry the factor arm has, where the two
    directions of a dimension are read separately and each one's top genes are the ones it
    drives UP.

    The DE runs on the 2,000-gene HVG panel, not on the 25,133-gene object, and that is a
    decision: it is the feature set DRVI was trained on, so it is the only background both
    arms can share. Giving the cluster arm the whole transcriptome would give it genes the
    decoder could never have proposed and a different hypergeometric universe, and the
    comparison would stop being one.
    """
    adata.obs[key] = pd.Categorical(labels.reindex(adata.obs_names.astype(str)).values)
    t0 = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.rank_genes_groups(adata, groupby=key, method="wilcoxon", use_raw=False,
                                n_genes=adata.n_vars, key_added="_de")
        df = sc.get.rank_genes_groups_df(adata, group=None, key="_de")
    keep = (df["logfoldchanges"] > min_lfc) & (df["pvals_adj"] < de_fdr)
    up = df[keep]
    out, stats = {}, []
    for g in adata.obs[key].cat.categories:
        g = str(g)
        grp = up[up["group"].astype(str) == g]
        out[g] = grp.nlargest(n_top, "scores")["names"].tolist()
        stats.append({"cluster": g, "n_de_genes": int(len(grp)),
                      "n_genes_tested": int(len(out[g]))})
    stats = pd.DataFrame(stats)
    empty = int((stats["n_genes_tested"] == 0).sum())
    print(f"    DE: {len(out)} clusters, FDR < {de_fdr} and lfc > {min_lfc}, "
          f"median {stats['n_de_genes'].median():.0f} markers per cluster, "
          f"{empty} with none, truncated at {n_top} ({time.time() - t0:.0f}s)", flush=True)
    return out, stats


def axis_top_genes(embed, gene_names, n_top: int):
    """{dimension-direction: [genes]} off the DRVI decoder. 05_7's own list, rebuilt here.

    Deliberately recomputed rather than read from `C.top_genes_tsv`: the point of this step
    is that both arms go through ONE code path from the gene list onward, and a list read
    from a file written by another run is a second provenance nobody can check from here.
    It is the same call 05_7 makes, so the lists are identical when the run ids match.
    """
    scores = C.interpretability_scores(embed, gene_names, key=C.SCORE_KEY)
    out = {d: scores[d].nlargest(n_top).index.tolist() for d in scores.columns}
    print(f"    decoder: {len(out)} dimension-directions x top {n_top} genes", flush=True)
    return out



# --------------------------------------------------------------------------- #
# The cache: a partition and its DE lists depend on the SPACE, never on the collection
# --------------------------------------------------------------------------- #
#
# Three collections run over the same cells, and the Leiden partitions and their one-vs-rest
# DE are identical in all three - the collection enters only at the ORA. Recomputing them per
# collection is three times the work for a byte-identical result and three chances for two of
# them to disagree about what "resolution 1.0 on the Harmony space" means. The cache is keyed
# on the parameters that CAN change it (k, list depth, the two DE filters) and is discarded
# whole when any of them moves, so a stale list cannot survive a flag change.

def cache_paths(arm: str, run: str):
    tum = CS.tum_dir()
    return (tum / f"pipeline_lists_{arm}_{run}.json",
            tum / f"pipeline_partitions_{arm}_{run}.csv.gz")


def cache_params(args, arm: str, description: str = "") -> dict:
    return {"k": args.k, "n_top_genes": args.n_top_genes, "de_fdr": args.de_fdr,
            "de_min_lfc": args.de_min_lfc, "space": SPACE_VERSION.get(arm, ""),
            "description": description}


def load_cache(arm: str, run: str, args, cells):
    """({resolution: (lists, stats)}, {resolution: labels}, description), or empties.

    The partitions are realigned on cell NAME and the cache is dropped if a single cell of
    this object is missing from it - the same guard every other cross-file read in the phase
    uses, and for the same reason.
    """
    jpath, cpath = cache_paths(arm, run)
    if args.overwrite or not jpath.exists() or not cpath.exists():
        return {}, {}, ""
    try:
        blob = json.loads(jpath.read_text())
    except Exception:                                   # noqa: BLE001 - a corrupt cache is a
        return {}, {}, ""                               # miss, never a crash
    want = cache_params(args, arm)
    have = dict(blob.get("params", {}))
    desc = str(have.pop("description", ""))
    have.setdefault("space", "")                        # caches written before the key
    if {k: v for k, v in want.items() if k != "description"} != have:
        print(f"    [cache] {jpath.name}: parameters changed, recomputing", flush=True)
        return {}, {}, ""
    parts_df = pd.read_csv(cpath, index_col=0)
    parts_df.index = parts_df.index.astype(str)
    if pd.Index(cells).difference(parts_df.index).size:
        print(f"    [cache] {cpath.name}: does not cover these cells, recomputing", flush=True)
        return {}, {}, ""
    lists, parts = {}, {}
    for res_s, payload in blob.get("resolutions", {}).items():
        if res_s not in parts_df.columns:
            continue
        res = float(res_s)
        lists[res] = (payload["lists"], pd.DataFrame(payload["stats"]))
        parts[res] = parts_df[res_s].reindex(cells).astype(str)
    if lists:
        print(f"    [cache] {jpath.name}: {len(lists)} resolutions reused "
              f"({', '.join(str(r) for r in sorted(lists))})", flush=True)
    return lists, parts, desc


def save_cache(arm: str, run: str, args, description: str, lists: dict, parts: dict):
    jpath, cpath = cache_paths(arm, run)
    blob = {"params": cache_params(args, arm, description),
            "resolutions": {str(r): {"lists": l, "stats": st.to_dict(orient="records")}
                            for r, (l, st) in sorted(lists.items())}}
    jpath.write_text(json.dumps(blob))
    pd.DataFrame({str(r): parts[r] for r in sorted(parts)}).to_csv(cpath)
    print(f"    [cache] wrote {jpath.name} and {cpath.name}", flush=True)


# --------------------------------------------------------------------------- #
# The DRVI states: Leiden on the live DRs, drawn as DR heatmaps
# --------------------------------------------------------------------------- #
#
# The DRs used as PCs: the k-NN graph of `drvi_leiden` (live dimensions, k = 15, no
# rescaling), Leiden at `--state-resolutions`, and one DR x cluster heatmap per resolution.
# Same graph and seed as the recall arm, so a resolution in both grids gives the same labels.
#
# CIRCULAR BY CONSTRUCTION. The clusters are regions of the same DR space the heatmap is drawn
# in, so every cluster WILL be high or low on some DRs. The heatmap describes the states as
# combinations of DRs; it is not evidence that the DRs are disentangled or biological. It does
# not depend on the collection, so it is written once, outside the collection folders.

def states_paths(run: str):
    tum = CS.tum_dir()
    fdir = C.OUT_ROOT / f"figures{C.OUT_TAG}" / STEP / "drvi_states" / run
    tdir = C.OUT_ROOT / f"tables{C.OUT_TAG}" / "drvi_states" / run
    fdir.mkdir(parents=True, exist_ok=True)
    tdir.mkdir(parents=True, exist_ok=True)
    return tum / f"drvi_state_partitions_{run}.csv.gz", fdir, tdir


def drvi_states(hvg, args, harmony_run: str, drvi_run: str):
    """Partitions, size table and DR heatmaps of the live-DR Leiden scan."""
    import drvi  # noqa: F401 - only for the heatmap, and only in this function

    C.banner("DRVI states: Leiden on the live DRs, one heatmap per resolution")
    ppath, fdir, tdir = states_paths(drvi_run)
    cells = hvg.obs_names.astype(str)
    parts = pd.DataFrame(index=cells)
    if ppath.exists() and not args.overwrite:
        have = pd.read_csv(ppath, index_col=0)
        have.index = have.index.astype(str)
        if not pd.Index(cells).difference(have.index).size:
            parts = have.reindex(cells)
    todo = [r for r in args.state_resolutions if str(r) not in parts.columns]
    if todo:
        coords, desc = load_coords("drvi", harmony_run, drvi_run)
        print(f"    {desc}", flush=True)
        for res, lab in leiden_labels(hvg, coords, todo, args.k, C.SEED,
                                      "drvi_states").items():
            parts[str(res)] = lab.values
        parts.to_csv(ppath)
        print(f"    wrote {ppath.name}", flush=True)
    else:
        print(f"    [cache] {ppath.name}: every resolution reused", flush=True)

    rows = []
    for res in args.state_resolutions:
        n = parts[str(res)].astype(str).value_counts()
        rows.append({"resolution": res, "n_clusters": len(n), "min_cells": int(n.min()),
                     "median_cells": float(n.median()), "max_cells": int(n.max())})
    sizes = pd.DataFrame(rows)
    tpath = tdir / f"drvi_state_clusters_{drvi_run}.csv"
    sizes.to_csv(tpath, index=False)
    print(f"[table] {tpath}")
    print(sizes.to_string(index=False), flush=True)

    embed = ad.read_h5ad(CS.tum_dir() / f"embed_{drvi_run}.h5ad")
    embed.obs_names = embed.obs_names.astype(str)
    for res in args.state_resolutions:
        key = f"drvi_leiden_{res}"
        lab = parts[str(res)].reindex(embed.obs_names).astype(int)
        embed.obs[key] = pd.Categorical(lab.astype(str),
                                        categories=[str(c) for c in sorted(lab.unique())])
        for sort in (False, True):
            name = f"dr_heatmap_{key}{'_sorted' if sort else ''}_{drvi_run}.png"
            drvi.utils.pl.plot_latent_dims_in_heatmap(embed, key, title_col="title",
                                                      sort_by_categorical=sort, show=False)
            fig = plt.gcf()
            fig.suptitle(f"DRVI {drvi_run}: live DRs by Leiden on the same DRs, resolution "
                         f"{res} ({embed.obs[key].nunique()} clusters)\nclusters built on "
                         f"these DRs - descriptive, circular by construction", fontsize=9,
                         y=1.02)
            fig.savefig(fdir / name, dpi=300, bbox_inches="tight")
            plt.close(fig)
            print(f"[fig] {fdir / name}")

    # The picture is 05_3's UMAP (all 64 dimensions, vanished included), the one every other
    # figure of the phase is drawn on. Only the PICTURE: the clusters are still cut from the
    # 56 live DRs. The exact k-NN graph is identical with or without the vanished dimensions
    # (100% neighbour overlap), so the clusters belong on this layout as much as on any.
    if "X_umap" not in embed.obsm:
        sys.exit(f"embed_{drvi_run}.h5ad has no X_umap: run 05_3's neighbors/umap step")
    umap = (pd.DataFrame(np.asarray(embed.obsm["X_umap"]), index=embed.obs_names)
            .reindex(cells).values)
    print("    UMAP: 05_3's, on all 64 dimensions", flush=True)

    view = ad.AnnData(obs=pd.DataFrame(index=cells))
    view.obsm["X_umap"] = umap
    for res in args.state_resolutions:
        lab = parts[str(res)].astype(int)
        view.obs[f"drvi_leiden_{res}"] = pd.Categorical(
            lab.astype(str), categories=[str(c) for c in sorted(lab.unique())])

    def _panel(ax, res, fontsize):
        key = f"drvi_leiden_{res}"
        sc.pl.umap(view, color=key, ax=ax, show=False, frameon=False,
                   legend_loc="on data", legend_fontsize=fontsize,
                   legend_fontoutline=1.5, size=4,
                   title=f"resolution {res} - {view.obs[key].nunique()} clusters")

    for res in args.state_resolutions:
        fig, ax = plt.subplots(figsize=(7, 6.5))
        _panel(ax, res, 8)
        name = f"umap_drvi_leiden_{res}_{drvi_run}.png"
        fig.savefig(fdir / name, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"[fig] {fdir / name}")

    n = len(args.state_resolutions)
    ncol = 4 if n > 4 else n
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 4.4 * nrow), squeeze=False)
    for ax, res in zip(axes.flat, args.state_resolutions):
        _panel(ax, res, 5)
    for ax in list(axes.flat)[n:]:
        ax.axis("off")
    n_live = embed.n_vars - C.n_vanished(embed)
    fig.suptitle(f"DRVI {drvi_run}: Leiden on the {n_live} live DRs, drawn on 05_3's UMAP "
                 f"(all {embed.n_vars} dimensions)", fontsize=11)
    name = f"umap_drvi_leiden_all_resolutions_{drvi_run}.png"
    fig.savefig(fdir / name, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {fdir / name}")

    if args.smi_resolution is not None:
        drvi_states_smi(embed, view, args.smi_resolution, args.smi_threshold, drvi_run,
                        fdir, tdir)


def drvi_states_smi(embed, view, res: float, threshold: float, drvi_run: str, fdir, tdir):
    """SMI between every live DR direction and the clusters of ONE state resolution.

    The screening of the 05_3 notebook (the DRVI tutorial "Identifying cell types of DRVI
    factors"), unchanged, with the new partition as the target: every dimension split into its
    two directions, each clipped at zero, vanished directions dropped, and the Scaled Mutual
    Information of each direction against the one-hot of each cluster. Pairs at or above
    `threshold` (0.4, the notebook's value for a leiden target) are drawn one at a time.

    Circular in the same way as the heatmaps: the clusters were cut from these directions, so
    a high SMI is expected. What it adds is WHICH direction owns WHICH cluster one to one, and
    which clusters no single direction explains (combinations).
    """
    import networkx as nx
    import drvi
    from drvi.utils.metrics import DiscreteDisentanglementBenchmark

    key = f"drvi_leiden_{res}"
    if key not in view.obs:
        print(f"[skip] SMI: resolution {res} is not among --state-resolutions", flush=True)
        return
    C.banner(f"DRVI states: SMI of the live DR directions against {key}")

    pos = embed[:, ~embed.var["vanished_positive_direction"].astype(bool)].copy()
    neg = embed[:, ~embed.var["vanished_negative_direction"].astype(bool)].copy()
    pos.var.index = pos.var["title"].astype(str) + "+"
    neg.var.index = neg.var["title"].astype(str) + "-"
    pos.X = np.asarray(pos.X).clip(min=0)
    neg.X = -np.asarray(neg.X).clip(max=0)
    directional = pd.concat([pos.to_df(), neg.to_df()], axis=1).loc[view.obs_names]
    print(f"    {directional.shape[1]} live directions x {view.obs[key].nunique()} clusters",
          flush=True)

    target = view.obs[key]
    bench = DiscreteDisentanglementBenchmark(directional.values,
                                             dim_titles=directional.columns,
                                             discrete_target=target, metrics=["SMI"],
                                             aggregation_methods=[])
    bench.evaluate()
    smi = bench.get_results_details()["SMI"]
    smi.index.name = "direction"
    smi = smi[[c for c in target.cat.categories if c in smi.columns]]
    top = (smi.reset_index().melt(id_vars="direction", var_name="cluster", value_name="smi")
           .query("smi >= @threshold").sort_values("smi", ascending=False)
           .reset_index(drop=True))

    smi.to_csv(tdir / f"smi_matrix_{key}_{drvi_run}.csv")
    top.to_csv(tdir / f"smi_matches_{key}_{drvi_run}.csv", index=False)
    print(f"[table] {tdir / f'smi_matrix_{key}_{drvi_run}.csv'}")
    print(f"[table] {tdir / f'smi_matches_{key}_{drvi_run}.csv'}  ({len(top)} pairs >= "
          f"{threshold})")
    best = smi.max(axis=0)
    print("    best SMI per cluster: " + ", ".join(f"{c} {v:.2f}" for c, v in best.items()),
          flush=True)

    # 1. The whole matrix. Clusters as rows in their own order; directions ordered by the
    #    cluster they match best and then by how well, so one-to-one pairs fall on a diagonal
    #    and a cluster carried by several directions shows as a run of cells on its row.
    order = (pd.DataFrame({"best": smi.values.argmax(axis=1), "val": smi.values.max(axis=1)},
                          index=smi.index).sort_values(["best", "val"], ascending=[True, False]))
    mat = smi.loc[order.index].T
    fig, ax = plt.subplots(figsize=(max(12, 0.16 * mat.shape[1]), 0.32 * mat.shape[0] + 2))
    im = ax.imshow(mat.values, aspect="auto", cmap="Reds", vmin=0, vmax=1,
                   interpolation="nearest")
    ax.set_yticks(range(mat.shape[0]), mat.index, fontsize=7)
    ax.set_xticks(range(mat.shape[1]), mat.columns, rotation=90, fontsize=5)
    ax.set_ylabel(f"cluster ({key})")
    ax.set_xlabel("DR direction (live only), ordered by best-matching cluster")
    for (i, j), v in np.ndenumerate(mat.values):
        if v >= threshold:
            ax.text(j, i, f"{v:.2f}".lstrip("0"), ha="center", va="center", fontsize=4,
                    color="white" if v > 0.6 else INK)
    fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01, label="SMI")
    ax.set_title(f"DRVI {drvi_run}: SMI of each DR direction against each cluster, "
                 f"resolution {res}\nvalues shown where SMI >= {threshold}; clusters built on "
                 f"these DRs - circular by construction", fontsize=9)
    name = f"smi_heatmap_{key}_{drvi_run}.png"
    fig.savefig(fdir / name, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {fdir / name}")

    if top.empty:
        print(f"[skip] SMI: no pair >= {threshold}, no network and no pair panels", flush=True)
        return

    # 2. The matches as a bipartite graph, as in the notebook: a component that is not a
    #    single edge is one cluster carried by several directions, or one direction shared.
    G = nx.from_pandas_edgelist(top, "direction", "cluster", edge_attr="smi")
    layout = {}
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    for i, nodes in enumerate(comps):
        sub = nx.spring_layout(G.subgraph(nodes), weight="smi", k=0.5, seed=C.SEED)
        r, c = divmod(i, 5)
        for n, (x, y) in sub.items():
            layout[n] = (x + c * 3, y - r * 3)
    dirs = set(top["direction"])
    w = [d["smi"] for _, _, d in G.edges(data=True)]
    fig = plt.figure(figsize=(16, 3.2 * (int(np.ceil(len(comps) / 5)) + 1)))
    nx.draw(G, layout, with_labels=True, font_size=7, font_weight="bold", node_size=500,
            node_color=["#A0CBE2" if n in dirs else "#FF9E9E" for n in G.nodes()],
            width=[x * 4 for x in w], edge_color=w, edge_cmap=plt.cm.Oranges, alpha=0.7)
    nx.draw_networkx_edge_labels(G, layout, font_size=6,
                                 edge_labels={(u, v): f"{d['smi']:.2f}"
                                              for u, v, d in G.edges(data=True)})
    plt.title(f"SMI >= {threshold}: DR directions (blue) and clusters of {key} (red)",
              fontsize=10)
    plt.axis("off")
    name = f"smi_network_{key}_{drvi_run}.png"
    fig.savefig(fdir / name, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {fdir / name}")

    # 3. Every pair, as the notebook's two-panel figure on the states UMAP: the cluster on
    #    the left, the direction on the right (projection on that half, DRVI's colour range).
    pdir = fdir / f"smi_pairs_{key}"
    pdir.mkdir(exist_ok=True)
    titles = embed.var["title"].astype(str).to_numpy()
    X = np.asarray(embed.X)
    cells = pd.Index(embed.obs_names.astype(str))
    for row in top.itertuples(index=False):
        dim, sign = row.direction[:-1], row.direction[-1]
        j = int(np.flatnonzero(titles == dim)[0])
        s = 1.0 if sign == "+" else -1.0
        vals = pd.Series(s * X[:, j], index=cells).reindex(view.obs_names).to_numpy()
        lo, hi = float(embed.var["min"].iloc[j]), float(embed.var["max"].iloc[j])
        if s < 0:
            lo, hi = -hi, -lo
        panel = ad.AnnData(obs=pd.DataFrame({row.direction: vals}, index=view.obs_names))
        panel.obsm["X_umap"] = view.obsm["X_umap"]
        fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
        sc.pl.umap(view, color=key, groups=[row.cluster], ax=axes[0], show=False,
                   legend_loc="none", na_in_legend=False, size=4, frameon=False,
                   title=f"cluster {row.cluster} ({key})")
        sc.pl.umap(panel, color=row.direction, ax=axes[1], show=False, size=4,
                   frameon=False, cmap=drvi.utils.pl.cmap.saturated_sky_cmap,
                   vmin=min(lo, -1.0), vcenter=0, vmax=max(hi, 1.0),
                   title=f"{row.direction}   (SMI {row.smi:.2f})")
        name = (f"dim_vs_cluster_{row.cluster}_{row.direction.replace(' ', '')}_{key}_"
                f"{drvi_run}.png")
        fig.savefig(pdir / name, dpi=200, bbox_inches="tight")
        plt.close(fig)
    print(f"[fig] {len(top)} pair panels in {pdir}")


def live_directions(embed) -> pd.DataFrame:
    """Cells x live DR directions, each clipped at zero: `DR k+` = max(x, 0), `DR k-` =
    max(-x, 0). Vanished directions dropped, as in the SMI and the AUROC."""
    X = np.asarray(embed.X)
    titles = embed.var["title"].astype(str).to_numpy()
    cols, vals = [], []
    for j, t in enumerate(titles):
        if not bool(embed.var["vanished_positive_direction"].iloc[j]):
            cols.append(f"{t}+"); vals.append(X[:, j].clip(min=0))
        if not bool(embed.var["vanished_negative_direction"].iloc[j]):
            cols.append(f"{t}-"); vals.append((-X[:, j]).clip(min=0))
    return pd.DataFrame(np.column_stack(vals), index=embed.obs_names.astype(str), columns=cols)


def _direction_order(d: str):
    """`DR 3+` < `DR 3-` < `DR 10+`; anything that is not a direction goes last."""
    try:
        return (int(d.split()[1][:-1]), d[-1] == "-")
    except (IndexError, ValueError):
        return (10 ** 6, False)


def argmax_tag(args) -> str:
    """'' for the gated assignment, '_plain' with both gates off (every cell to its top
    direction): the two runs write side by side and never overwrite each other."""
    plain = args.argmax_min_activity <= 0 and args.argmax_max_ratio >= 1.0
    return "_plain" if plain else ""


def argmax_paths(run: str, tag: str = ""):
    return CS.tum_dir() / f"drvi_state_argmax{tag}_{run}.csv.gz"


def drvi_states_argmax(args, drvi_run: str):
    """Every cell to the live DR direction with its largest value: the max(DR) groups.

    A cell's top three directions and their values are kept, and two gates decide whether
    the winner is an assignment or a coin toss:
      inactive   top value < --argmax-min-activity: no direction is on in this cell
      ambiguous  second value > --argmax-max-ratio x top value: two directions are on about
                 equally, and a winner-take-all label would hide the second
    Measured on drvi_tum_64_nomt before choosing the defaults: the median second/top ratio is
    0.80, so the dominance gate is the one that decides (~50% of cells pass it) and the
    activity gate barely moves anything. `argmax_hard` is the ungated winner, kept beside it.

    Unlike Leiden this needs no graph and no resolution, and a group has a name before any
    test is run. It is more circular than Leiden, not less: a group IS the cells where its
    direction wins, so its own direction tops the DR enrichment by construction. What the
    groups add is the DE on a unit defined by one direction, and the directions that
    co-activate with it (the second-ranked ones).
    """
    C.banner("DRVI states: max(DR) assignment on the live directions")
    ppath, fdir, tdir = states_paths(drvi_run)
    tag = argmax_tag(args)
    adir = fdir / f"argmax{tag}"
    adir.mkdir(exist_ok=True)
    embed = ad.read_h5ad(CS.tum_dir() / f"embed_{drvi_run}.h5ad")
    embed.obs_names = embed.obs_names.astype(str)
    V = live_directions(embed)
    cols = np.array(V.columns)
    order = np.argsort(-V.values, axis=1)[:, :3]
    top = np.take_along_axis(V.values, order, axis=1)
    ratio = np.divide(top[:, 1], top[:, 0], out=np.ones(len(top)), where=top[:, 0] > 0)

    status = np.where(top[:, 0] < args.argmax_min_activity, "inactive",
                      np.where(ratio > args.argmax_max_ratio, "ambiguous", "assigned"))
    cells = pd.DataFrame({
        "top1": cols[order[:, 0]], "top2": cols[order[:, 1]], "top3": cols[order[:, 2]],
        "value1": top[:, 0], "value2": top[:, 1], "value3": top[:, 2],
        "ratio_2_to_1": ratio, "status": status, "argmax_hard": cols[order[:, 0]],
    }, index=V.index)
    cells["group"] = np.where(status == "assigned", cells["top1"], "unassigned_" + status)
    n = cells.loc[cells["status"] == "assigned", "group"].value_counts()
    cells["tested"] = cells["group"].map(n).fillna(0).ge(args.argmax_min_cells)
    cells.to_csv(argmax_paths(drvi_run, tag))
    print(f"    {len(V.columns)} live directions; gates: top >= {args.argmax_min_activity}, "
          f"second/top <= {args.argmax_max_ratio}")
    print("    " + ", ".join(f"{k} {v:,} ({v / len(cells):.0%})"
                             for k, v in cells["status"].value_counts().items()))
    print(f"    {len(n)} groups, {int((n >= args.argmax_min_cells).sum())} with >= "
          f"{args.argmax_min_cells} cells (tested); {len(V.columns) - len(n)} directions never "
          f"win a cell")
    print(f"[table] {argmax_paths(drvi_run, tag)}")

    # ---- one row per direction: how many cells it wins, gated and not, and with whom
    rows = []
    hard = cells["argmax_hard"].value_counts()
    for d in sorted(cols, key=_direction_order):
        a = cells[cells["group"] == d]
        co = a["top2"].value_counts(normalize=True)
        rows.append({"direction": d, "n_assigned": len(a), "n_argmax_hard": int(hard.get(d, 0)),
                     "tested": len(a) >= args.argmax_min_cells,
                     "median_value1": float(a["value1"].median()) if len(a) else np.nan,
                     "median_ratio": float(a["ratio_2_to_1"].median()) if len(a) else np.nan,
                     "top_second_direction": co.index[0] if len(co) else "",
                     "top_second_share": float(co.iloc[0]) if len(co) else np.nan})
    groups = pd.DataFrame(rows).set_index("direction")
    groups.to_csv(tdir / f"argmax_groups{tag}_{drvi_run}.csv")
    print(f"[table] {tdir / f'argmax_groups{tag}_{drvi_run}.csv'}")

    # ---- against the Leiden clusters of --smi-resolution: the same cells cut two ways
    if ppath.exists():
        from sklearn.metrics import adjusted_rand_score
        parts = pd.read_csv(ppath, index_col=0)
        parts.index = parts.index.astype(str)
        if str(args.smi_resolution) in parts.columns:
            leid = parts[str(args.smi_resolution)].reindex(cells.index).astype(int).astype(str)
            ct = pd.crosstab(cells["group"], leid)
            ct = ct.reindex(sorted(ct.index, key=_direction_order),
                            columns=sorted(ct.columns, key=int))
            ct.to_csv(tdir / f"argmax{tag}_vs_leiden_{args.smi_resolution}_{drvi_run}.csv")
            m = (cells["status"] == "assigned").values
            print(f"    ARI with drvi_leiden_{args.smi_resolution}, assigned cells only: "
                  f"{adjusted_rand_score(cells['group'][m], leid[m]):.2f}")
            print(f"[table] {tdir / f'argmax{tag}_vs_leiden_{args.smi_resolution}_{drvi_run}.csv'}")

    # ---- figures: the gates, the group sizes, the groups on 05_3's UMAP
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
    axes[0].hist(cells["value1"], bins=80, color="#2a78d6")
    axes[0].axvline(args.argmax_min_activity, color=INK, lw=1, ls="--")
    axes[0].set_xlabel("top direction value (raw DR coordinate)")
    axes[1].hist(cells["ratio_2_to_1"], bins=80, color="#2a78d6")
    axes[1].axvline(args.argmax_max_ratio, color=INK, lw=1, ls="--")
    axes[1].set_xlabel("second value / top value")
    for ax in axes:
        ax.set_ylabel("cells")
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
    st = cells["status"].value_counts(normalize=True)
    fig.suptitle(f"max(DR) gates, {drvi_run}: " + ", ".join(f"{k} {v:.0%}" for k, v in
                                                          st.items()), fontsize=10)
    fig.tight_layout()
    fig.savefig(adir / f"argmax_gates_{drvi_run}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {adir / f'argmax_gates_{drvi_run}.png'}")

    g = groups[groups["n_assigned"] > 0].sort_values("n_assigned", ascending=False)
    fig, ax = plt.subplots(figsize=(0.16 * len(g) + 2, 3.8))
    ax.bar(range(len(g)), g["n_argmax_hard"], color="#d0d0d0", width=0.8,
           label="ungated argmax")
    ax.bar(range(len(g)), g["n_assigned"], width=0.8, label="assigned (both gates)",
           color=np.where(g["tested"], "#2a78d6", "#9a9a9a"))
    ax.axhline(args.argmax_min_cells, color=INK, lw=0.8, ls="--")
    ax.set_yscale("log")
    ax.set_xticks(range(len(g)), g.index, rotation=90, fontsize=5)
    ax.set_xlim(-0.8, len(g) - 0.2)
    ax.set_ylabel("cells (log)")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.legend(fontsize=7, frameon=False)
    ax.set_title(f"cells per max(DR) group; blue = tested (>= {args.argmax_min_cells} "
                 f"assigned cells)", fontsize=9)
    fig.tight_layout()
    fig.savefig(adir / f"argmax_group_sizes_{drvi_run}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {adir / f'argmax_group_sizes_{drvi_run}.png'}")

    if "X_umap" in embed.obsm:
        # Unassigned and untested cells are NaN: drawn grey underneath, with no label on data.
        tested = cells["group"].where(cells["tested"])
        cats = sorted(tested.dropna().unique(), key=_direction_order)
        view = ad.AnnData(obs=pd.DataFrame(index=cells.index))
        view.obsm["X_umap"] = np.asarray(embed.obsm["X_umap"])
        view.obs["group"] = pd.Categorical(tested.values, categories=cats)
        pal = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors) + \
            list(plt.get_cmap("tab20c").colors)
        view.uns["group_colors"] = [matplotlib.colors.to_hex(pal[i % len(pal)])
                                    for i in range(len(cats))]
        fig, ax = plt.subplots(figsize=(8, 7))
        sc.pl.umap(view, color="group", ax=ax, show=False, frameon=False, size=4,
                   legend_loc="on data", legend_fontsize=5, legend_fontoutline=1.5,
                   na_color="#e6e6e6",
                   title=f"max(DR) groups on 05_3's UMAP ({len(cats)} tested)")
        fig.savefig(adir / f"umap_argmax_groups_{drvi_run}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"[fig] {adir / f'umap_argmax_groups_{drvi_run}.png'}")
        view.obs["status"] = pd.Categorical(cells["status"].values)
        fig, ax = plt.subplots(figsize=(7, 6.5))
        sc.pl.umap(view, color="status", ax=ax, show=False, frameon=False, size=4,
                   title="max(DR) gate status")
        fig.savefig(adir / f"umap_argmax_status_{drvi_run}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"[fig] {adir / f'umap_argmax_status_{drvi_run}.png'}")


def state_partition(args, drvi_run: str, argmax: bool) -> dict:
    """The groups `--states-de` and `--states-dr` read: Leiden clusters or max(DR) groups.

      key      the name the tables and figures are written under
      lab      cell -> group, str
      groups   every group, in display order
      picks    the groups characterised by the DE step (Leiden: STATE_PICKS; argmax: the
               tested groups)
      tested   the groups the DE is run on and the DR step reads (Leiden: all of them)
      links    group -> the directions it is tied to (Leiden: the SMI picks; argmax: itself)
    """
    ppath, _, _ = states_paths(drvi_run)
    if argmax:
        tag = argmax_tag(args)
        path = argmax_paths(drvi_run, tag)
        if not path.exists():
            sys.exit(f"missing {path}: run --states-argmax")
        cells = pd.read_csv(path, index_col=0)
        cells.index = cells.index.astype(str)
        lab = cells["group"].astype(str)
        tested = sorted(cells.loc[cells["tested"], "group"].unique(), key=_direction_order)
        return {"key": f"drvi_argmax{tag}", "lab": lab, "unit": "group",
                "groups": sorted(lab.unique(), key=_direction_order),
                "picks": tested, "tested": tested, "links": {g: [g] for g in tested},
                "desc": (f"plain max(DR) groups, every cell (>= {args.argmax_min_cells} cells)"
                         if tag else
                         f"max(DR) groups (top >= {args.argmax_min_activity}, second/top <= "
                         f"{args.argmax_max_ratio}, >= {args.argmax_min_cells} cells)")}
    if not ppath.exists():
        sys.exit(f"missing {ppath}: run with --states-only first")
    parts = pd.read_csv(ppath, index_col=0)
    parts.index = parts.index.astype(str)
    res = args.smi_resolution
    if str(res) not in parts.columns:
        sys.exit(f"resolution {res} not in {ppath.name}: have {list(parts.columns)}")
    lab = parts[str(res)].astype(int).astype(str)
    groups = [str(c) for c in sorted(lab.astype(int).unique())]
    return {"key": f"drvi_leiden_{res}", "lab": lab, "unit": "cluster", "groups": groups,
            "picks": list(STATE_PICKS), "tested": groups, "links": STATE_PICKS,
            "desc": f"Leiden on the live DRs, resolution {res}"}


def drvi_states_de(hvg, args, drvi_run: str, argmax: bool = False):
    """Characterise the STATE_PICKS clusters: what is in them, and is it their DR's programme.

    1. DE, one cluster against every other cell, Wilcoxon, on ALL genes of the log-normalised
       05_2 object (genes expressed in >= 1% of cells). This is the characterisation, so the
       whole transcriptome - unlike `de_top_genes`, whose HVG-only design serves the recall
       comparison. Every cluster of the resolution is tested, so the table is complete; the
       figures and summaries read the picked ones.
    2. PER-PATIENT CONSISTENCY of the top markers. A Leiden cluster in this object can be one
       patient (the Harmony/PCA clusters of this step were, NMI 0.81), and then its "markers"
       are that patient's genes. For each marker, in each patient with >= STATE_DE_MIN_CELLS
       cells both inside and outside the cluster: mean inside minus mean outside. A marker
       that is up in most patients is the state's; one up in one patient is not.
    3. COMPOSITION: patient, treatment, response, phase, plus depth / mt / doublet / CNV.
    4. THE LINK BACK TO THE DR, from the gene side. The cluster's markers, restricted to the
       HVG panel and cut at the top N like a decoder list, against the top-N decoder genes of
       every live direction (hypergeometric, HVG background). The SMI tied the direction to
       the cluster through the CELLS; this asks whether the direction's GENES are the
       cluster's markers - the same two-route logic as 05_8, on a data-driven unit.
    5. ORA of the same HVG-restricted list against scie, emt, gavish_tnbc and Hallmark, BH per
       library across the picked clusters.

    With `argmax` the same five steps run on the max(DR) groups (`state_partition`): every
    tested group is a pick, linked to its own direction, and the DE tests only those groups -
    unassigned cells and groups below --argmax-min-cells stay in "rest".
    """
    from scipy.stats import hypergeom

    P = state_partition(args, drvi_run, argmax)
    key, picks, links, unit = P["key"], P["picks"], P["links"], P["unit"]
    _, fdir, tdir = states_paths(drvi_run)
    ddir = fdir / f"state_de_{key}"
    ddir.mkdir(exist_ok=True)
    C.banner(f"DRVI states: characterising {len(picks)} {unit}s of {key}")

    # ---------------------------------------------------------------- the cells
    full = ad.read_h5ad(C.FULL_H5AD)
    full.obs_names = full.obs_names.astype(str)
    embed = ad.read_h5ad(CS.tum_dir() / f"embed_{drvi_run}.h5ad", backed="r")
    meta = embed.obs.copy()
    meta.index = meta.index.astype(str)
    embed.file.close()
    meta = meta.reindex(full.obs_names)
    lab = P["lab"].reindex(full.obs_names).astype(str)
    full.obs[key] = pd.Categorical(lab.values, categories=P["groups"])
    sc.pp.filter_genes(full, min_cells=int(0.01 * full.n_obs))
    print(f"    {full.n_obs:,} cells x {full.n_vars:,} genes (expressed in >= 1% of cells)",
          flush=True)

    # ---------------------------------------------------------------- 1. DE
    t0 = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.rank_genes_groups(full, groupby=key, groups=P["tested"], method="wilcoxon",
                                use_raw=False, reference="rest", pts=True, key_added="_de")
    de = sc.get.rank_genes_groups_df(full, group=None, key="_de")
    de["group"] = de["group"].astype(str)
    de = de.rename(columns={"group": "cluster", "names": "gene"})
    print(f"    DE: {len(P['tested'])} {unit}s vs rest in {time.time() - t0:.0f}s", flush=True)
    sig = de[(de["pvals_adj"] < 0.05) & (de["logfoldchanges"].abs() > STATE_DE_MIN_LFC)]
    up = sig[(sig["logfoldchanges"] > 0) & (sig["pct_nz_group"] >= STATE_DE_MIN_PCT)]
    up = up.sort_values(["cluster", "scores"], ascending=[True, False])
    sig.to_csv(tdir / f"state_de_{key}_{drvi_run}.csv.gz", index=False)
    print(f"[table] {tdir / f'state_de_{key}_{drvi_run}.csv.gz'}  (all clusters, "
          f"|lfc| > {STATE_DE_MIN_LFC}, FDR < 0.05)")

    # ---------------------------------------------------------------- 2. per patient
    cohort = meta["cohort"].astype(str)
    rows = []
    for cl in picks:
        genes = up[up["cluster"] == cl]["gene"].head(STATE_DE_N_CHECK).tolist()
        if not genes:
            continue
        X = full[:, genes].X
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        inside = (lab == cl).values
        n_ok = 0
        deltas = []
        for pat in cohort.unique():
            pm = (cohort == pat).values
            a, b = pm & inside, pm & ~inside
            if a.sum() < STATE_DE_MIN_CELLS or b.sum() < STATE_DE_MIN_CELLS:
                continue
            n_ok += 1
            deltas.append(X[a].mean(0) - X[b].mean(0))
        deltas = np.array(deltas) if deltas else np.zeros((0, len(genes)))
        for j, g in enumerate(genes):
            rows.append({"cluster": cl, "gene": g, "rank": j + 1, "n_patients_tested": n_ok,
                         "n_patients_up": int((deltas[:, j] > 0).sum()) if n_ok else 0})
    cons = pd.DataFrame(rows)
    cons["share_patients_up"] = cons["n_patients_up"] / cons["n_patients_tested"].replace(0, np.nan)
    marker_tbl = up[up["cluster"].isin(picks)].merge(cons, on=["cluster", "gene"], how="left")
    marker_tbl.to_csv(tdir / f"state_markers_{key}_{drvi_run}.csv", index=False)
    print(f"[table] {tdir / f'state_markers_{key}_{drvi_run}.csv'}")

    # ---------------------------------------------------------------- 3. composition
    comp = []
    for cl in picks:
        m = (lab == cl).values
        pat = cohort[m].value_counts(normalize=True)
        row = {"cluster": cl, "directions": ", ".join(links[cl]), "n_cells": int(m.sum()),
               "n_patients": int((cohort[m].value_counts() >= 1).sum()),
               "top_patient": pat.index[0], "top_patient_share": float(pat.iloc[0]),
               "patients_for_80pct": int((pat.cumsum() < 0.8).sum() + 1)}
        for col in ("treatment", "response", "phase"):
            for k, v in meta.loc[m, col].astype(str).value_counts(normalize=True).items():
                row[f"{col}_{k}"] = float(v)
        for col in ("n_genes_by_counts", "pct_counts_mt", "doublet_score", "cnv_score"):
            row[f"median_{col}"] = float(meta.loc[m, col].median())
            row[f"median_{col}_all"] = float(meta[col].median())
        comp.append(row)
    comp = pd.DataFrame(comp).set_index("cluster")
    comp.to_csv(tdir / f"state_composition_{key}_{drvi_run}.csv")
    print(f"[table] {tdir / f'state_composition_{key}_{drvi_run}.csv'}")

    # ---------------------------------------------------------------- 4. link to the DRs
    hvg_genes = list(hvg.var_names)
    hv = set(hvg_genes)
    n_top = args.n_top_genes
    dec = pd.read_csv(C.top_genes_tsv(n_top), sep="\t", index_col=0)
    emb_var = ad.read_h5ad(CS.tum_dir() / f"embed_{drvi_run}.h5ad", backed="r").var
    live = set(emb_var.loc[~emb_var["vanished_positive_direction"].astype(bool), "title"] + "+") | \
        set(emb_var.loc[~emb_var["vanished_negative_direction"].astype(bool), "title"] + "-")
    dec = dec[[c for c in dec.columns if c in live]]
    lists = {cl: [g for g in up[up["cluster"] == cl]["gene"] if g in hv][:n_top] for cl in picks}
    ov = []
    for cl, genes in lists.items():
        gs = set(genes)
        for d in dec.columns:
            dg = set(dec[d].dropna())
            k = len(gs & dg)
            p = float(hypergeom.sf(k - 1, len(hv), len(dg), len(gs))) if k else 1.0
            ov.append({"cluster": cl, "direction": d, "n_markers_hvg": len(gs), "overlap": k,
                       "p": p, "linked": d in links[cl]})
    ov = pd.DataFrame(ov)
    ov["fdr_bh"] = multipletests(ov["p"], method="fdr_bh")[1]
    ov["rank_in_cluster"] = ov.groupby("cluster")["p"].rank(method="min")
    ov.to_csv(tdir / f"state_marker_vs_decoder_{key}_{drvi_run}.csv", index=False)
    print(f"[table] {tdir / f'state_marker_vs_decoder_{key}_{drvi_run}.csv'}")

    # ---------------------------------------------------------------- 5. ORA
    libs = {c.name: C.read_gmt(C.gmt_path(c))
            for c in (SC.get("scie"), SC.get("emt"), SC.resolve(argparse.Namespace(
                collection="gavish", all_metaprograms=False)))}
    hm = CS.tum_dir() / "msigdb_hallmark_2020.json"
    if hm.exists():
        libs["hallmark"] = json.loads(hm.read_text())
    ora = []
    for lib, sets in libs.items():
        sets_bg = {k: [g for g in v if g in hv] for k, v in sets.items()}
        sets_bg = {k: v for k, v in sets_bg.items() if len(v) >= C.MIN_SIGNATURE_GENES}
        long = run_ora({cl: g for cl, g in lists.items() if g}, sets_bg, hvg_genes, lib, key)
        long, _ = apply_bh(long, sum(1 for g in lists.values() if g), len(sets_bg))
        long["library"] = lib
        ora.append(long)
    ora = pd.concat(ora, ignore_index=True)
    ora = ora.drop(columns=["Genes"], errors="ignore")
    ora.to_csv(tdir / f"state_ora_{key}_{drvi_run}.csv", index=False)
    print(f"[table] {tdir / f'state_ora_{key}_{drvi_run}.csv'}")

    # ---------------------------------------------------------------- report
    C.banner("per cluster")
    for cl in picks:
        c = comp.loc[cl]
        mk = marker_tbl[marker_tbl["cluster"] == cl]
        o = ov[ov["cluster"] == cl].sort_values("p")
        link = o[o["linked"]]
        sig_ora = ora[(ora["gene_list"] == cl) & ora["significant"]].sort_values("fdr_bh")
        print(f"\n{unit} {cl}  ({c['n_cells']:,} cells, {c['directions']})")
        print(f"  patients: {c['n_patients']}, top {c['top_patient']} "
              f"{c['top_patient_share']:.0%}, {c['patients_for_80pct']} make 80%")
        tr = ", ".join(f"{a} {c.get(f'treatment_{a}', 0):.0%}" for a in ("BASE", "PD1", "RTPD1"))
        print(f"  treatment: {tr}")
        print(f"  median genes {c['median_n_genes_by_counts']:.0f} (all "
              f"{c['median_n_genes_by_counts_all']:.0f}), mt {c['median_pct_counts_mt']:.1f}% "
              f"(all {c['median_pct_counts_mt_all']:.1f}%), doublet "
              f"{c['median_doublet_score']:.2f} (all {c['median_doublet_score_all']:.2f})")
        print(f"  {len(mk)} up markers; top 15: " + ", ".join(
            f"{r.gene}({r.logfoldchanges:.1f}"
            + (f";{int(r.n_patients_up)}/{int(r.n_patients_tested)}p" if r.n_patients_tested == r.n_patients_tested else "")
            + ")" for r in mk.head(15).itertuples()))
        print("  decoder overlap, linked: " + "; ".join(
            f"{r.direction} {r.overlap}/{r.n_markers_hvg} p={r.p:.1e} (rank {int(r.rank_in_cluster)})"
            for r in link.itertuples()))
        print("  decoder overlap, best 3: " + "; ".join(
            f"{r.direction} {r.overlap} p={r.p:.1e}" for r in o.head(3).itertuples()))
        print("  ORA (FDR<0.05): " + ("; ".join(
            f"{r.library}:{r.Term} {r.fdr_bh:.1e}" for r in sig_ora.head(6).itertuples())
            or "none"))

    # ---------------------------------------------------------------- figures
    # (a) dot plot, top 6 per picked cluster, all clusters of the resolution shown as columns
    #     would drown the picked ones: the others are pooled as 'other'
    #     With the max(DR) groups (~60 picks) the depth drops to 2 markers per group.
    n_dot = 6 if len(picks) <= 12 else 2
    top6 = {cl: [g for g in marker_tbl[marker_tbl["cluster"] == cl]["gene"].head(n_dot)]
            for cl in picks}
    grp = lab.where(lab.isin(picks), "other")
    full.obs["_pick"] = pd.Categorical(grp.values, categories=picks + ["other"])
    dp = sc.pl.dotplot(full, var_names={f"{k}" if argmax else f"cl {k}": v
                                        for k, v in top6.items() if v},
                       groupby="_pick", standard_scale="var", show=False, return_fig=True,
                       title=f"top markers (one-vs-rest, all genes) of the picked {unit}s, {key}")
    dp.savefig(ddir / f"dotplot_markers_{key}_{drvi_run}.png", dpi=300, bbox_inches="tight")
    plt.close("all")
    print(f"[fig] {ddir / f'dotplot_markers_{key}_{drvi_run}.png'}")

    # (b) composition: patient and treatment shares per picked cluster, all cells for scale
    fig, axes = plt.subplots(1, 2, figsize=(13, (0.45 if len(picks) <= 12 else 0.26)
                                            * (len(picks) + 1) + 1.8),
                             gridspec_kw={"width_ratios": [3, 1.2]})
    rows_ = picks + ["all"]
    pats = cohort.value_counts().index.tolist()
    pal = dict(zip(pats, plt.get_cmap("tab20")(np.linspace(0, 1, len(pats)))))
    for ax, col, order, colours in (
            (axes[0], cohort, pats, pal),
            (axes[1], meta["treatment"].astype(str), ["BASE", "PD1", "RTPD1"],
             {"BASE": "#2a78d6", "PD1": "#eb6834", "RTPD1": "#1baf7a"})):
        for i, cl in enumerate(rows_):
            m = np.ones(len(lab), bool) if cl == "all" else (lab == cl).values
            share = col[m].value_counts(normalize=True).reindex(order).fillna(0)
            left = 0.0
            for k, v in share.items():
                ax.barh(i, v, left=left, color=colours[k], height=0.7,
                        edgecolor="white", lw=0.5, label=k if i == 0 else None)
                left += v
        ax.set_yticks(range(len(rows_)), [f"{unit} {c}" if c != "all" else "all cells"
                                          for c in rows_], fontsize=8 if len(rows_) <= 13 else 6)
        ax.invert_yaxis()
        ax.set_xlim(0, 1)
        ax.set_xlabel("share of cells", fontsize=8)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
    axes[0].set_title("patient", fontsize=9)
    axes[1].set_title("treatment", fontsize=9)
    axes[0].legend(fontsize=6, ncol=1, frameon=False, bbox_to_anchor=(1.0, 1.0), loc="upper left")
    axes[1].legend(fontsize=7, frameon=False, bbox_to_anchor=(1.0, 1.0), loc="upper left")
    fig.suptitle(f"Composition of the picked {unit}s, {key}", fontsize=10)
    fig.tight_layout()
    fig.savefig(ddir / f"composition_{key}_{drvi_run}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {ddir / f'composition_{key}_{drvi_run}.png'}")

    # (c) marker vs decoder overlap: picked clusters x the picked directions, -log10 p
    dirs = list(dict.fromkeys(d for cl in picks for d in links[cl]))
    mat = (ov[ov["direction"].isin(dirs)].pivot(index="cluster", columns="direction", values="p")
           .reindex(index=picks, columns=dirs))
    cnt = (ov[ov["direction"].isin(dirs)].pivot(index="cluster", columns="direction",
                                                values="overlap").reindex(index=picks, columns=dirs))
    val = -np.log10(mat.clip(lower=1e-300))
    big = len(dirs) > 12       # the max(DR) groups: a ~60 x 60 matrix, its diagonal boxed
    cw, fs = (0.26, 4) if big else (0.8, 7)
    fig, ax = plt.subplots(figsize=(cw * len(dirs) + 2.5, (0.26 if big else 0.55) * len(picks)
                                    + 1.8))
    im = ax.imshow(val.values, cmap="Blues", aspect="auto", vmin=0,
                   vmax=max(10.0, float(np.nanpercentile(val.values, 95))))
    for (i, j), v in np.ndenumerate(val.values):
        ax.text(j, i, f"{int(cnt.values[i, j])}", ha="center", va="center", fontsize=fs,
                color="white" if v > 0.6 * im.get_clim()[1] else INK)
        if dirs[j] in links[picks[i]]:
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="k",
                                       lw=1.0 if big else 1.5))
    ax.set_xticks(range(len(dirs)), dirs, rotation=45 if not big else 90,
                  ha="right" if not big else "center", fontsize=fs + 1)
    ax.set_yticks(range(len(picks)), [f"{unit} {c}" for c in picks], fontsize=fs + 1)
    fig.colorbar(im, ax=ax, fraction=0.04, label="-log10 p, hypergeometric (HVG background)")
    pair = "the group's own direction" if argmax else "the SMI pair"
    ax.set_title(f"{unit.capitalize()} markers (HVG, top {n_top}) vs decoder top {n_top} genes "
                 f"of the linked directions\nnumbers = shared genes; boxed = {pair}", fontsize=9)
    fig.tight_layout()
    fig.savefig(ddir / f"markers_vs_decoder_{key}_{drvi_run}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {ddir / f'markers_vs_decoder_{key}_{drvi_run}.png'}")


def _auroc(values: np.ndarray, inside: np.ndarray) -> float:
    """AUROC of `values` separating `inside` from the rest, via ranks (ties averaged)."""
    from scipy.stats import rankdata
    n1, n0 = int(inside.sum()), int((~inside).sum())
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(values)
    return float((r[inside].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def direction_interpretation(drvi_run: str, live: list[str]) -> pd.DataFrame:
    """One row per live direction: everything the phase already says about it.

    Nothing is recomputed. Read, in order of how much weight a reader should give them:
      name         05_8's naming bar (both routes, same family, rho >= 0.30), consensus table
      depth        the DIMENSION's |rho| with log UMI >= 0.30 (dr_naming's DIM_DEPTH_FLAG)
      route_b      Route B hits (FDR < 0.05) per collection, whatever Route A said
      hallmark     the best Hallmark term of 05_7 on this direction
      top_genes    the decoder's top 10 genes
    and `interpretation` is the first of those that exists, spelled out.
    """
    run = drvi_run
    pr = C.OUT_ROOT / "tables_pruned"
    naming = pd.read_csv(pr / "consensus" / run / f"dr_naming_consensus_{run}.csv",
                         comment="#", index_col=0)
    progs = pd.read_csv(pr / "consensus" / run / f"programme_dimensions_consensus_{run}.csv",
                        comment="#")
    named = {}
    for r in progs.itertuples():
        for d in str(r.dimensions).split(","):
            named.setdefault(d.strip(), []).append(r.programme)
    route_b = {}
    for coll in ("scie", "emt", "gavish_tnbc"):
        cv = pd.read_csv(pr / coll / run / f"convergence_{coll}_{run}.csv", comment="#",
                         index_col=0)
        for d, r in cv[cv["B_significant"]].iterrows():
            route_b.setdefault(d, []).append((float(r["B_fdr"]),
                                              f"{r['B_best_signature']} (FDR {r['B_fdr']:.0e}, "
                                              f"A rho {r['A_rho']:.2f})"))
    hm = pd.read_csv(pr / "scie" / run / f"factor_first_hallmark_significant_scie_{run}.csv",
                     comment="#")
    hm_best = (hm.sort_values("fdr_bh_global").groupby("dimension").head(2)
               .groupby("dimension")["Term"].agg(lambda s: "; ".join(s)))
    dec = pd.read_csv(C.top_genes_tsv(N_TOP_GENES), sep="\t", index_col=0)

    rows = []
    for d in live:
        dim = d[:-1].strip()
        cov = naming.loc[dim] if dim in naming.index else None
        depth = cov is not None and abs(float(cov["rho_log_umi"])) >= 0.30
        nm = named.get(d, [])
        rb = [t for _, t in sorted(route_b.get(d, []))]     # strongest Route B hit first
        hmk = hm_best.get(d, "")
        genes = ", ".join(dec[d].head(10)) if d in dec.columns else ""
        if nm:
            interp = "named: " + ", ".join(nm)
        elif depth:
            interp = f"depth axis (dimension rho log UMI {float(cov['rho_log_umi']):+.2f})"
        elif rb:
            interp = "Route B only: " + rb[0].split(" (")[0]
        elif hmk:
            interp = "Hallmark only: " + hmk.split(";")[0]
        else:
            interp = "unnamed - read the genes"
        rows.append({
            "direction": d, "interpretation": interp,
            "named_programmes": "; ".join(nm), "route_b_hits": "; ".join(rb),
            "hallmark_top": hmk, "top10_decoder_genes": genes,
            "dim_rho_log_umi": float(cov["rho_log_umi"]) if cov is not None else np.nan,
            "dim_rho_S": float(cov["rho_S_score"]) if cov is not None else np.nan,
            "dim_rho_G2M": float(cov["rho_G2M_score"]) if cov is not None else np.nan,
            "dim_rho_pct_mt": float(cov["rho_pct_mt"]) if cov is not None else np.nan,
        })
    return pd.DataFrame(rows).set_index("direction")


def drvi_states_dr(args, drvi_run: str, argmax: bool = False):
    """DR enrichment per cluster, for EVERY cluster of one state resolution, and what it means.

    ENRICHMENT is an AUROC: how well a direction's value separates the cluster's cells from
    all the others (0.5 = not at all, 1 = perfectly). Computed on the signed coordinate - the
    `+` direction is the coordinate, the `-` direction its negative - so each dimension gives
    one AUROC per side and they sum to 1. Vanished directions are dropped, as in the SMI.
    Unlike the SMI it has a direction of effect and a scale everyone reads the same way.

    WITHIN PATIENT, for the top STATE_DR_TOP directions of each cluster: the same AUROC inside
    every patient with >= STATE_DE_MIN_CELLS cells both in and out of the cluster, then the
    median. A cluster made of one patient has one such patient at most, and the column says so
    rather than hiding it: an enrichment that holds in one patient describes that patient.

    Descriptive and circular, like everything on these clusters: they were cut from these
    directions, so enrichment is expected. What it adds is the READING - which programme,
    confounder or patient each cluster is made of.

    With `argmax` the rows are the tested max(DR) groups. There the circularity is total for
    the group's own direction (rank 1 by construction); the information is in ranks 2-3, the
    directions that come up WITH it, and in the within-patient column.
    """
    P = state_partition(args, drvi_run, argmax)
    key, unit = P["key"], P["unit"]
    _, fdir, tdir = states_paths(drvi_run)
    odir = fdir / f"state_dr_enrichment_{key}"
    odir.mkdir(exist_ok=True)
    de_path = tdir / f"state_de_{key}_{drvi_run}.csv.gz"
    if not de_path.exists():
        sys.exit(f"missing {de_path}: run --states-de first")
    de = pd.read_csv(de_path)
    de["cluster"] = de["cluster"].astype(str)
    C.banner(f"DRVI states: DR enrichment and interpretation, every {unit} of {key}")

    embed = ad.read_h5ad(CS.tum_dir() / f"embed_{drvi_run}.h5ad")
    embed.obs_names = embed.obs_names.astype(str)
    meta = embed.obs
    lab = P["lab"].reindex(embed.obs_names).astype(str)
    clusters = list(P["tested"])
    X = np.asarray(embed.X)
    titles = embed.var["title"].astype(str).to_numpy()
    cols, vals = [], []
    for j, t in enumerate(titles):
        if not bool(embed.var["vanished_positive_direction"].iloc[j]):
            cols.append(f"{t}+"); vals.append(X[:, j])
        if not bool(embed.var["vanished_negative_direction"].iloc[j]):
            cols.append(f"{t}-"); vals.append(-X[:, j])
    V = np.column_stack(vals)
    print(f"    {len(clusters)} {unit}s x {len(cols)} live directions", flush=True)

    # ---------------------------------------------------------------- AUROC, pooled
    from scipy.stats import rankdata
    ranks = np.column_stack([rankdata(V[:, k]) for k in range(V.shape[1])])
    au = pd.DataFrame(index=clusters, columns=cols, dtype=float)
    smd = pd.DataFrame(index=clusters, columns=cols, dtype=float)
    sd = V.std(0)
    for cl in clusters:
        m = (lab == cl).values
        n1, n0 = m.sum(), (~m).sum()
        au.loc[cl] = (ranks[m].sum(0) - n1 * (n1 + 1) / 2) / (n1 * n0)
        smd.loc[cl] = (V[m].mean(0) - V[~m].mean(0)) / np.where(sd > 0, sd, 1)
    smi_path = tdir / f"smi_matrix_{key}_{drvi_run}.csv"
    smi = pd.read_csv(smi_path, index_col=0) if smi_path.exists() else None

    interp = direction_interpretation(drvi_run, cols)
    interp.to_csv(tdir / f"state_dr_interpretation_{key}_{drvi_run}.csv")
    print(f"[table] {tdir / f'state_dr_interpretation_{key}_{drvi_run}.csv'}")

    # ---------------------------------------------------------------- long table
    cohort = meta["cohort"].astype(str).values
    long = []
    for cl in clusters:
        m = (lab == cl).values
        order = au.loc[cl].sort_values(ascending=False)
        for rank, (d, a) in enumerate(order.items(), start=1):
            row = {"cluster": cl, "direction": d, "rank": rank, "auroc": float(a),
                   "smd": float(smd.loc[cl, d]),
                   "smi": float(smi.loc[d, cl]) if smi is not None and d in smi.index
                   and cl in smi.columns else np.nan,
                   "enriched": bool(a >= STATE_DR_AUROC)}
            if rank <= STATE_DR_TOP:
                k = cols.index(d)
                within = []
                for pat in np.unique(cohort):
                    pm = cohort == pat
                    if (pm & m).sum() >= STATE_DE_MIN_CELLS and (pm & ~m).sum() >= STATE_DE_MIN_CELLS:
                        within.append(_auroc(V[pm, k], m[pm]))
                row["n_patients_within"] = int(len(within))
                row["median_auroc_within"] = float(np.median(within)) if within else np.nan
            long.append(row)
    long = pd.DataFrame(long).merge(interp[["interpretation"]], left_on="direction",
                                    right_index=True, how="left")
    long.to_csv(tdir / f"state_dr_enrichment_{key}_{drvi_run}.csv", index=False)
    print(f"[table] {tdir / f'state_dr_enrichment_{key}_{drvi_run}.csv'}")

    # ---------------------------------------------------------------- one row per cluster
    up = de[(de["logfoldchanges"] > STATE_DE_MIN_LFC) & (de["pct_nz_group"] >= STATE_DE_MIN_PCT)]
    up = up.sort_values(["cluster", "scores"], ascending=[True, False])
    summ = []
    for cl in clusters:
        m = (lab == cl).values
        pat = pd.Series(cohort[m]).value_counts(normalize=True)
        tr = meta.loc[m, "treatment"].astype(str).value_counts(normalize=True)
        top = long[(long["cluster"] == cl) & (long["rank"] <= STATE_DR_TOP)]
        summ.append({
            "cluster": cl, "n_cells": int(m.sum()),
            "n_patients_ge20": int((pd.Series(cohort[m]).value_counts() >= 20).sum()),
            "top_patient": pat.index[0], "top_patient_share": float(pat.iloc[0]),
            "BASE": float(tr.get("BASE", 0)), "PD1": float(tr.get("PD1", 0)),
            "RTPD1": float(tr.get("RTPD1", 0)),
            "median_genes": float(meta.loc[m, "n_genes_by_counts"].median()),
            "median_pct_mt": float(meta.loc[m, "pct_counts_mt"].median()),
            "median_doublet": float(meta.loc[m, "doublet_score"].median()),
            "n_enriched_directions": int(((long["cluster"] == cl) & long["enriched"]).sum()),
            "top_directions": "; ".join(
                f"{r.direction} (AUROC {r.auroc:.2f}, within {r.median_auroc_within:.2f} "
                f"in {int(r.n_patients_within)}p)" for r in top.itertuples()),
            "top_direction_reading": "; ".join(f"{r.direction}: {r.interpretation}"
                                               for r in top.itertuples()),
            "top_markers": ", ".join(up[up["cluster"] == cl]["gene"].head(10)),
        })
    summ = pd.DataFrame(summ).set_index("cluster")
    summ.to_csv(tdir / f"state_cluster_summary_{key}_{drvi_run}.csv")
    print(f"[table] {tdir / f'state_cluster_summary_{key}_{drvi_run}.csv'}")

    C.banner("per cluster")
    med_g = float(meta["n_genes_by_counts"].median())
    for cl, r in summ.iterrows():
        print(f"\n{unit} {cl}: {r.n_cells:,} cells; {r.top_patient} {r.top_patient_share:.0%} "
              f"({r.n_patients_ge20} patients with >= 20 cells); BASE {r.BASE:.0%} PD1 "
              f"{r.PD1:.0%} RTPD1 {r.RTPD1:.0%}; genes {r.median_genes:.0f} (all {med_g:.0f}), "
              f"mt {r.median_pct_mt:.1f}%")
        print(f"  enriched directions (AUROC >= {STATE_DR_AUROC}): {r.n_enriched_directions}")
        for t in long[(long["cluster"] == cl) & (long["rank"] <= STATE_DR_TOP)].itertuples():
            print(f"    {t.direction:7s} AUROC {t.auroc:.2f} (within {t.median_auroc_within:.2f}, "
                  f"{int(t.n_patients_within)}p)  {t.interpretation}")
        print(f"  markers: {r.top_markers}")

    # ---------------------------------------------------------------- figures (new names)
    # (a) clusters x the directions enriched somewhere, AUROC
    keep = [d for d in cols if (au[d] >= STATE_DR_AUROC).any()]
    pos = {c: i for i, c in enumerate(clusters)}
    best_cl = au[keep].idxmax(axis=0).map(pos)
    keep = sorted(keep, key=lambda d: (best_cl[d], -au[d].max()))
    mat = au[keep]
    short = {d: interp.loc[d, "interpretation"].replace("named: ", "").replace(
        "Route B only: ", "B: ").replace("Hallmark only: ", "HM: ")[:34] for d in keep}
    fig, ax = plt.subplots(figsize=(0.32 * len(keep) + 3.5, 0.34 * len(clusters) + 3.4))
    im = ax.imshow(mat.values, cmap="Blues", vmin=0.5, vmax=1.0, aspect="auto",
                   interpolation="nearest")
    for (i, j), v in np.ndenumerate(mat.values):
        if v >= STATE_DR_AUROC:
            ax.text(j, i, f"{v:.2f}".lstrip("0"), ha="center", va="center", fontsize=5,
                    color="white" if v > 0.9 else INK)
    ax.set_yticks(range(len(clusters)), [
        f"{c}  ({summ.loc[c, 'n_cells']:,}, {summ.loc[c, 'top_patient'].replace('Patient', 'P')} "
        f"{summ.loc[c, 'top_patient_share']:.0%})" for c in clusters], fontsize=7)
    ax.set_xticks(range(len(keep)), [f"{d} · {short[d]}" for d in keep], rotation=90,
                  fontsize=6)
    ax.set_ylabel(f"{unit} ({key})  (cells, top patient share)", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01, label="AUROC, cluster vs rest")
    ax.set_title(f"DR enrichment per {unit}, {drvi_run}, {P['desc']}: the {len(keep)} "
                 f"live directions with AUROC >= {STATE_DR_AUROC} in at least one cluster\n"
                 "labels: the phase's reading of each direction (name / depth / B = Route B "
                 "only / HM = Hallmark only); circular by construction", fontsize=9)
    fig.tight_layout()
    fig.savefig(odir / f"dr_enrichment_heatmap_{key}_{drvi_run}.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {odir / f'dr_enrichment_heatmap_{key}_{drvi_run}.png'}")

    # (b) composition of every cluster: patient and treatment
    fig, axes = plt.subplots(1, 2, figsize=(13, 0.34 * len(clusters) + 2.2),
                             gridspec_kw={"width_ratios": [3, 1.2]})
    pats = pd.Series(cohort).value_counts().index.tolist()
    pal = dict(zip(pats, plt.get_cmap("tab20")(np.linspace(0, 1, len(pats)))))
    rows_ = clusters + ["all"]
    for ax, col, order, colours in (
            (axes[0], pd.Series(cohort), pats, pal),
            (axes[1], meta["treatment"].astype(str).reset_index(drop=True),
             ["BASE", "PD1", "RTPD1"], {"BASE": "#2a78d6", "PD1": "#eb6834", "RTPD1": "#1baf7a"})):
        for i, cl in enumerate(rows_):
            m = np.ones(len(lab), bool) if cl == "all" else (lab == cl).values
            share = col[m].value_counts(normalize=True).reindex(order).fillna(0)
            left = 0.0
            for k, v in share.items():
                ax.barh(i, v, left=left, color=colours[k], height=0.7, edgecolor="white",
                        lw=0.4, label=k if i == 0 else None)
                left += v
        ax.set_yticks(range(len(rows_)), [f"{unit} {c}" if c != "all" else "all cells"
                                          for c in rows_], fontsize=7)
        ax.invert_yaxis()
        ax.set_xlim(0, 1)
        ax.set_xlabel("share of cells", fontsize=8)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
    axes[0].set_title("patient", fontsize=9)
    axes[1].set_title("treatment", fontsize=9)
    axes[0].legend(fontsize=6, frameon=False, bbox_to_anchor=(1.0, 1.0), loc="upper left")
    axes[1].legend(fontsize=7, frameon=False, bbox_to_anchor=(1.0, 1.0), loc="upper left")
    fig.suptitle(f"Composition of every {unit}, {key}", fontsize=10)
    fig.tight_layout()
    fig.savefig(odir / f"composition_all_clusters_{key}_{drvi_run}.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {odir / f'composition_all_clusters_{key}_{drvi_run}.png'}")

    # Only the partitions with a curated reading: the Leiden clusters and the plain groups.
    if not argmax:
        fig_cluster_factors(summ, long, key, drvi_run, odir)
    elif key == "drvi_argmax_plain":
        fig_cluster_factors(summ, long, key, drvi_run, odir, reading=ARGMAX_PLAIN_READING,
                            unit="group", what="plain max(DR) group")


def fig_cluster_factors(summ: pd.DataFrame, long: pd.DataFrame, key: str, drvi_run: str,
                        odir, reading: dict | None = None, unit: str = "cluster",
                        what: str = "DRVI-Leiden cluster") -> None:
    """One bar per cluster: its size, coloured by the factor that explains it, with its DRs.

    Bars are grouped by factor (programme first) and sorted by size inside a group. The bar
    length is the cell count on a log axis - the clusters span 13 to 6,844 cells - and the
    label carries the top enriched directions (pooled AUROC, and the within-patient median
    where there is more than one patient to take it over), the curated reading, and the
    patient make-up that decides most of these calls.
    """
    reading = CLUSTER_READING if reading is None else reading
    order = list(FACTOR_STYLE)
    df = summ.copy()
    df["factor"] = [reading.get(c, ("none", ""))[0] for c in df.index]
    df["reading"] = [reading.get(c, ("none", "not curated"))[1] for c in df.index]
    df["rank"] = df["factor"].map({f: i for i, f in enumerate(order)})
    df = df.sort_values(["rank", "n_cells"], ascending=[True, False])
    tdir = states_paths(drvi_run)[2]
    df.drop(columns="rank").to_csv(tdir / f"state_cluster_factors_{key}_{drvi_run}.csv")
    print(f"[table] {tdir / f'state_cluster_factors_{key}_{drvi_run}.csv'}")

    def dr_text(cl):
        top = long[(long["cluster"] == cl) & (long["rank"] <= 2)]
        parts = []
        for r in top.itertuples():
            w = (f", {r.median_auroc_within:.2f} within" if r.n_patients_within
                 and r.n_patients_within > 1 else "")
            parts.append(f"{r.direction} {r.auroc:.2f}{w}")
        return "; ".join(parts)

    n = len(df)
    y = np.arange(n)[::-1]
    fig, ax = plt.subplots(figsize=(13.5, 0.42 * n + 2.4))
    for yi, (cl, r) in zip(y, df.iterrows()):
        colour = FACTOR_STYLE[r["factor"]][0]
        ax.barh(yi, r["n_cells"], height=0.6, color=colour)
        pat = (f"{r['top_patient'].replace('Patient', 'P')} {r['top_patient_share']:.0%}, "
               f"{int(r['n_patients_ge20'])} pts")
        ax.text(r["n_cells"] * 1.12, yi,
                f"{dr_text(cl)}   |   {r['reading']}   |   {pat}",
                va="center", ha="left", fontsize=7, color="0.25")
    ax.set_xscale("log")
    ax.set_xlim(8, df["n_cells"].max() * 400)
    # The right half of the axis is room for the labels, not data: ticks stop at the data.
    ticks = [t for t in (10, 100, 1000, 10000, 100000) if t <= df["n_cells"].max() * 2]
    ax.set_xticks(ticks, [f"{t:,}" for t in ticks])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.spines["bottom"].set_bounds(8, ticks[-1])
    ax.set_yticks(y, [f"{unit} {c}" for c in df.index], fontsize=8)
    ax.set_xlabel(f"cells in the {unit} (log scale)", fontsize=9)
    for s_ in ("top", "right", "left"):
        ax.spines[s_].set_visible(False)
    ax.spines["bottom"].set_color("0.7")
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="0.9", lw=0.6)
    ax.set_axisbelow(True)
    handles = [plt.Rectangle((0, 0), 1, 1, color=FACTOR_STYLE[f][0]) for f in order
               if f in set(df["factor"])]
    labels = [FACTOR_STYLE[f][1] for f in order if f in set(df["factor"])]
    ax.legend(handles, labels, fontsize=7.5, frameon=False, loc="lower right")
    ax.set_title(f"What each {what} is ({key}, {drvi_run})", fontsize=12,
                 loc="left", pad=40)
    ax.text(0, 1.005, f"colour = the one factor that best explains the {unit} (curated from "
            "composition, QC, DE markers and DR enrichment);\nlabel = top two DR directions "
            f"(AUROC {unit} vs rest, and median within patient) | reading | top patient, "
            "patients with >= 20 cells",
            transform=ax.transAxes, fontsize=8, color="0.35", va="bottom", linespacing=1.4)
    fig.tight_layout()
    fig.savefig(odir / f"cluster_factors_barplot_{key}_{drvi_run}.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {odir / f'cluster_factors_barplot_{key}_{drvi_run}.png'}")


# --------------------------------------------------------------------------- #
# The ORA, identical for every arm
# --------------------------------------------------------------------------- #

def run_ora(lists: dict[str, list[str]], sets: dict[str, list[str]], background: list[str],
            arm: str, resolution) -> pd.DataFrame:
    """Every gene list against every signature: offline hypergeometric, declared background.

    `gp.enrich`, as 05_7 - never Enrichr's implicit all-human-genes universe. Benjamini-
    Hochberg is NOT applied here: it is applied once per arm-and-resolution by `apply_bh`,
    across every pair actually tested including the ones with no overlap, which is the only
    denominator that makes two arms of different widths comparable.
    """
    records = []
    for name, genes in lists.items():
        try:
            res = gp.enrich(gene_list=genes, gene_sets=sets, background=background,
                            outdir=None)
        except Exception as exc:                       # noqa: BLE001 - reported, not raised
            print(f"      {name}: ORA FAILED ({exc})", flush=True)
            continue
        res_df = getattr(res, "results", None)
        if res_df is None or len(res_df) == 0:
            continue                                   # no overlap at all: p = 1, added below
        df = pd.DataFrame(res_df).copy()
        df.insert(0, "gene_list", str(name))
        records.append(df)
    if not records:
        return pd.DataFrame(columns=["gene_list", "Term", "P-value", "arm", "resolution"])
    long = pd.concat(records, ignore_index=True)
    long["arm"] = arm
    long["resolution"] = resolution
    return long


def apply_bh(long: pd.DataFrame, n_lists: int, n_sets: int) -> tuple[pd.DataFrame, int]:
    """BH across all `n_lists x n_sets` pairs of ONE arm at ONE resolution.

    The pairs gseapy left out are the ones with no overlap. They are p = 1 and can never
    become significant, but they belong in the denominator, so they are added back - exactly
    as 05_7 does it. Corrected WITHIN an arm, because that is the correction each pipeline
    would actually apply to itself: an arm is not penalised for the other one's width, and
    an arm that buys its hits with more tests pays for them in its own denominator.
    """
    n_pairs = n_lists * n_sets
    if long.empty:
        return long.assign(fdr_bh=pd.Series(dtype=float),
                           significant=pd.Series(dtype=bool)), n_pairs
    pvals = np.concatenate([long["P-value"].values, np.ones(max(0, n_pairs - len(long)))])
    _, padj, _, _ = multipletests(pvals, alpha=FDR, method="fdr_bh")
    out = long.copy()
    out["fdr_bh"] = padj[:len(long)]
    out["significant"] = out["fdr_bh"] < FDR
    return out, n_pairs


# --------------------------------------------------------------------------- #
# The diagnostic that explains the gap
# --------------------------------------------------------------------------- #

def eta_squared(y: np.ndarray, codes: np.ndarray, n_groups: int) -> float:
    """Share of a per-cell score's variance that lies BETWEEN clusters. One-way eta^2.

    This is the number the whole argument rests on, and it is computed on 05_6's cached
    within-cohort z scores - the same values 05_9 measured presence on, never re-scored
    here. A programme with a high eta^2 is a programme the partition has a group for; a
    programme with a low one is a gradient ACROSS the groups, and a gradient is what a
    one-vs-rest DE has nothing to test. If the signatures the cluster arm misses are the
    low-eta^2 ones, the miss is explained rather than merely counted.
    """
    y = np.asarray(y, dtype=np.float64)
    grand = y.mean()
    ss_tot = float(((y - grand) ** 2).sum())
    if ss_tot <= 0:
        return float("nan")
    counts = np.bincount(codes, minlength=n_groups).astype(np.float64)
    sums = np.bincount(codes, weights=y, minlength=n_groups)
    ok = counts > 0
    means = np.zeros_like(counts)
    means[ok] = sums[ok] / counts[ok]
    ss_between = float((counts[ok] * (means[ok] - grand) ** 2).sum())
    return ss_between / ss_tot


def cluster_profile(labels: pd.Series, cohorts: pd.Series) -> pd.DataFrame:
    """Per cluster: size, and how much of it is one patient.

    `cohort_dominance` is here because the standard pipeline's other well-known failure is
    not about gradients at all: a cluster that is 90% one cohort makes its DE a patient
    contrast, and any programme it enriches for is that patient's. It is reported rather
    than filtered - filtering it would be deciding for the pipeline which of its clusters
    were allowed to count.
    """
    df = pd.DataFrame({"cluster": labels.values, "cohort": cohorts.reindex(labels.index).values})
    rows = []
    for cl, grp in df.groupby("cluster", observed=True):
        frac = grp["cohort"].value_counts(normalize=True)
        rows.append({"cluster": str(cl), "n_cells": len(grp),
                     "n_cohorts": int(grp["cohort"].nunique()),
                     "cohort_dominance": float(frac.iloc[0]),
                     "top_cohort": str(frac.index[0])})
    return pd.DataFrame(rows).sort_values("n_cells", ascending=False)


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #

def write_table(df: pd.DataFrame, name: str, coll, run_id: str, index: bool = False):
    """C.table_path for the path, this step's own header for the first lines.

    Not `C.write_table`, for 05_9's reason: that helper stamps "<Space> run <run_id>" from
    the EMBEDDINGS registry, which is right for a table measured in one space and wrong for
    one that compares pipelines across two.
    """
    path = C.table_path(name, coll, run_id)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"# {name} | collection {coll.name} ({coll.title}) | run {run_id} "
                 f"| {HEADER_NOTE} | {STEP}\n")
        for line in C.CAVEAT.split(". "):
            if line.strip():
                fh.write(f"# CAVEAT: {line.strip().rstrip('.')}.\n")
        df.to_csv(fh, index=index)
    print(f"[table] {path}  ({df.shape[0]} x {df.shape[1]})", flush=True)
    return path


def fig_recall_curve(summary: pd.DataFrame, coll, run_id: str, n_testable: int, args):
    """Recall against the resolution, and against the list budget that bought it.

    Two panels because the budget objection is real: the left one is the scan as it was run,
    the right one puts the same points on the axis the objection is about, so a reader can
    check whether the cluster arms are merely under-tested.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    leiden = summary[summary["unit"] == "cluster"]
    flat = summary[summary["unit"] == "axis"]

    for ax, xcol, xlabel in ((axes[0], "resolution", "Leiden resolution"),
                             (axes[1], "n_gene_lists", "gene lists handed to the ORA")):
        for arm, grp in leiden.groupby("arm", observed=True):
            grp = grp.sort_values(xcol)
            ax.plot(grp[xcol], grp["recall_named"], "-o", ms=4, lw=1.8,
                    color=ARM_COLOUR.get(arm, INK), label=ARM_LABEL.get(arm, arm))
            ax.plot(grp[xcol], grp["recall_any"], ":", lw=1.1, alpha=0.75,
                    color=ARM_COLOUR.get(arm, INK))
        for _, row in flat.iterrows():
            ax.axhline(row["recall_named"], color=ARM_COLOUR.get(row["arm"], INK), lw=1.8,
                       ls="--", label=ARM_LABEL.get(row["arm"], row["arm"]))
            ax.axhline(row["recall_any"], color=ARM_COLOUR.get(row["arm"], INK), lw=1.0,
                       ls=":", alpha=0.75)
            if xcol == "n_gene_lists":
                ax.plot([row["n_gene_lists"]], [row["recall_named"]], "D", ms=6,
                        color=ARM_COLOUR.get(row["arm"], INK))
        null = summary[summary["arm"] == "null_partition"]
        if len(null):
            ax.axhline(float(null["recall_named"].iloc[0]), color=ARM_COLOUR["null_partition"],
                       lw=1.4, ls="-.", label=ARM_LABEL["null_partition"])
        ax.set_xlabel(xlabel)
        ax.set_ylim(-0.03, 1.03)
        ax.grid(alpha=0.25, lw=0.6)
        if xcol == "n_gene_lists":
            ax.set_xscale("log")

    axes[0].set_ylabel(f"fraction of the {n_testable} testable signatures")
    handles, labels = axes[0].get_legend_handles_labels()
    seen, h, l = set(), [], []
    for hh, ll in zip(handles, labels):
        if ll not in seen:
            seen.add(ll); h.append(hh); l.append(ll)
    axes[0].legend(h, l, fontsize=7.5, loc="lower right", frameon=False)
    fig.suptitle(f"{coll.title}: what each pipeline recovers\n"
                 f"identical ORA downstream - top {args.n_top_genes} genes per list, "
                 f"hypergeometric against the {n_testable}-set collection on the 2,000-HVG "
                 f"background, BH within arm at FDR < {FDR}.\n"
                 "SOLID: named, i.e. some unit is identified by it.   "
                 "DOTTED: any significant hit - the saturated metric, drawn to be dismissed",
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    C.savefig("recall_by_resolution", STEP, coll, fig, run_id=run_id)
    plt.close(fig)


def fig_recovery_matrix(rec: pd.DataFrame, coll, run_id: str, arms: list[str], ref_res):
    """Signature x arm, coloured by -log10 FDR. The table the recall number summarises."""
    cols = [f"neglog_fdr__{a}" for a in arms if f"neglog_fdr__{a}" in rec.columns]
    if not cols:
        return
    order = coll.order(rec.index.tolist(), for_figure=True)
    order = order + [s for s in rec.index if s not in order]
    M = rec.loc[order, cols].astype(float)
    fig, ax = plt.subplots(figsize=(C.fig_span(len(cols), 1.5, 3.4),
                                    C.fig_span(len(order), 0.26, 2.2)))
    thr = -np.log10(FDR)
    im = ax.imshow(M.values, aspect="auto", cmap="Blues",
                   vmin=0, vmax=max(thr * 2, float(np.nanmax(M.values)) or thr * 2))
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.values[i, j]
            named = bool(rec.loc[order[i], f"named__{cols[j].split('__', 1)[1]}"]) \
                if f"named__{cols[j].split('__', 1)[1]}" in rec.columns else False
            if np.isfinite(v) and v >= thr:
                ax.text(j, i, "*" if named else ".", ha="center", va="center",
                        fontsize=13 if named else 15,
                        color="#08306b" if named else "#7f7f7f")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([ARM_LABEL.get(c.split("__", 1)[1], c) for c in cols],
                       rotation=25, ha="right", fontsize=8)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=7)
    ax.set_title(f"{coll.title}: which signature each pipeline reached\n"
                 f"Leiden arms at resolution {ref_res}. "
                 f"*  some unit is IDENTIFIED BY it (its best hit)    "
                 f".  significant but never the best hit of any unit", fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.4, label="-log10 FDR (best gene list)")
    C.savefig("recovery_matrix", STEP, coll, fig, run_id=run_id)
    plt.close(fig)


def fig_eta2(rec: pd.DataFrame, coll, run_id: str, arm: str, ref_res):
    """The explanation: what the cluster arm could not name, against how clustered it is.

    x is eta^2(cluster) of the signature's own per-cell score at the reference resolution -
    how much of it lies BETWEEN clusters rather than across them. One row per signature,
    sorted, because the identities are the point: a strip plot of seventeen metaprogram names
    at two y values is unreadable, and the first version of this figure was.

    A signature low on this axis that the arm could not name is a GRADIENT the partition had
    no group for. One high on it that the arm could not name is NOT explained by this figure
    and must not be claimed to be - MP3_CELL_CYCLE_HMG_RICH is the case on this run.
    """
    xcol, ycol = f"eta2_cluster__{arm}", f"named__{arm}"
    if xcol not in rec.columns or ycol not in rec.columns:
        return
    d = rec[[xcol, ycol]].dropna(subset=[xcol]).sort_values(xcol)
    if d.empty:
        return
    named = d[ycol].astype(bool).values
    y = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(8.0, C.fig_span(len(d), 0.30, 2.4)))
    ax.hlines(y, 0, d[xcol].values, color="0.85", lw=1.2, zorder=1)
    ax.scatter(d.loc[named, xcol], y[named], s=70, marker="o",
               color=ARM_COLOUR.get(arm, INK), edgecolor="white", linewidth=0.8, zorder=3,
               label=f"names a unit ({int(named.sum())})")
    ax.scatter(d.loc[~named, xcol], y[~named], s=80, marker="X", color="#c0392b",
               edgecolor="white", linewidth=0.8, zorder=3,
               label=f"never a unit's best hit ({int((~named).sum())})")
    med_hit = float(d.loc[named, xcol].median()) if named.any() else np.nan
    med_miss = float(d.loc[~named, xcol].median()) if (~named).any() else np.nan
    for v, colour, lab in ((med_hit, ARM_COLOUR.get(arm, INK), "median, named"),
                           (med_miss, "#c0392b", "median, not named")):
        if np.isfinite(v):
            ax.axvline(v, color=colour, lw=1.2, ls="--", alpha=0.7, zorder=2,
                       label=f"{lab}: {v:.3f}")
    ax.set_yticks(y)
    ax.set_yticklabels(d.index, fontsize=7.5)
    ax.set_ylim(-0.8, len(d) - 0.2)
    ax.set_xlim(left=0)
    ax.set_xlabel(r"$\eta^2$(cluster) of the per-cell score  "
                  "- share of its variance lying BETWEEN clusters")
    ax.grid(alpha=0.25, lw=0.6, axis="x")
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")
    ax.set_title(f"{coll.title}: {ARM_LABEL.get(arm, arm)} at resolution {ref_res}\n"
                 "a programme the partition has no group for is a programme a one-vs-rest "
                 "DE cannot propose", fontsize=9)
    fig.tight_layout()
    C.savefig(f"eta2_vs_recovery_{arm}", STEP, coll, fig, run_id=run_id)
    plt.close(fig)


def fig_hits_per_list(per_list: pd.DataFrame, coll, run_id: str, arms: list[str], ref_res):
    """How many signatures one gene list is significant for: the entanglement read.

    A unit that hits eight signatures at once has not resolved eight programmes, it has
    found the block they share. This is the same statement 05_9 makes with `best_dim` - H1
    being the best correlate of all five immunogenicity lists - counted on the discovery
    side instead of the correlation side.
    """
    d = per_list[per_list["arm"].isin(arms)]
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    order = [a for a in arms if a in set(d["arm"])]
    width = 0.8 / max(1, len(order))
    mx = int(d["n_signatures_hit"].max())
    bins = np.arange(0, mx + 2)
    for i, arm in enumerate(order):
        v = d.loc[d["arm"] == arm, "n_signatures_hit"].values
        counts, _ = np.histogram(v, bins=bins)
        ax.bar(bins[:-1] + i * width - 0.4 + width / 2, counts / max(1, counts.sum()),
               width=width, color=ARM_COLOUR.get(arm, INK),
               label=f"{ARM_LABEL.get(arm, arm)}  (n = {len(v)}, mean {v.mean():.2f})")
    ax.set_xticks(bins[:-1])
    ax.set_xlabel("signatures one gene list is significant for")
    ax.set_ylabel("fraction of that arm's gene lists")
    ax.legend(fontsize=7.5, frameon=False)
    ax.grid(alpha=0.25, lw=0.6, axis="y")
    ax.set_title(f"{coll.title}: how many programmes one unit of analysis carries\n"
                 f"Leiden arms at resolution {ref_res}", fontsize=9)
    fig.tight_layout()
    C.savefig("hits_per_gene_list", STEP, coll, fig, run_id=run_id)
    plt.close(fig)


# --------------------------------------------------------------------------- #

def main():
    args = parse_args()
    coll = SC.resolve(args)
    drvi_run = CS.run_id(method="drvi")
    harmony_run = CS.run_id(method="harmony")
    run_id = drvi_run          # names the CELLS and the gene sets, not a space; see 05_9

    CS.banner(f"{STEP} - what does a standard pipeline recover? ({coll.name})")
    print(f"collection : {coll.title}")
    print(f"question   : {coll.question}")
    print(f"arms       : {', '.join(args.arms)}")
    print(f"reference  : {drvi_run}   |   harmony: {harmony_run}")
    print(f"resolutions: {args.resolutions}   (tables at {args.ref_resolution})")
    print(f"list depth : top {args.n_top_genes} genes, both arms", flush=True)

    if args.ref_resolution not in args.resolutions:
        sys.exit(f"--ref-resolution {args.ref_resolution} is not in --resolutions "
                 f"{args.resolutions}: the per-programme tables would have no partition")

    # ---- the gene sets, and the denominator neither arm chose ---------------
    gmt = C.gmt_path(coll)
    if not gmt.exists():
        sys.exit(f"missing {gmt}: run 05_4_signatures/build_signatures_tum.py "
                 f"--collection {args.collection} first")
    sets_all = C.read_gmt(gmt)

    hvg = ad.read_h5ad(C.HVG_H5AD)
    hvg.obs_names = hvg.obs_names.astype(str)
    if args.states_argmax:
        drvi_states_argmax(args, drvi_run)
        drvi_states_de(hvg, args, drvi_run, argmax=True)
        drvi_states_dr(args, drvi_run, argmax=True)
        print("\ndone (--states-argmax).")
        return
    if args.states_dr:
        drvi_states_dr(args, drvi_run)
        print("\ndone (--states-dr).")
        return
    if args.states_de:
        drvi_states_de(hvg, args, drvi_run)
        print("\ndone (--states-de).")
        return
    if args.states_only:
        drvi_states(hvg, args, harmony_run, drvi_run)
        print("\ndone (--states-only).")
        return
    background = list(hvg.var_names)
    bg = set(background)
    sets_bg = {k: [g for g in v if g in bg] for k, v in sets_all.items()}
    testable = {k: v for k, v in sets_bg.items() if len(v) >= C.MIN_SIGNATURE_GENES}
    dropped = {k: len(v) for k, v in sets_bg.items() if k not in testable}

    C.banner("the denominator: every signature either arm could in principle enrich for")
    print(f"{hvg.n_obs:,} cells x {len(background):,} HVGs - the DRVI training feature set, "
          f"the DE feature set and the ORA background, all three")
    print(f"{len(sets_all)} signatures in the collection -> {len(testable)} with at least "
          f"{C.MIN_SIGNATURE_GENES} genes inside it")
    for k, n in dropped.items():
        print(f"  dropped: {k:32s} {n:3d} genes in background")

    if args.presence_filter:
        pres = C.table_path(f"signature_presence", coll, run_id)
        if not pres.exists():
            print(f"[skip] --presence-filter: {pres.name} is not on disk, the denominator "
                  f"stays the testable set", flush=True)
        else:
            t = pd.read_csv(pres, comment="#")
            keep = set(t.loc[(t["space"] == "harmony") & (t["morans_i_vs_random_z"] >= 3),
                             "readout"])
            before = len(testable)
            testable = {k: v for k, v in testable.items() if k in keep}
            print(f"--presence-filter: {before} -> {len(testable)} signatures 05_9 found "
                  f"present in the Harmony space at >= 3 sd over the matched random level")

    if not testable:
        sys.exit("no signature is testable on this background; nothing to measure")
    n_testable = len(testable)

    # ---- the cells, the cohorts, the cached per-cell scores -----------------
    cohorts = hvg.obs[C.BATCH_KEY].astype(str)
    scores_csv = C.scores_csv(coll)
    z = None
    if scores_csv.exists():
        s = pd.read_csv(scores_csv, index_col=0)
        s.index = s.index.astype(str)
        z = s[[c for c in s.columns if c.startswith("z_")]].rename(columns=lambda c: c[2:])
        z = z.reindex(hvg.obs_names)
        print(f"\n[read] {scores_csv.name}: {z.shape[1]} per-cell readouts, "
              f"05_6's within-cohort z - re-scored nowhere in this step", flush=True)
    else:
        print(f"\n[skip] {scores_csv.name} is not on disk, so eta^2(cluster) cannot be "
              f"computed; run 05_6 for the diagnostic half of this step", flush=True)

    # ---- the gene lists, arm by arm -----------------------------------------
    # {(arm, resolution): {list name: [genes]}}, plus the partitions the cluster arms used.
    gene_lists: dict[tuple, dict] = {}
    de_stats: dict[tuple, pd.DataFrame] = {}
    partitions: dict[tuple, pd.Series] = {}
    descriptions: dict[str, str] = {}

    for arm in args.arms:
        if arm not in LEIDEN_ARMS:
            continue
        C.banner(f"arm {arm}")
        cached, cparts, desc = load_cache(arm, drvi_run, args, hvg.obs_names)
        todo = [r for r in args.resolutions if r not in cached]
        if todo:
            coords, space_desc = load_coords(ARM_SPACE[arm], harmony_run, drvi_run)
            if coords is None:
                print(f"[skip] arm {arm}: its space could not be read", flush=True)
                continue
            absent = pd.Index(hvg.obs_names).difference(coords.index)
            if len(absent):
                sys.exit(f"{arm}: {len(absent)} cells of the object are not in its embedding "
                         f"(e.g. {list(absent[:3])}). Two pipelines compared on two cell "
                         f"sets is not a comparison")
            desc = space_desc
            print(f"    {desc}", flush=True)
            labs = leiden_labels(hvg, coords, todo, args.k, C.SEED, arm)
            for res, lab in labs.items():
                cparts[res] = lab.astype(str)
                cached[res] = de_top_genes(hvg, lab, args.n_top_genes, args.de_fdr,
                                           args.de_min_lfc)
            save_cache(arm, drvi_run, args, desc, cached, cparts)
        else:
            print(f"    {desc or ARM_LABEL.get(arm, arm)}", flush=True)
        descriptions[arm] = desc
        for res in args.resolutions:
            gene_lists[(arm, res)], de_stats[(arm, res)] = cached[res]
            partitions[(arm, res)] = cparts[res]

    if "drvi_leiden" in args.arms and args.state_resolutions:
        drvi_states(hvg, args, harmony_run, drvi_run)

    if "drvi_axes" in args.arms:
        embed = ad.read_h5ad(CS.tum_dir() / f"embed_{drvi_run}.h5ad")
        n_van = C.n_vanished(embed)
        desc = (f"05_3: DRVI decoder, {embed.n_vars} dimensions x 2 directions "
                + (f"({n_van} vanished PRUNED)" if C.PRUNE_VANISHED
                   else f"({n_van} flagged vanished and NOT pruned)"))
        descriptions["drvi_axes"] = desc
        C.banner(f"arm drvi_axes: {desc}")
        assert embed.varm[f"{C.SCORE_KEY}_positive"].shape[1] == len(background), \
            "the embedding and the DRVI input disagree on the gene axis"
        gene_lists[("drvi_axes", np.nan)] = axis_top_genes(embed, hvg.var_names,
                                                           args.n_top_genes)

    # THE BUDGET OBJECTION, CLOSED RATHER THAN ARGUED. The decoder brings 2 x n_latent lists
    # and the widest Leiden partition brings fewer. This arm is the decoder cut down to the
    # widest cluster count in the run, taking DRVI's OWN top dimensions - `var['order']`, the
    # reconstruction-effect ranking, both directions of each - so it is the subset the model
    # itself would have offered first, not one chosen to win. If it still names what the
    # cluster arms cannot at the same number of tests, the gap is not a test budget.
    if ("drvi_axes", np.nan) in gene_lists:
        widest = max((len(v) for (a, _), v in gene_lists.items() if a in LEIDEN_ARMS),
                     default=0)
        full = gene_lists[("drvi_axes", np.nan)]
        if 0 < widest < len(full):
            keep = list(full)[:widest]                 # already in `order`, +/- interleaved
            gene_lists[("drvi_axes_matched", np.nan)] = {k: full[k] for k in keep}
            descriptions["drvi_axes_matched"] = (
                f"the first {widest} of the {len(full)} decoder directions, in DRVI's own "
                f"reconstruction-effect order - the widest Leiden partition of this run")
            print(f"\n[budget control] drvi_axes_matched: {widest} directions, "
                  f"matching the widest cluster count in the run", flush=True)

    if not args.no_null and partitions:
        ref_arm = next((a for a in args.arms if (a, args.ref_resolution) in partitions), None)
        if ref_arm is not None:
            C.banner("the floor: the reference partition shuffled among the cells")
            lab = partitions[(ref_arm, args.ref_resolution)]
            rng = np.random.default_rng(C.SEED)
            shuffled = pd.Series(rng.permutation(lab.values), index=lab.index)
            print(f"    {shuffled.nunique()} groups, every group size preserved, "
                  f"membership destroyed", flush=True)
            descriptions["null_partition"] = (
                f"the {ref_arm} partition at resolution {args.ref_resolution}, shuffled")
            gene_lists[("null_partition", args.ref_resolution)], \
                de_stats[("null_partition", args.ref_resolution)] = de_top_genes(
                    hvg, shuffled, args.n_top_genes, args.de_fdr, args.de_min_lfc)
            partitions[("null_partition", args.ref_resolution)] = shuffled

    if not gene_lists:
        sys.exit("no arm produced a gene list; nothing to test")

    # ---- one ORA, every arm -------------------------------------------------
    C.banner(f"ORA: offline hypergeometric against the {len(background):,}-HVG background, "
             f"BH within arm")
    summary_rows, per_list_rows, long_all = [], [], []

    for (arm, res), lists in gene_lists.items():
        long = run_ora(lists, testable, background, arm, res)
        long, n_pairs = apply_bh(long, len(lists), n_testable)
        long_all.append(long)

        hits = long[long["significant"]] if len(long) else long
        recovered = sorted(set(hits["Term"])) if len(hits) else []
        # A unit is IDENTIFIED BY its strongest enrichment; the distinct identities are the
        # arm's vocabulary. This is 05_8's own definition of naming a dimension, applied
        # unchanged to a cluster.
        best_per_list = (hits.loc[hits.groupby("gene_list", observed=True)["fdr_bh"].idxmin()]
                         if len(hits) else pd.DataFrame(columns=["gene_list", "Term"]))
        named = sorted(set(best_per_list["Term"])) if len(best_per_list) else []
        per_list_hits = (hits.groupby("gene_list", observed=True)["Term"].nunique()
                         if len(hits) else pd.Series(dtype=int))
        per_list_hits = per_list_hits.reindex(list(lists), fill_value=0)

        prof = None
        if (arm, res) in partitions:
            prof = cluster_profile(partitions[(arm, res)], cohorts)

        list_sizes = np.array([len(v) for v in lists.values()])
        row = {
            "arm": arm,
            "description": descriptions.get(arm, ""),
            "unit": "axis" if arm.startswith("drvi_axes") else "cluster",
            "resolution": res,
            "n_gene_lists": len(lists),
            "n_gene_lists_empty": int((list_sizes == 0).sum()),
            "median_genes_per_list": float(np.median(list_sizes)),
            "n_tests": n_pairs,
            "n_signatures_testable": n_testable,
            "n_signatures_recovered": len(recovered),
            "recall_any": len(recovered) / n_testable,
            "n_signatures_named": len(named),
            "recall_named": len(named) / n_testable,
            "n_units_with_a_name": int(len(best_per_list)),
            "n_significant_pairs": int(len(hits)),
            "mean_hits_per_list": float(per_list_hits.mean()),
            "max_hits_per_list": int(per_list_hits.max()) if len(per_list_hits) else 0,
            "frac_lists_with_no_hit": float((per_list_hits == 0).mean()),
            "median_best_fdr": (float(hits.groupby("Term")["fdr_bh"].min().median())
                                if len(hits) else float("nan")),
            "recovered": "; ".join(recovered),
            "named": "; ".join(named),
        }
        if prof is not None:
            row.update({
                "min_cluster_cells": int(prof["n_cells"].min()),
                "median_cluster_cells": float(prof["n_cells"].median()),
                "max_cohort_dominance": float(prof["cohort_dominance"].max()),
                "n_clusters_one_cohort_over_80pct":
                    int((prof["cohort_dominance"] > 0.8).sum()),
            })
        summary_rows.append(row)

        n_de = ({r["cluster"]: r["n_de_genes"]
                 for _, r in de_stats[(arm, res)].iterrows()}
                if (arm, res) in de_stats else {})
        for name, n in per_list_hits.items():
            per_list_rows.append({"arm": arm, "resolution": res, "gene_list": str(name),
                                  "n_genes_in_list": len(lists[name]),
                                  "n_de_genes": n_de.get(str(name), np.nan),
                                  "n_signatures_hit": int(n),
                                  "n_cells": (int((partitions[(arm, res)] == str(name)).sum())
                                              if (arm, res) in partitions else np.nan)})

        print(f"  {arm:16s} res {str(res):<5} {len(lists):>4} lists  "
              f"named {len(named):>3}/{n_testable} (recall_named {row['recall_named']:.3f})  "
              f"any {len(recovered):>3}/{n_testable}  "
              f"mean hits/list {row['mean_hits_per_list']:.2f}", flush=True)

    summary = pd.DataFrame(summary_rows).sort_values(["unit", "arm", "resolution"])
    # `recall` is kept as an alias of the headline so nothing downstream has to know which
    # of the two it wanted; both columns are in the table and the figures name them.
    summary["recall"] = summary["recall_named"]
    per_list = pd.DataFrame(per_list_rows)
    long_all = pd.concat(long_all, ignore_index=True) if long_all else pd.DataFrame()

    # ---- the union over the whole scan, the most generous reading -----------
    #
    # A practitioner picks ONE resolution. This row gives the cluster arms every resolution
    # at once - the union of everything they name anywhere in the scan - which is an upper
    # bound no single run of that pipeline could reach, and is labelled as one. A programme
    # absent from this union is a programme the standard pipeline does not reach on these
    # cells at any granularity it was given.
    union_rows = []
    for arm, grp in summary.groupby("arm", observed=True):
        named_union = set()
        for v in grp["named"].dropna():
            named_union |= {t for t in str(v).split("; ") if t}
        union_rows.append({
            "arm": arm,
            "n_settings": int(len(grp)),
            "n_signatures_named_union": len(named_union),
            "recall_named_union": len(named_union) / n_testable,
            "best_single_setting": float(grp["recall_named"].max()),
            "named_union": "; ".join(sorted(named_union)),
            "never_named": "; ".join(sorted(set(testable) - named_union)),
        })
    union = pd.DataFrame(union_rows).sort_values("recall_named_union", ascending=False)

    # ---- the per-programme table, at the reference resolution ---------------
    ref = {}
    extra = (["drvi_axes_matched"] if ("drvi_axes_matched", np.nan) in gene_lists else []) \
        + (["null_partition"] if not args.no_null else [])
    for arm in args.arms + extra:
        key = (arm, np.nan) if arm.startswith("drvi_axes") else (arm, args.ref_resolution)
        if key in gene_lists:
            ref[arm] = key

    rec = pd.DataFrame(index=sorted(testable))
    rec.index.name = "signature"
    rec["n_genes_in_collection"] = [len(sets_all[s]) for s in rec.index]
    rec["n_genes_in_background"] = [len(testable[s]) for s in rec.index]

    for arm, key in ref.items():
        sub = long_all[(long_all["arm"] == key[0]) &
                       ((long_all["resolution"] == key[1])
                        if not pd.isna(key[1]) else long_all["resolution"].isna())]
        best = (sub.loc[sub.groupby("Term")["fdr_bh"].idxmin()].set_index("Term")
                if len(sub) else pd.DataFrame())
        fdr = best["fdr_bh"].reindex(rec.index) if len(best) else pd.Series(index=rec.index,
                                                                           dtype=float)
        rec[f"recovered__{arm}"] = (fdr < FDR).fillna(False)
        rec[f"best_fdr__{arm}"] = fdr
        rec[f"neglog_fdr__{arm}"] = -np.log10(fdr.clip(lower=1e-300))
        rec[f"best_list__{arm}"] = (best["gene_list"].reindex(rec.index)
                                    if len(best) else pd.NA)
        sig = sub[sub["significant"]] if len(sub) else sub
        rec[f"n_lists_hitting__{arm}"] = (
            sig.groupby("Term")["gene_list"].nunique().reindex(rec.index).fillna(0)
            .astype(int) if len(sig) else 0)
        bpl = (sig.loc[sig.groupby("gene_list", observed=True)["fdr_bh"].idxmin()]
               if len(sig) else pd.DataFrame(columns=["gene_list", "Term"]))
        n_naming = (bpl.groupby("Term")["gene_list"].nunique().reindex(rec.index).fillna(0)
                    .astype(int) if len(bpl) else 0)
        rec[f"n_lists_naming__{arm}"] = n_naming
        rec[f"named__{arm}"] = (n_naming > 0) if len(bpl) else False

    # eta^2(cluster), the diagnostic
    eta_rows = []
    if z is not None:
        for (arm, res), lab in partitions.items():
            codes, uniq = pd.factorize(lab.reindex(hvg.obs_names).values)
            for sig in rec.index:
                if sig not in z.columns:
                    continue
                y = z[sig].values
                ok = np.isfinite(y)
                eta_rows.append({"arm": arm, "resolution": res, "signature": sig,
                                 "n_clusters": len(uniq),
                                 "eta2_cluster": eta_squared(y[ok], codes[ok], len(uniq))})
        eta = pd.DataFrame(eta_rows)
        for arm, key in ref.items():
            if arm.startswith("drvi_axes"):
                continue
            sub = eta[(eta["arm"] == key[0]) & (eta["resolution"] == key[1])]
            rec[f"eta2_cluster__{arm}"] = sub.set_index("signature")["eta2_cluster"].reindex(
                rec.index)
    else:
        eta = pd.DataFrame(columns=["arm", "resolution", "signature", "n_clusters",
                                    "eta2_cluster"])

    # ---- tables --------------------------------------------------------------
    C.banner("tables")
    write_table(summary, "pipeline_recall", coll, run_id, index=False)
    write_table(union, "pipeline_recall_union", coll, run_id, index=False)
    write_table(rec.reset_index(), "programme_recovery", coll, run_id, index=False)
    write_table(per_list, "gene_list_hits", coll, run_id, index=False)
    if len(eta):
        write_table(eta, "programme_cluster_eta2", coll, run_id, index=False)

    comp_rows = []
    for arm, key in ref.items():
        if key not in partitions:
            continue
        prof = cluster_profile(partitions[key], cohorts).assign(arm=arm, resolution=key[1])
        hits = long_all[(long_all["arm"] == key[0]) & (long_all["resolution"] == key[1])
                        & long_all["significant"]]
        by_cl = hits.groupby("gene_list", observed=True)["Term"].apply(
            lambda s: "; ".join(sorted(set(s)))) if len(hits) else pd.Series(dtype=str)
        if key in de_stats:
            prof = prof.merge(de_stats[key], on="cluster", how="left")
        prof["signatures_hit"] = prof["cluster"].map(by_cl).fillna("")
        prof["n_signatures_hit"] = prof["cluster"].map(
            by_cl.str.count("; ").add(1) if len(by_cl) else pd.Series(dtype=float)).fillna(0)
        comp_rows.append(prof)
    if comp_rows:
        write_table(pd.concat(comp_rows, ignore_index=True), "cluster_composition", coll,
                    run_id, index=False)

    if len(long_all):
        write_table(long_all.drop(columns=["Genes"], errors="ignore"), "ora_all_pairs", coll,
                    run_id, index=False)

    # ---- figures -------------------------------------------------------------
    C.banner("figures")
    fig_recall_curve(summary, coll, run_id, n_testable, args)
    fig_recovery_matrix(rec, coll, run_id, list(ref), args.ref_resolution)
    for arm in ref:
        if arm in LEIDEN_ARMS:
            fig_eta2(rec, coll, run_id, arm, args.ref_resolution)
    fig_hits_per_list(per_list[(per_list["resolution"] == args.ref_resolution)
                               | per_list["resolution"].isna()],
                      coll, run_id, list(ref), args.ref_resolution)

    # ---- the read ------------------------------------------------------------
    C.banner("what came out")
    print(f"  denominator: {n_testable} testable signatures of {len(sets_all)}\n")
    best = (summary[summary["unit"] == "cluster"].groupby("arm")["recall_named"].max()
            if (summary["unit"] == "cluster").any() else pd.Series(dtype=float))
    for arm, r in best.items():
        row = summary[(summary["arm"] == arm) & (summary["recall_named"] == r)].iloc[0]
        print(f"  {ARM_LABEL.get(arm, arm):34s} named {int(row['n_signatures_named']):>3}"
              f"/{n_testable} (best, at resolution {row['resolution']} with "
              f"{int(row['n_gene_lists'])} clusters)   any {int(row['n_signatures_recovered'])}"
              f"/{n_testable}   {row['mean_hits_per_list']:.2f} signatures per cluster")
    for _, row in summary[summary["unit"] == "axis"].iterrows():
        print(f"  {ARM_LABEL.get(row['arm'], row['arm']):34s} named "
              f"{int(row['n_signatures_named']):>3}/{n_testable} (with "
              f"{int(row['n_gene_lists'])} directions)   any "
              f"{int(row['n_signatures_recovered'])}/{n_testable}   "
              f"{row['mean_hits_per_list']:.2f} signatures per direction")

    print("\n  union over the whole scan - every resolution at once, an upper bound no "
          "single run\n  of the clustering pipeline could reach:")
    for _, r in union.iterrows():
        never = r["never_named"]
        print(f"    {ARM_LABEL.get(r['arm'], r['arm']):34s} "
              f"{int(r['n_signatures_named_union']):>3}/{n_testable} "
              f"(best single setting {r['best_single_setting']:.3f})")
        if never and r["arm"] in LEIDEN_ARMS:
            print(f"        never named at any resolution: {never}")

    if "drvi_axes" in ref and any(a in ref for a in LEIDEN_ARMS):
        axis_hit = set(rec.index[rec["named__drvi_axes"].astype(bool)])
        for arm in [a for a in LEIDEN_ARMS if a in ref]:
            cl_hit = set(rec.index[rec[f"named__{arm}"].astype(bool)])
            only_axis = sorted(axis_hit - cl_hit)
            only_cl = sorted(cl_hit - axis_hit)
            print(f"\n  {ARM_LABEL.get(arm, arm)} at resolution {args.ref_resolution}:")
            print(f"    named by the decoder and NOT by the clusters ({len(only_axis)}): "
                  f"{', '.join(only_axis) or '-'}")
            print(f"    named by the clusters and NOT by the decoder ({len(only_cl)}): "
                  f"{', '.join(only_cl) or '-'}")
            col = f"eta2_cluster__{arm}"
            if col in rec.columns and rec[col].notna().any():
                miss = rec.loc[[s for s in only_axis if s in rec.index], col].dropna()
                hit = rec.loc[[s for s in sorted(cl_hit) if s in rec.index], col].dropna()
                if len(miss) and len(hit):
                    print(f"    median eta^2(cluster): {miss.median():.4f} for the ones it "
                          f"could not name, {hit.median():.4f} for the ones it named")

    print("\ndone.")


if __name__ == "__main__":
    main()
