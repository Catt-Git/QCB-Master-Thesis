#!/usr/bin/env python3
"""05_1 step 1: build the per-patient inferCNV inputs from shiao.h5ad.

inferCNV calls copy-number variation by comparing the smoothed expression of the cells
under test against a REFERENCE set of cells assumed to be karyotypically normal. This
dataset has no such set given to us: `treatment` is BASE/PD1/RTPD1, i.e. three timepoints
of the same tumour, and there is no normal or adjacent tissue anywhere in the 34 cohorts.
The reference therefore has to come from inside each tumour, and the standard choice
applies - immune cells, which are diploid and are here in abundance (T + myeloid never
drop below 1,272 cells in any cohort).

What each per-patient run contains, and why:

  reference  ref_tcell    up to N_REF_PER_GROUP T/NK cells      diploid baseline
             ref_myeloid  up to N_REF_PER_GROUP myeloid cells   diploid baseline
  observed   epi          every epithelial cell of the cohort   what we want to call
             stromal      up to N_STROMAL fibro/vascular cells  internal negative control

Two reference groups rather than one pooled group is inferCNV's own recommendation: with
`ref_group_names` of length > 1 the residual of a gene is taken against the *bounds* of
the per-group means, so a gene that is simply higher in myeloid cells than in T cells
cannot masquerade as a gain. The stromal block is not a reference: it is passed as an
OBSERVATION so that it goes through exactly the same smoothing and denoising as the
epithelium, and its CNV score distribution becomes a free specificity check - fibroblasts
and endothelium are not the malignant compartment in a carcinoma, so if they score like
the epithelium the call is measuring something other than aneuploidy.

Per patient, not pooled. Every published application of inferCNV to a multi-patient
cohort runs one patient at a time, for the same reason 04_4 scores CytoTRACE2 per
patient: the residual is defined against the reference cells present in the run, so a
pooled run would compare patient A's epithelium against patient B's immune cells and
read the batch difference as copy number. It also keeps every run small enough to hold
in memory.

Gene ordering file: hg38 / GENCODE v27, the file the inferCNV authors distribute
(https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt), downloaded once
into $DATA_DIR/05_tum/ and reused. It is keyed by gene symbol, matches the GRCh38
symbols in `var_names` on 30,226 of the 30,869 genes (97.9%), and has no duplicated
symbol on the main chromosomes. Some of its lines carry a trailing tab, hence the
explicit `usecols`.

Input : $DATA_DIR/shiao.h5ad  (read backed; `layers['counts']` = raw integer counts)
Output: $DATA_DIR/05_tum/gene_order_hg38_gencode_v27.txt
        $DATA_DIR/05_tum/input/<cohort>/counts.mtx      genes x cells, integer
        $DATA_DIR/05_tum/input/<cohort>/genes.tsv
        $DATA_DIR/05_tum/input/<cohort>/barcodes.tsv
        $DATA_DIR/05_tum/input/<cohort>/annotations.tsv cell id -> group
        $DATA_DIR/05_tum/cohort_census.csv              what went into each run

Nothing under $DATA_DIR/04_epi is read or written. Phase 04 is a finished, frozen branch off
`shiao.h5ad` that predates any CNV call; phase 05 is a second branch off the same object, and
its outputs live in their own directory.

Local usage (benchmark-py-r):
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    python3 prepare_infercnv_input.py
    python3 prepare_infercnv_input.py --cohorts Patient52 Patient16   # a subset
    python3 prepare_infercnv_input.py --force                        # rewrite existing
"""

from __future__ import annotations

import argparse
import os
import sys
import urllib.request
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

# --------------------------------------------------------------------------------------
# Lineage sets. Verbatim the sets of 01_4/fraction_reassignment.py, split further into the
# four groups the inferCNV run needs. Their union must be the whole label vocabulary, and
# that is asserted below rather than assumed: a label missing from all four would silently
# vanish from the run.
# --------------------------------------------------------------------------------------
EPITHELIAL = {
    "LummHR-SCGB", "LummHR-active", "LummHR-major",
    "Lumsec-HLA", "Lumsec-KIT", "Lumsec-basal", "Lumsec-lac", "Lumsec-major",
    "Lumsec-myo", "Lumsec-prol",
    "basal",
}
STROMAL = {
    "Fibro-SFRP4", "Fibro-major", "Fibro-matrix", "Fibro-prematrix",
    "Lymph-immune", "Lymph-major", "Lymph-valve1", "Lymph-valve2",
    "Vas-arterial", "Vas-capillary", "Vas-venous", "pericytes", "vsmc",
}
T_NK = {
    "CD4-Tem", "CD4-Th", "CD4-Th-like", "CD4-Treg", "CD4-activated", "CD4-naive",
    "CD8-Tem", "CD8-Trm", "CD8-activated", "GD", "NK", "NK-ILCs", "NKT", "T_prol",
}
MYELOID = {
    "Macro-IFN", "Macro-lipo", "Macro-m1", "Macro-m1-CCL", "Macro-m2", "Macro-m2-CXCL",
    "Mast", "Mono-classical", "Mono-non-classical", "Neutrophil",
    "cDC1", "cDC2", "mDC", "pDC", "mye-prol",
}
# B / plasma cells are deliberately in NEITHER reference group. Plasma cells carry a huge,
# clonally skewed immunoglobulin transcriptome (IGH/IGK/IGL sit on chr14/chr2/chr22) and
# there are 62,074 of them here; using them as a diploid baseline would put a spurious
# structure on three chromosomes. They are simply left out of every run.
B_PLASMA = {"b_naive", "bmem_switched", "bmem_unswitched", "plasma_IgA", "plasma_IgG"}

GROUP_EPI = "epi"
GROUP_STROMAL = "stromal"
GROUP_REF_T = "ref_tcell"
GROUP_REF_MYE = "ref_myeloid"
REFERENCE_GROUPS = (GROUP_REF_T, GROUP_REF_MYE)

BATCH_KEY = "cohort"
LABEL_KEY = "cell_type"

N_REF_PER_GROUP = 1000   # T/NK and myeloid cells drawn per cohort, each
N_STROMAL = 1000         # stromal cells carried as the internal negative control
MIN_EPI_CELLS = 50       # a cohort with fewer epithelial cells is not run at all
MIN_CELLS_PER_GENE = 3   # drop genes detected in < this many cells OF THIS COHORT
SEED = 0

GENE_ORDER_URL = "https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt"
GENE_ORDER_NAME = "gene_order_hg38_gencode_v27.txt"
MAIN_CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--cohorts", nargs="+", default=None,
                   help="only these cohorts (default: every cohort with enough epithelium)")
    p.add_argument("--force", action="store_true",
                   help="rewrite inputs for cohorts that already have them")
    return p.parse_args()


def fetch_gene_order(cnv_dir: Path) -> pd.DataFrame:
    """Download the inferCNV hg38 gene ordering file once, return it parsed."""
    path = cnv_dir / GENE_ORDER_NAME
    if not path.exists():
        print(f"Downloading gene order file -> {path}", flush=True)
        tmp = path.with_suffix(".tmp")
        urllib.request.urlretrieve(GENE_ORDER_URL, tmp)
        os.replace(tmp, path)
    go = pd.read_csv(path, sep="\t", header=None,
                     names=["gene", "chr", "start", "stop", "_trailing"],
                     usecols=[0, 1, 2, 3])
    go = go[go["chr"].isin(MAIN_CHROMS)]
    assert not go["gene"].duplicated().any(), "duplicated gene symbol in the ordering file"
    # inferCNV wants the file sorted by position; sort explicitly rather than trusting it.
    go["chr"] = pd.Categorical(go["chr"], categories=MAIN_CHROMS, ordered=True)
    go = go.sort_values(["chr", "start", "stop"]).set_index("gene")
    print(f"Gene order: {len(go):,} genes on {go['chr'].nunique()} chromosomes", flush=True)
    return go


def assign_group(labels: pd.Series) -> pd.Series:
    """Map `cell_type` onto the four inferCNV groups; B/plasma and anything else -> NaN."""
    out = pd.Series(np.nan, index=labels.index, dtype=object)
    out[labels.isin(EPITHELIAL)] = GROUP_EPI
    out[labels.isin(STROMAL)] = GROUP_STROMAL
    out[labels.isin(T_NK)] = GROUP_REF_T
    out[labels.isin(MYELOID)] = GROUP_REF_MYE
    return out


def main() -> int:
    args = parse_args()
    rng = np.random.default_rng(SEED)

    data_dir = Path(os.environ["DATA_DIR"]).expanduser().resolve()
    in_path = data_dir / "shiao.h5ad"
    cnv_dir = data_dir / "05_tum"
    input_dir = cnv_dir / "input"
    cnv_dir.mkdir(parents=True, exist_ok=True)
    input_dir.mkdir(exist_ok=True)

    print(f"DATA_DIR: {data_dir}")
    print(f"Input   : {in_path} ({in_path.stat().st_size / 1024**3:.2f} GB)")
    print(f"Output  : {input_dir}")
    print()

    gene_order = fetch_gene_order(cnv_dir)

    # Backed, because only .obs and one cohort's slice of the counts are ever needed at once.
    adata = ad.read_h5ad(in_path, backed="r")
    assert "counts" in adata.layers, "expected raw counts in layers['counts']"

    labels = adata.obs[LABEL_KEY].astype(str)
    unclassified = set(labels.unique()) - EPITHELIAL - STROMAL - T_NK - MYELOID - B_PLASMA
    assert not unclassified, (
        f"cell_type value(s) in no lineage set: {sorted(unclassified)}; "
        "update the sets at the top of this script before running"
    )
    groups = assign_group(labels)

    # Genes are intersected with the ordering file ONCE, on the full object, so every
    # cohort is run on the same feature space and the CNV scores stay comparable.
    keep_genes = adata.var_names.intersection(gene_order.index)
    keep_genes = pd.Index(gene_order.index[gene_order.index.isin(keep_genes)])  # positional order
    gene_pos = pd.Series(np.arange(adata.n_vars), index=adata.var_names)[keep_genes].to_numpy()
    print(f"Genes: {adata.n_vars:,} in the object, {len(keep_genes):,} also in the gene order file")
    print()

    cohorts = list(adata.obs[BATCH_KEY].cat.categories) if args.cohorts is None else args.cohorts
    census_rows = []

    for cohort in cohorts:
        in_cohort = (adata.obs[BATCH_KEY].astype(str) == cohort).to_numpy()
        grp = groups[in_cohort]
        n_epi = int((grp == GROUP_EPI).sum())

        if n_epi < MIN_EPI_CELLS:
            print(f"[skip] {cohort}: {n_epi} epithelial cells (< {MIN_EPI_CELLS})")
            census_rows.append(dict(cohort=cohort, n_epi=n_epi, status="skipped_too_few_epi"))
            continue

        out_dir = input_dir / cohort
        annot_path = out_dir / "annotations.tsv"
        if annot_path.exists() and not args.force:
            # Still census it, reading the counts back off the annotations that are already
            # there: a resumed run that only appended the cohorts it rewrote would drop the
            # ones it skipped out of the census, and stage 2 of infercnv_all.sh reads the
            # census to decide what to run.
            have = pd.read_csv(annot_path, sep="\t", header=None, names=["cell", "group"])
            have_n = have["group"].value_counts()
            print(f"[have] {cohort}: inputs already written, skipping")
            census_rows.append(dict(
                cohort=cohort, status="prepared",
                n_cells=len(have),
                n_genes=sum(1 for _ in open(out_dir / "genes.tsv")),
                **{f"n_{g}": int(have_n.get(g, 0)) for g in
                   (GROUP_REF_T, GROUP_REF_MYE, GROUP_STROMAL, GROUP_EPI)},
            ))
            continue

        # Pick the cells: all epithelium, a capped stromal control, a capped reference.
        picked = {}
        for group, cap in ((GROUP_EPI, None), (GROUP_STROMAL, N_STROMAL),
                           (GROUP_REF_T, N_REF_PER_GROUP), (GROUP_REF_MYE, N_REF_PER_GROUP)):
            idx = np.flatnonzero((grp == group).to_numpy())
            if cap is not None and len(idx) > cap:
                idx = rng.choice(idx, size=cap, replace=False)
            picked[group] = np.sort(idx)

        # Positions are relative to the cohort slice; lift them back to the full object.
        cohort_pos = np.flatnonzero(in_cohort)
        order = np.concatenate([picked[g] for g in
                                (GROUP_REF_T, GROUP_REF_MYE, GROUP_STROMAL, GROUP_EPI)])
        rows = cohort_pos[order]
        cell_groups = np.concatenate([[g] * len(picked[g]) for g in
                                      (GROUP_REF_T, GROUP_REF_MYE, GROUP_STROMAL, GROUP_EPI)])
        cell_ids = adata.obs_names[rows]

        # Slice the counts. `.layers['counts']` on a backed object needs an increasing
        # index, which `rows` is by construction (cohort_pos sorted, each block sorted,
        # blocks not interleaved only if the picks happen to be ordered) - so sort here and
        # reorder the labels to match, rather than relying on that.
        srt = np.argsort(rows, kind="stable")
        counts = adata.layers["counts"][rows[srt], :]
        counts = sp.csr_matrix(counts)[:, gene_pos]
        cell_ids = cell_ids[srt]
        cell_groups = cell_groups[srt]

        # Per-cohort gene filter: a gene seen in fewer than MIN_CELLS_PER_GENE cells of THIS
        # patient carries no CNV signal and only lengthens the matrix inferCNV has to smooth.
        detected = np.asarray((counts > 0).sum(axis=0)).ravel()
        gene_mask = detected >= MIN_CELLS_PER_GENE
        counts = counts[:, gene_mask]
        genes_here = keep_genes[gene_mask]

        assert counts.shape == (len(cell_ids), len(genes_here))
        n_by_group = pd.Series(cell_groups).value_counts()

        out_dir.mkdir(parents=True, exist_ok=True)
        # inferCNV wants genes as ROWS and cells as COLUMNS.
        sio.mmwrite(str(out_dir / "counts.mtx"), counts.T.astype(np.int32), field="integer")
        pd.Series(genes_here).to_csv(out_dir / "genes.tsv", index=False, header=False)
        pd.Series(cell_ids).to_csv(out_dir / "barcodes.tsv", index=False, header=False)
        # annotations.tsv is inferCNV's own format: cell id, tab, group. No header.
        pd.DataFrame({"cell": cell_ids, "group": cell_groups}).to_csv(
            out_dir / "annotations.tsv", sep="\t", index=False, header=False)

        mtx_mb = (out_dir / "counts.mtx").stat().st_size / 1024**2
        print(f"[ok]   {cohort}: {counts.shape[0]:,} cells x {counts.shape[1]:,} genes "
              f"({mtx_mb:.0f} MB) | " +
              " ".join(f"{g}={n_by_group.get(g, 0)}" for g in
                       (GROUP_REF_T, GROUP_REF_MYE, GROUP_STROMAL, GROUP_EPI)), flush=True)

        census_rows.append(dict(
            cohort=cohort, status="prepared",
            n_cells=counts.shape[0], n_genes=counts.shape[1],
            # n_epi is the GROUP_EPI count: the epithelium is never capped, so the two
            # are the same number and the skipped rows below use the same column name.
            **{f"n_{g}": int(n_by_group.get(g, 0)) for g in
               (GROUP_REF_T, GROUP_REF_MYE, GROUP_STROMAL, GROUP_EPI)},
        ))

    census = pd.DataFrame(census_rows)
    census_path = cnv_dir / "cohort_census.csv"
    if census_path.exists() and args.cohorts is not None:
        # A partial re-run must not truncate the census of the cohorts it did not touch.
        old = pd.read_csv(census_path)
        census = pd.concat([old[~old["cohort"].isin(census["cohort"])], census])
    census.sort_values("cohort").to_csv(census_path, index=False)
    print()
    print(f"Wrote {census_path}")
    prepared = census[census["status"] == "prepared"] if len(census) else census
    print(f"{len(prepared)} cohort(s) ready for run_infercnv.R")
    return 0


if __name__ == "__main__":
    sys.exit(main())
