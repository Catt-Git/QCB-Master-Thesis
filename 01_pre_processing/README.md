# 01_pre_processing

Second phase of the thesis: from the concatenated raw-count object produced in
`00_raw_data_processing` to a single **quality-controlled, normalized, annotated and
dimensionally-reduced** unintegrated AnnData (`shiao.h5ad`), ready for the integration
benchmark. Each step reads the previous step's `.h5ad` and writes the next one.

Dataset: Shiao et al., *Cancer Cell* 2024 (DOI 10.1016/j.ccell.2023.12.012),
BioProject **PRJNA1032700** / GEO **GSE246613**.

## Execution order

| # | Folder | Script | What it does | Where |
|---|--------|--------|--------------|-------|
| 1 | `01_1_load_and_check` | `load_and_check.ipynb` | Sanity-check the concatenated object: raw-count `.X`/`counts`, unique gene symbols, metadata fields; `cohort × response/treatment/fraction` crosstabs (the `batch_key` decision), NaN checks. No output file. | local (notebook) |
| 2 | `01_2_qc_scrublet_filtering` | `qc_scrublet_filter.ipynb` | QC metrics (mt/ribo), per-biopsy Scrublet, cell filter (doublets, top-1% UMI, `pct_mt<10`, `n_genes>100`), drop zero-count genes. QC violins/scatters. | local (notebook) |
| 3 | `01_3_normalization` | `scran_norm.py` | scran normalization via scib (`quickCluster → computeSumFactors → logNormCounts`); adds `size_factors`, clears `.raw`, casts `.X` to float32. | HPC (SLURM) |
| 3 | `01_3_normalization` | `submit_scran_norm.slurm` | SLURM wrapper for the scran step. | HPC (SLURM) |
| 4 | `01_4_cc_and_annotation` | `cell_cycle_score.py` | Tirosh/Regev cell-cycle scoring → `S_score`, `G2M_score`, `phase`. | terminal |
| 5 | `01_4_cc_and_annotation` | `celltypist_annotation.py` | CellTypist annotation (`Cells_Adult_Breast.pkl`, majority voting) → `cell_type`, `celltypist_predicted`. | terminal |
| 6 | `01_4_cc_and_annotation` | `fraction_reassignment.py` | Recode `fraction` CD45+/CD45- → `imm`/`non_imm` from the CellTypist lineage (in-place); `dataset_origin` left untouched. | terminal |
| 7 | `01_5_scib_pp` | `scib_reduce_data.py` | Batch-aware HVG (`cohort`, 2000) + PCA(50) + neighbors + UMAP; also writes the selected HVG list. | terminal |
| 8 | `01_5_scib_pp` | `scib_clustering.py` | Leiden optimal-resolution sweep (0.1–1.0) vs `cell_type` NMI → `shiao.h5ad`. | terminal |
| 9 | `01_6_visualization` | `visualization_unintegrated.ipynb` | Figures on the final unintegrated object: count barplots, annotation/sort-purity check, normalization, QC & Leiden UMAPs. Read-only. | local (notebook) |

> **Note on step 6.** `fraction_reassignment.py` must run *after* `celltypist_annotation.py`
> (it needs `cell_type`) and *before* `01_5` (so the recoded `fraction` propagates into
> `shiao.h5ad`). It rewrites the 01_4 annotated file in place.

## Data location (`DATA_DIR`)

As in phase 00, the heavy objects live **outside** the repo; every script resolves
`DATA_DIR` from the environment and never hardcodes a default. The repo holds only code
and lightweight QC figures (`01_pre_processing/figures/`).

The `.h5ad` chain under `$DATA_DIR` (each step consumes the previous file):

    $DATA_DIR/
    ├── all_samples_combined.h5ad                                  # input, from 00_6
    ├── all_samples_combined_scrublet.h5ad                         # 01_2  (QC + Scrublet + filters)
    ├── all_samples_combined_scrublet_norm.h5ad                    # 01_3  (scran log-norm)
    ├── all_samples_combined_scrublet_norm_cc.h5ad                 # 01_4  (+ cell cycle)
    ├── all_samples_combined_scrublet_norm_cc_annotated.h5ad       # 01_4  (+ CellTypist; fraction recoded in place)
    ├── all_samples_combined_scrublet_norm_cc_annotated_reduced.h5ad  # 01_5  (+ HVG/PCA/neighbors/UMAP)
    ├── shiao_hvg_2k_unintegrated.csv                              # 01_5  (selected HVG symbols)
    └── shiao.h5ad                                                 # 01_5  (+ Leiden) — definitive unintegrated object

Auxiliary inputs also expected under `$DATA_DIR`:

| File | Used by | Via env var |
|---|---|---|
| `regev_lab_cell_cycle_genes.txt` (97 genes: 43 S + 54 G2M) | `cell_cycle_score.py` | `CC_GENES` |
| `Cells_Adult_Breast.pkl` (CellTypist model, Kumar et al. 2023) | `celltypist_annotation.py` | — (loaded from `DATA_DIR`) |

**Local steps** — everything except 01_3 was run locally: the notebooks (01_1, 01_2, 01_6)
interactively, the `.py` scripts (01_4, 01_5) with `python3` after `export DATA_DIR=...`
(`cell_cycle_score.py` also needs `export CC_GENES=...`):

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
python3 01_4_cc_and_annotation/fraction_reassignment.py
```

**HPC step** — only `01_3` (scran) was run on the cluster via its SLURM wrapper (`DATA_DIR`
set at the top; submit from the folder that contains it, so the job anchors on
`SLURM_SUBMIT_DIR`):

```bash
cd 01_3_normalization && sbatch submit_scran_norm.slurm
```

## Object conventions (carried through every step)

- `.X` = **scran log1p-normalized** expression from 01_3 onward (already in log space — do
  **not** re-`log1p`). `.layers['counts']` = raw integer counts, kept unchanged throughout.
- `var_names` = gene symbols (Ensembl in `.var['gene_ids']`).
- Integration `batch_key = 'cohort'` (patient, technical batch); biological `label_key = 'cell_type'`.
- `dataset_origin` (`immune`/`non_immune`) = technical **CD45 sort**, never modified.
  `fraction` (`imm`/`non_imm`) = **biological** immune/non-immune class from CellTypist (recoded in 01_4).

## Key parameters (verbatim, ready for Materials & Methods)

**QC & filtering (01_2)**
- Scrublet per biopsy (`batch_key='sample'`), `expected_doublet_rate=0.10`,
  `n_prin_comps` retried `[30, 20, 10, 5]`, `random_state=0`.
- Doublet call: fixed cutoff `doublet_score > 0.25` (not Scrublet's per-batch auto threshold,
  unreliable on small/homogeneous biopsies); NaN scores → not a doublet.
- Presumed doublets: top **1% by `total_counts` per biopsy**.
- Cell filter (keep iff): not a Scrublet doublet **and** not top-1% UMI **and**
  `pct_counts_mt < 10` **and** `n_genes_by_counts > 100`.
- Gene filter: `sc.pp.filter_genes(min_cells=1)` (zero-count genes break scran).

**Normalization (01_3)**
- `scib.preprocessing.normalize`: `min_mean=0.1`, `log=True`, `precluster=True`
  (leiden `quickCluster`), `sparsify=False`; `SEED=0`.
- Output `.X` = scran log-normalized float32; `.obs['size_factors']`; `.raw` cleared
  (scib sets a full duplicate that nothing downstream reads).

**Cell cycle (01_4)**
- `sc.tl.score_genes_cell_cycle`, Tirosh/Regev 97-gene signature (43 S + 54 G2M),
  `random_state=0`. Phase = dominant score; G1 when both scores non-positive.

**Annotation (01_4)**
- CellTypist model `Cells_Adult_Breast.pkl` (Kumar et al. 2023 adult breast atlas, 58 labels),
  `majority_voting=True`, run on a **temporary** CP10K+log1p matrix built from raw counts
  (the main scran `.X` is never modified). `cell_type` = majority-voting label;
  `celltypist_predicted` = raw per-cell prediction.

**Fraction reassignment (01_4)**
- `fraction` recoded CD45+/CD45- → `imm`/`non_imm` by CellTypist lineage (immune vs
  non-immune label sets); `dataset_origin` preserved so the sort-vs-lineage mismatch stays
  inspectable. Written back atomically (temp file + `os.replace`).

**Feature selection + reduction (01_5)**
- `scib.preprocessing.reduce_data`: `batch_key='cohort'`, `flavor='cell_ranger'`,
  `n_top_genes=2000`, `n_bins=20`, PCA 50 comps (`svd_solver='arpack'`), neighbors + UMAP on
  `X_pca`. HVG selection is per-patient then merged (`overwrite_hvg=True`).

**Clustering (01_5)**
- `scib.clustering.cluster_optimal_resolution`, `label_key='cell_type'`, leiden
  (`flavor='igraph'`, `n_iterations=2`), resolution grid **0.1–1.0** (capped from scib's
  0.1–2.0), selected by max NMI; every per-resolution column
  (`optscib_unintegrated_leiden_<res>`) is kept for the 01_6 grids.

## Established numbers

- Final unintegrated object `shiao.h5ad`: **619,693 cells × 30,869 genes**.
- CD45 sort (`dataset_origin`): 489,577 immune (CD45+) + 130,116 non-immune (CD45-).
- CellTypist: **48 of 58** atlas labels observed.
- Sort-vs-lineage mismatch: **10.68%** of CD45+ cells (52,275 / 489,577) are annotated as a
  non-immune type; 4.44% of CD45- cells are annotated immune.
- After the biological recode: `fraction` = 443,083 `imm` / 176,610 `non_imm`.

## Critical methodological notes

- `dataset_origin` is the **library-level** CD45 magnetic-bead sort (from phase 00) and is
  never recomputed. `fraction` is the only field re-derived here (biological, from CellTypist).
- CellTypist runs on a throwaway CP10K+log1p copy: the training normalization it expects
  differs from scran, so scoring on scran `.X` would be miscalibrated.
- `01_6` is purely diagnostic (figures) and does not modify the data.
- In `Cells_Adult_Breast.pkl` the `Lymph-*` labels are **lymphatic endothelial** subtypes
  (not lymphocytes) and are therefore classified non-immune.
