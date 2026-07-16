# 00_raw_data_processing

This is the first phase of the thesis: from raw ENA/SRA FASTQ files to a single `.h5ad` anndata object comprehensive of the trial's metadata like cohort, treatment, response, ecc...
Dataset: Shiao et al., *Cancer Cell* 2024 (DOI 10.1016/j.ccell.2023.12.012),
BioProject **PRJNA1032700** / GEO **GSE246613**

## Execution order

| # | Folder | Script | What it does | Where |
|---|--------|--------|--------------|-------|
| 1 | `00_1_metadata` | `fetch_metadata.sh` | Download ENA metadata (incl. `fastq_md5`) for the BioProject | terminal |
| 2 | `00_2_sample_mapping` | `build_sample_map.sh` | SRR to patient/timepoint/fraction; exclude TCR/BCR and mouse | terminal |
| 3 | `00_3_download` | `download_fastq.sh` | Download FASTQ via `prefetch` + `fasterq-dump` | terminal (I/O) |
| 3b | `00_3_download` | `verify_integrity.sh` | Checks: `gzip -t`, R1/R2 lengths, 5' TSO | terminal |
| 4 | `00_4_rename` | `rename_for_cellranger.sh` | Rename to `*_S1_L001_R[12]_001.fastq.gz`; build `unique_samples.txt` | terminal |
| 5 | `00_5_cellranger` | `run_cellranger_batch.sh` | `cellranger count` in parallel (5x16 cores) | HPC |
| 5 | `00_5_cellranger` | `submit_cellranger.slurm` | SLURM wrapper for the Cell Ranger batch | HPC |
| 6 | `00_6_build_h5ad` | `build_requested_metadata.py` | Build `sample_metadata_final.tsv` | terminal |
| 6 | `00_6_build_h5ad` | `build_combined_h5ad.py` | Concatenate the `.h5` files into one AnnData | HPC |
| 6 | `00_6_build_h5ad` | `submit_h5ad_concat.slurm` | SLURM wrapper for the concatenation | HPC |
| - | `utils` | `check_download_completeness.sh` | Check download completeness (R1 & R2) | terminal |
| - | `utils` | `fetch_missing_sample.sh` | Re-download a single corrupted/missing SRR | terminal |
| - | `utils` | `check_cellranger_completeness.sh` | Check for missing sample matrices | terminal |


## Data location (`DATA_DIR`)

The FASTQ and the Cell Ranger outputs (~2.4 TB) live **outside** the repo. Every script that touches
them requires `DATA_DIR` and exits immediately with an explanatory message if it is unset. This is
deliberate: a hardcoded default would only ever be correct on one machine, and would otherwise fail
late and cryptically, or silently write terabytes into the clone.

Expected layout under `$DATA_DIR` (the scripts create it as they go):

    $DATA_DIR/
    ├── fastq_raw/                 # step 3: SRR*_1/_2.fastq.gz as downloaded (~2.4 TB)
    ├── renamed_fastq/             # step 4: <sample>_S1_L001_R[12]_001.fastq.gz + unique_samples.txt
    ├── cellranger_out/            # step 5: <sample>/outs/filtered_feature_bc_matrix.h5
    └── all_samples_combined.h5ad  # step 6: final AnnData

The small TSVs (`ena_full_metadata.tsv`, `sample_map_gex_final.tsv`, `sample_metadata_final.tsv`)
stay in the repo, next to the script that writes them. Those scripts resolve their paths from their
own location, so they can be run from any working directory.

**Terminal steps** — export `DATA_DIR` once:

```bash
export DATA_DIR=/path/with/2.4TB/free
bash 00_3_download/download_fastq.sh
```

**HPC steps** — `DATA_DIR` and `REF` are set at the top of the two SLURM wrappers. To replicate,
edit those lines, then submit **from the folder containing the wrapper**: SLURM copies the script to
a spool directory, so `$0` cannot locate the repo and the job anchors on `SLURM_SUBMIT_DIR` instead.

```bash
cd 00_5_cellranger && sbatch submit_cellranger.slurm
```

Values used for the thesis run (recorded here for Methods; intentionally *not* defaults in the code):

| Variable | Value |
|---|---|
| `DATA_DIR` | `/users/genomics/albertoc/Tesi/hopes_and_dreams` |
| `REF` | `/users/genomics/albertoc/Tesi/cell_ranger/refdata-cellranger-GRCh38-3.0.0` |


## Key parameters (verbatim, ready for Materials & Methods)

**Download**
- SRA Toolkit: `prefetch --max-size 100G` (the 20 GB default skips the 37-40 GB files)
  + `fasterq-dump --include-technical --split-files` — `--include-technical` is
  critical: without it the barcode read (R1) is dropped as a "technical read".

**Cell Ranger**
- Version 4.0.0 (via `module load cellranger/4.0.0`).
- Reference refdata-cellranger-GRCh38-3.0.0 (same as the paper).
- `--r1-length=26` (R1 = 16 bp barcode + 10 bp UMI + padding to 101/151 bp).
- `--chemistry=SC5P-R2` (10x Chromium 5' v1, confirmed by the TSO `AAGCAGTGGTATCAAC`).
- `--localcores=16 --localmem=64`, 5 samples in parallel (`xargs -P 5`) on a 90c/470G node.

**AnnData construction**
- Input: `filtered_feature_bc_matrix.h5` (Cell Ranger cell-calling already applied).
- `var_names` = gene symbols made unique (`var_names_make_unique()`); Ensembl in `.var['gene_ids']`.
- Raw counts in `.layers['counts']`.
- `ad.concat(join='outer', merge='same', index_unique=None)`, barcodes made unique per sample.
- Metadata: `srr, sample_id, batch, cohort, treatment, fraction, dataset_origin, response`.

## Established numbers

- 150 mapped human GEX runs to 149 processable (P30_C_P / SRR26541168 excluded, empty).
- 104 CD45+ (P, immune) + 46 CD45- (N, non_immune).
- Only 16/34 patients have CD45- libraries deposited.
- Response (author-provided, GSE246613): 9 R1 / 14 R2 / 11 NR across 34 patients.

## Critical methodological notes

- `dataset_origin` (immune/non_immune) is library-level, obtained using magnetic beads.
- `response` is patient-level, not recomputed: transcribed from GSE246613.

