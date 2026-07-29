# 00_raw_data_processing

This is the first phase of the thesis: from raw ENA/SRA FASTQ files to a single `.h5ad` anndata object comprehensive of the trial's metadata like cohort, treatment, response, ecc...
Dataset: Shiao et al., *Cancer Cell* 2024 (DOI 10.1016/j.ccell.2023.12.012),
BioProject **PRJNA1032700** / GEO **GSE246613**

## Execution order

**The whole phase runs on the cluster**, one SLURM job per step: the 2.4 TB of FASTQ is downloaded
straight onto the HPC storage and never touches the local machine, and every step after it reads
that storage. Each folder therefore holds the script that does the work plus the `.slurm` wrapper
that is the actual entry point.

| # | Folder | Script | What it does | Where |
|---|--------|--------|--------------|-------|
| 1 | `00_1_metadata` | `fetch_metadata.sh` | Download ENA metadata for the BioProject | HPC (SLURM, `submit_metadata.slurm`) |
| 2 | `00_2_sample_mapping` | `build_sample_map.sh` | SRR to patient/timepoint/fraction; exclude TCR/BCR and mouse | HPC (SLURM, `submit_sample_map.slurm`) |
| 3 | `00_3_download` | `download_fastq.sh` | Download FASTQ via `prefetch` + `fasterq-dump` | HPC (SLURM, `submit_download.slurm`; array, 1 task = 1 SRR) |
| 3b | `00_3_download` | `verify_integrity.sh` | Checks: `gzip -t`, R1/R2 lengths, 5' TSO | HPC (SLURM, `submit_verify_integrity.slurm`) |
| 4 | `00_4_rename` | `rename_for_cellranger.sh` | Rename to `*_S1_L001_R[12]_001.fastq.gz`; build `unique_samples.txt` | HPC (SLURM, `submit_rename.slurm`) |
| 5 | `00_5_cellranger` | `run_cellranger_batch.sh` | `cellranger count` in parallel (5x16 cores) | HPC (SLURM, `submit_cellranger.slurm`) |
| 6 | `00_6_build_h5ad` | `build_requested_metadata.py` + `build_combined_h5ad.py` | Build `sample_metadata_final.tsv`, then concatenate the `.h5` into one AnnData | HPC (SLURM, `submit_h5ad_concat.slurm`, both in order) |
| - | `utils` | `check_download_completeness.sh` | Check download completeness (R1 & R2) | HPC (shell, seconds — no allocation) |
| - | `utils` | `check_cellranger_completeness.sh` | Check for missing sample matrices | HPC (shell, seconds — no allocation) |
| - | `utils` | `fetch_missing_sample.sh` | Re-download a single corrupted/missing SRR | local (shell) — superseded on the cluster, see below |

The two `utils` completeness checks are `stat` loops over 150 paths, no allocation
needed. `fetch_missing_sample.sh` is the local way to repair a single run; on the cluster the same
job is done by resubmitting that run's index of the download array (below), so the FASTQ is
re-fetched inside a proper allocation.


## Data location (`DATA_DIR`)

The FASTQ and the Cell Ranger outputs (~2.4 TB) live **outside** the repo. Every script that touches
them requires `DATA_DIR` and exits immediately with an explanatory message if it is unset.

Expected layout under `$DATA_DIR` (the scripts create it as they go):

    $DATA_DIR/
    ├── fastq_raw/                 # step 3: SRR*_1/_2.fastq.gz as downloaded (~2.4 TB)
    ├── renamed_fastq/             # step 4: <sample>_S1_L001_R[12]_001.fastq.gz + unique_samples.txt
    ├── cellranger_out/            # step 5: <sample>/outs/filtered_feature_bc_matrix.h5
    └── all_samples_combined.h5ad  # step 6: final AnnData

The small TSVs (`ena_full_metadata.tsv`, `sample_map_gex_final.tsv`, `sample_metadata_final.tsv`)
stay in the repo, next to the script that writes them. Those scripts resolve their paths from their
own location, so they can be run from any working directory.

Steps 1, 2 and `build_requested_metadata.py` need no `DATA_DIR` at all: they read and write small
TSVs kept in the repo. They still run as cluster jobs, because their outputs are the inputs of the
steps that do touch the terabytes, and having them there means the repo on the cluster is always
self-consistent, no TSV copied by hand from a laptop.

Values used for the thesis run (Methods), these are the defaults written in the wrappers, as
`${VAR:-default}`, so a value passed at submission always wins:

| Variable | Value | Set in |
|---|---|---|
| `DATA_DIR` | `/users/genomics/albertoc/Tesi/hopes_and_dreams` | `submit_download.slurm`, `submit_verify_integrity.slurm`, `submit_rename.slurm`, `submit_cellranger.slurm`, `submit_h5ad_concat.slurm` |
| `REF` | `/users/genomics/albertoc/Tesi/cell_ranger/refdata-cellranger-GRCh38-3.0.0` | `submit_cellranger.slurm` |

The scripts themselves keep no default and abort if `DATA_DIR` is unset: a hardcoded path inside
them would only ever be correct on one machine and would otherwise fail late and cryptically, or
silently write terabytes into the clone. The wrappers are the one place where the site is named.


## Running the phase

Submit **from the folder containing the wrapper**: SLURM copies the script to a spool directory, so
`$0` cannot locate the repo and the jobs anchor on `SLURM_SUBMIT_DIR` instead.

```bash
(cd 00_1_metadata      && sbatch submit_metadata.slurm)
(cd 00_2_sample_mapping && sbatch submit_sample_map.slurm)
(cd 00_3_download      && sbatch submit_download.slurm)          # array, hours per task
(cd 00_3_download      && sbatch submit_verify_integrity.slurm)  # only when the array is done
(cd 00_4_rename        && sbatch submit_rename.slurm)
(cd 00_5_cellranger    && sbatch submit_cellranger.slurm)
(cd 00_6_build_h5ad    && sbatch submit_h5ad_concat.slurm)
```

To run somewhere else, override the paths at submission instead of editing the files:

```bash
(cd 00_3_download && sbatch --export=ALL,DATA_DIR=/other/root submit_download.slurm)
```

Steps 3b to 6 can be chained so each starts only if the previous one succeeded:

```bash
cd 00_3_download && dl=$(sbatch --parsable submit_download.slurm)
vf=$(sbatch --parsable --dependency=afterok:$dl submit_verify_integrity.slurm)
cd ../00_4_rename && rn=$(sbatch --parsable --dependency=afterok:$vf submit_rename.slurm)
cd ../00_5_cellranger && cr=$(sbatch --parsable --dependency=afterok:$rn submit_cellranger.slurm)
cd ../00_6_build_h5ad && sbatch --dependency=afterok:$cr submit_h5ad_concat.slurm
```

Chaining is convenient but unforgiving: `afterok` on the download array requires **all 150** tasks to
succeed, and task 44 (`SRR26541168`, the empty run) never will. Either drop that index
(`--array=1-43,45-150`) or submit the steps one at a time, checking each with `sacct` and the `utils`
scripts before moving on — which is how the thesis run was actually done.

Two steps reach the internet from a compute node: step 1 (ENA API) and step 3 (`prefetch` to NCBI).
If the compute nodes are firewalled, those two have to run where outbound access exists (the login
node for step 1, which is seconds of `wget`); everything else is storage-local.

Allocations requested by the wrappers (for Methods; all CPU-only, the cluster has no GPU). The
partition sets the wall-clock limit — `long` is 30 days — so the choice between `normal` and `long`
is the only thing to get right.

| Step | Partition | CPUs | Mem | Note |
|---|---|---|---|---|
| 1 metadata | `normal` | 1 | 2 G | one `wget` |
| 2 sample map | `normal` | 1 | 2 G | needs GNU awk (checked by the wrapper) |
| 3 download | `long` | 8 | 16 G | array `1-150%4`, hours per task |
| 3b verify | `long` | 2 | 8 G | `gzip -t` streams all 2.4 TB |
| 4 rename | `long` | 2 | 8 G | `rsync` copy, needs another ~2.4 TB free |
| 5 Cell Ranger | `long` | 90 | 470 G | whole node, 5 samples x 16 cores |
| 6 h5ad | `normal` | 8 | 150 G | 149 matrices held in memory at concat |


## Key parameters (verbatim, ready for Materials & Methods)

**Download**
- SRA Toolkit: `prefetch --max-size 100G` (the 20 GB default skips the 37-40 GB files)
  + `fasterq-dump --include-technical --split-files` — `--include-technical` is
  critical: without it the barcode read (R1) is dropped as a "technical read".
- Run as a SLURM array, one task per run (`--array=1-150%4`, the index being the line of
  `sample_map_gex_final.tsv`), 8 cores and 16 GB per task on partition `long`.
  Only 4 tasks download at a time: the step is network- and I/O-bound and NCBI
  throttles parallel transfers. Environment: `catalano_env` (from
  `environments/benchmark-hpc.yml`, which declares `sra-tools`); override with
  `DOWNLOAD_ENV`.
- The array is resumable: `download_fastq.sh` skips any SRR whose `_2.fastq.gz` is
  already present, so a task killed on wall clock (or a whole failed batch) is
  resubmitted with `sbatch --array=<idx>,<idx> submit_download.slurm` and picks up
  where it stopped. Task 44 (`SRR26541168`, the empty run) is expected to fail.
- A task interrupted mid-`pigz` leaves a truncated `.fastq.gz` that the skip test still
  accepts — this is what step 3b (`gzip -t`) catches afterwards. To repair, delete the
  two files of that run and resubmit its index; the index of a given accession is
  `grep -n <SRR> 00_2_sample_mapping/sample_map_gex_final.tsv | cut -d: -f1`.

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

