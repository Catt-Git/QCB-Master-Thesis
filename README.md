# QCB---Master-Thesis

Repository with the full code used and curated for the analysis carried out during my traineeship abroad at the Hospital del Mar Research Institute (IMIM).

Dataset: Shiao et al., *Cancer Cell* 2024 (DOI 10.1016/j.ccell.2023.12.012), BioProject
**PRJNA1032700** / GEO **GSE246613** - 149 human 10x 5' GEX libraries from 34 breast cancer
patients, from FASTQ to an integration benchmark and a DRVI interpretation of the non-immune
and epithelial compartments.

## Reproducing the analysis

**1. Clone anywhere.** No script contains the repository's own path: the `.sh`/`.py`/`.R` files
locate themselves through `SCRIPT_DIR`, the SLURM wrappers through `SLURM_SUBMIT_DIR`. The clone
can be `~/Desktop/QCB-Master-Thesis` locally and `~/Tesi/QCB-Master-Thesis` on the cluster, and
nothing needs editing.

**2. Create the environment** for the machine you are on - `benchmark-py-r` locally,
`benchmark-hpc` on the cluster, plus `scgen-py` for the scGen method only. See
[environments/README.md](environments/README.md).
The SLURM wrappers activate this cluster environment under the name it was built with on the HPC
used for the thesis, `catalano_env`; create it under that name, or override it per phase
(`DOWNLOAD_ENV`, `PREPROC_ENV`, `METRICS_ENV`).

**3. Set `DATA_DIR`.** This is the only thing to configure, and it is **one value for the whole
pipeline** - every phase reads and writes the same root, in the per-phase subfolders documented
in [datasets/README.md](datasets/README.md).

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets                  # local
export DATA_DIR=/users/genomics/albertoc/Tesi/hopes_and_dreams/datasets   # cluster
```

On the cluster the root sits outside the clone (2.4 TB does not belong in a git working tree,
and the HPC home is the quota'd filesystem); symlink it into the clone if you want the export
line to be literally identical on both machines. The scripts never default this path: they abort
with an explanatory message if `DATA_DIR` is unset, so the terabytes can never land in the clone
by accident. The SLURM wrappers do carry the thesis-run value as `${DATA_DIR:-default}`, which a
value passed at submission always overrides:

```bash
sbatch --export=ALL,DATA_DIR=/other/root submit_download.slurm
```

**4. Run the phases in order.** Each folder's README is the authority on its own steps, on where
they run (local vs SLURM) and on the parameters used for the thesis run.

| Phase | What it produces | Where |
|---|---|---|
| [00_raw_data_processing](00_raw_data_processing/) | FASTQ → Cell Ranger → `all_samples_combined.h5ad` | cluster (SLURM) |
| [01_pre_processing](01_pre_processing/) | QC, Scrublet, scran, cell cycle, CellTypist, HVG/PCA/UMAP → `shiao.h5ad` | mixed |
| [02_integration_benchmark](02_integration_benchmark/) | 10 integration methods × 13 metrics | integration local, metrics on SLURM |
| [03_drvi_non_immune](03_drvi_non_immune/) | non-immune subset + DRVI latent space | local |
| [04_drvi_epithelial](04_drvi_epithelial/) | epithelial subset (one lineage deeper than 03) | local |
| [05_drvi_tumoral_epi](05_drvi_tumoral_epi/) | inferCNV → `malignant` / `non_malignant`, CellTypist re-run on the non-malignant, then 04's procedure redone on that basis | local |
| [06_epi_treatment](06_epi_treatment/) | the epithelial subset cut by treatment (BASE / PD1 / RTPD1) + one DRVI run each | local (or a SLURM array) |

Phase 00 is skippable if the counts are obtained another way: everything downstream needs only
`$DATA_DIR/all_samples_combined.h5ad`, and phases 01→06 chain from there. 03, 04 and 05 are three
independent branches off `shiao.h5ad`: none of them reads what another wrote. 04 and 05 apply the
same procedure to two different definitions of the cell set - 04 to the epithelial compartment as
CellTypist labelled it, 05 to the same compartment after inferCNV has separated the malignant
cells - so their step numbering runs in parallel and their results are meant to be read side by
side. 06 is the one branch that continues another: it splits 05's object by treatment.

**What is not in git.** Everything under `datasets/` except three small non-regenerable inputs
(the CellTypist model, the Tirosh/Regev cell cycle list, and that README), plus the run logs.
Every other artefact is rebuilt by re-running the phase that writes it.
