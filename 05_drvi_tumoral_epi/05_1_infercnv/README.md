# 05_1_infercnv

First step of phase 05, and the one the whole phase rests on. It answers the question phase 04
states in its own README that it cannot answer: **which epithelial cells are malignant**.

Phase 04 is not modified by any of this. It stays a frozen branch off `shiao.h5ad`, with its
results and its honest caveat about what it could not know; phase 05 is a second, independent
branch off the same object that starts here.

## Why

`cell_type` comes from CellTypist with `Cells_Adult_Breast.pkl` (Kumar et al. 2023), a model
trained on the **normal** adult breast atlas. It has no malignant class. A TNBC cell is
therefore not left unassigned - it is assigned to whichever normal state it resembles most, and
it then votes in the majority-voting step of every neighbourhood it sits in. Two consequences,
and 04_1 - 04_9 inherit both:

1. the epithelial subset mixes normal and malignant epithelium under labels that cannot tell
   them apart, so every "epithelial state" reported downstream is a state of that mixture;
2. because voting is neighbourhood-based, the normal cells' labels are contaminated too.

inferCNV infers copy-number variation from the smoothed expression of each cell against a
diploid reference. Aneuploidy is the thing that actually separates carcinoma cells from the
normal epithelium around them, and it is measurable from this data without any new experiment.

Once the call exists, both problems are fixed in one move: the malignant cells collapse to a
single label `malignant`, and the 01_4 annotation is re-run on what is left - a normal-breast
model against a normal-breast population, which is the regime it was trained for.

## The design, and the two choices it rests on

**The reference is immune, and it is per patient.** This dataset has no normal or adjacent
tissue: `treatment` is BASE / PD1 / RTPD1, three timepoints of the same tumour, in all 34
cohorts. The diploid baseline therefore has to come from inside each tumour, and the standard
choice applies - T/NK and myeloid cells, which are here in quantity (never below 1,272 per
cohort). Two reference *groups* rather than one pooled group, so that inferCNV takes the
residual against the bounds of the per-group means and a gene that is simply higher in myeloid
than in T cells cannot masquerade as a gain. B and plasma cells are in neither: 62,074 plasma
cells with a clonally skewed immunoglobulin transcriptome would put spurious structure on
chr2, chr14 and chr22.

One run per patient, never pooled - the residual is defined against the reference cells present
in the run, so a pooled run would compare patient A's epithelium against patient B's immune
cells and read the batch difference as copy number. Same reasoning as 04_4 scoring CytoTRACE2
per patient.

**The stromal block is the internal control.** Fibroblasts and endothelium are carried into
every run as *observations*, not as reference. They go through identical smoothing and
denoising, they are not the malignant compartment of a carcinoma, and they are therefore a free
null distribution - one per cohort, on the same scale as that cohort's epithelium. The
thresholds in `call_malignant.ipynb` are quantiles of that null.

That leaves the immune reference cells unused by the rule, which makes their crossing rate a
genuine specificity check rather than a tautology. It is printed next to the epithelial rate.

## Execution order

| # | File | What it does | Where |
|---|------|--------------|-------|
| 1 | `prepare_infercnv_input.py` | `shiao.h5ad` → per-cohort sparse counts + `annotations.tsv`; downloads the hg38/GENCODE v27 gene ordering file once | local, `benchmark-py-r` |
| 2 | `run_infercnv.R` | `infercnv::run()` on one cohort → a per-cell table (`cnv_score`, `cnv_corr`, per-chromosome means) + the heatmap; deletes the working directory | local, `infercnv-r` |
| 3 | `infercnv_all.sh` | Driver for 1 and 2: prepares everything, then loops step 2 over every prepared cohort, resuming, logging to `logs/` | local |
| 4 | `call_malignant.ipynb` | The decision. Per-cohort thresholds off the stromal null, specificity check, figures → `cnv_status.csv` | local (notebook) |
| 5 | `recelltypist_nonmalignant.py` | CellTypist re-run on the non-malignant cells only → `cell_annotation_cnv.csv` | local, `benchmark-py-r` |

```bash
conda env create -f ../../environments/infercnv-r.yml     # once
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
./infercnv_all.sh --threads 12                            # steps 1-3, ~2 min per small cohort
#   then run call_malignant.ipynb
conda activate benchmark-py-r && python3 recelltypist_nonmalignant.py
```

## What it writes, and what it does not

Everything lands in **`$DATA_DIR/05_tum/`**. Nothing under `$DATA_DIR/04_epi/` is read or
written, and `shiao.h5ad` is only ever opened read-only and backed, so every result phase 04
already has is untouched by construction rather than by care.

```
$DATA_DIR/05_tum/
├── gene_order_hg38_gencode_v27.txt   downloaded once
├── cohort_census.csv                 what went into each run
├── input/<cohort>/                   counts.mtx, genes.tsv, barcodes.tsv, annotations.tsv
├── work/<cohort>/                    inferCNV's working directory - DELETED on success
├── summary/<cohort>_cnv.csv          one row per cell: group, cnv_score, cnv_corr, chr1..chr22
├── cnv_status.csv                    the call, one row per cell of shiao.h5ad
└── cell_annotation_cnv.csv           the re-annotation, one row per cell of shiao.h5ad
```

The deliverable is **two per-cell tables, not a new `.h5ad`**. A CNV call is metadata; writing a
second 4 GB copy of `shiao.h5ad` to carry two columns would cost more disk than this machine has
(see below) and would fork the object every later step has to choose between. 05_2 joins
`cell_annotation_cnv.csv` onto the object it builds.

Figures go to `../figures/05_1_infercnv/`, one `infercnv_<cohort>.png` heatmap per cohort plus
the notebook's; the per-cohort thresholds table goes to `../tables/05_1_infercnv/`.

> **Disk.** inferCNV's working directory is 1-3 GB *per cohort* and this machine had ~29 GB free
> when 05_1 was written - 33 cohorts do not fit. `run_infercnv.R` therefore reduces each run to
> its per-cell table in memory and then deletes the working directory, and `infercnv_all.sh`
> deletes it again after a cohort that failed. `--keep-work` exists for debugging one cohort;
> do not pass it to a full run. This is also why the counts are handed to R as Matrix Market
> rather than as the dense tab-delimited matrix inferCNV documents: 20 MB - 260 MB per cohort
> instead of ~1 GB.

## The three values of `cnv_status`

| value | which cells | meaning |
|---|---|---|
| `malignant` | epithelial cells above the cohort's stromal null on **both** axes | aneuploid |
| `non_malignant` | everything a run covered and did not call | diploid as far as this test can tell |
| `not_tested` | B/plasma cells, and every cell of a cohort no run covered | **no evidence either way** |

`not_tested` is deliberately not folded into `non_malignant`. For an immune cell the distinction
is academic; for an *epithelial* cell it is exactly the case this step exists to stop treating as
normal. One cohort is skipped upfront for having too little epithelium to define anything
(`MIN_EPI_CELLS = 50`; Patient06 has 14 epithelial cells), and its epithelium lands here.

## What the rest of phase 05 does with this

The call is metadata, and this step stops there. `05_2_subsetting` is what turns it into an
object: it joins `cell_annotation_cnv.csv` onto `shiao.h5ad` and re-runs the phase-01
pre-processing on the subset, the way `04_1_subsetting` did for the epithelial compartment.

One thing that does **not** carry over from 04 and has to be decided in `05_2`: `04_1` picks the
leiden resolution by maximising NMI against `cell_type`. On a malignant-only subset the new
`cell_type_cnv` is the constant `malignant`, so that criterion has no target left. The 01_4
labels the malignant cells were carrying are not constant - 18,943 of them are `Lumsec-prol`,
12,028 `Lumsec-basal`, 3,975 `LummHR-SCGB`, 1,174 `Lumsec-KIT` - so they stay usable as the NMI
target even though they are the contaminated labels this phase exists to replace. That is a
choice `05_2` has to make explicitly and write down, not inherit.

Phase 04's own `CAVEAT` machinery is untouched and stays true of phase 04: those tables and
figures really were produced without a CNV call. Phase 05 needs its own wording, not a patch to
04's.

## How to read `cnv_score` (it is smaller than it looks)

`cnv_score` has a median of ~0.002 on the malignant cells and ~0.0005 on the immune
reference, which reads as tiny until you remember it is a **mean of squares**. Its square
root is the quantity with units:

| | `cnv_score` | RMS residual |
|---|---|---|
| epithelial | 0.0019 | **4.4 %** |
| stromal (control) | 0.0009 | 3.0 % |
| T/NK reference | 0.0007 | 2.7 % |
| myeloid reference | 0.0004 | 2.0 % |

A 4-5 % RMS deviation from the diploid baseline is not marginal on inferCNV's own scale. Its
`expr.data` on this dataset runs 0.79 to 1.46, with the 1st and 99th percentiles at 0.913 and
1.112 - i.e. **±11 %** - and `plot_cnv` auto-thresholds the heatmap at 0.888-1.112. Measured on
Patient64 by keeping the residual matrix in memory: epithelial RMS 0.0482, stromal 0.0304,
reference 0.0228, and 0.0482² = 0.00232, exactly the reported median.

The number looks small because it averages over the WHOLE genome, most of which is not
altered. Per cell, the count of chromosomes whose mean residual exceeds 0.02 is a median of
**7 of 22** for epithelial cells against 2 for the reference; the fifteen flat ones drag the
mean towards zero. On Patient64 the per-chromosome medians of the malignant cells are chr9
+0.055, chr17 -0.059, chr7 +0.042, chr18 -0.035, chr13 -0.030 - and those are already
chromosome-wide averages, so individual genes inside the altered arms deviate considerably
more. It is the same profile the heatmap shows by eye.

> **Do not compare this number with a published one.** inferCNV neither defines nor outputs a
> "CNV score"; the statistic comes from Puram et al. 2017 and Neftel et al. 2019, where it is
> computed on a **log2** scale rather than on the modified expression centred at 1 that
> inferCNV returns (a factor of ~1/ln2 on the residual, so ~2.1x on the square). Its magnitude
> also depends on the smoothing window, on `cutoff`, and above all on the denoising - here
> `sd_amplifier = 1.5`, which flattened everything between 0.957 and 1.044. What is
> interpretable is the contrast WITHIN a run: 4.4 % against 2.7 %, i.e. 2.1x in RMS and 4.5x
> in squared units. The second axis, `cnv_corr`, has no scale by construction and does not
> have this problem at all.

## The i6 HMM

`run_infercnv.R --hmm` turns on inferCNV's i6 HMM, which converts residuals into discrete CNV
states and is why `jags` and `r-rjags` are in the environment. It is **off by default**: it costs
hours per cohort, and the binary call this step produces is read off the two continuous axes, as
in Puram et al. 2017 and Neftel et al. 2019. Turn it on for a cohort whose heatmap is worth
reading in detail, not for the full sweep.
