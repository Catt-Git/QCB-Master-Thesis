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

That leaves the immune reference cells unused by the rule, and the first version of this
README claimed their crossing rate was therefore a genuine specificity check. **That was
wrong.** `cnv_score` is a residual taken against the mean of the reference cells, so they sit
at ~0 by arithmetic, not on the method's merit, and their low crossing rate measures nothing.
The fix is the `immune_heldout` group described in "The revised architecture" below: immune
cells drawn before the reference and removed from its pool, which centre nothing and set no
threshold. The rate is still printed next to the epithelial one, but the number to read is
the held-out group's.

## Execution order

| # | File | What it does | Where |
|---|------|--------------|-------|
| 1 | `prepare_infercnv_input.py` | `shiao.h5ad` → per-cohort sparse counts + `annotations.tsv`; downloads the hg38/GENCODE v27 gene ordering file once | local, `benchmark-py-r` |
| 2 | `run_infercnv.R` | `infercnv::run()` on one cohort → a per-cell table (`cnv_score`, `cnv_corr`, per-chromosome means) + the heatmap; deletes the working directory | local, `infercnv-r` |
| 3 | `infercnv_all.sh` | Driver for 1 and 2: prepares everything, then loops step 2 over every prepared cohort, resuming, logging to `logs/` | local |
| 4 | `call_malignant.ipynb` | The decision. Per-cohort thresholds off the stromal null, specificity check, figures → `cnv_status.csv` | local (notebook) |
| 5 | `recelltypist_nonmalignant.py` | CellTypist re-run on the non-malignant cells only → `cell_annotation_cnv.csv` | local, `benchmark-py-r` |
| S1 | `sensitivity_resolution.sh` + `.py` | 4 leiden resolutions × 4 cohorts → does the verdict depend on granularity? | local, both envs |
| S2 | `sensitivity_reference.sh` + `.py` + `.R` | 4 reference configurations × 2 cohorts → does it depend on the baseline? | local, both envs |

The two sensitivity steps are not optional extras: they are what the architecture below
rests on, and each answers one objection that cannot be answered by argument. They write to
`../tables/05_1_infercnv/sensitivity/` and never touch `input/` or `summary/`.

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

---

# The `newcnv` call

A second, independent version of step 4. It reads the **same** `summary/*_cnv.csv` tables -
inferCNV is not re-run, nothing under `input/` or `summary/` changes - and writes a parallel
set of outputs tagged `newcnv`. `call_malignant.ipynb`, `cnv_status.csv` and
`cell_annotation_cnv.csv` are untouched, so the two calls can be carried side by side and
compared rather than one replacing the other.

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   run call_malignant_newcnv.ipynb          (minutes; no inferCNV, no GPU)
conda activate benchmark-py-r && python3 recelltypist_nonmalignant.py --tag newcnv
```

| | default | `newcnv` |
|---|---|---|
| notebook | `call_malignant.ipynb` | `call_malignant_newcnv.ipynb` |
| call | `cnv_status.csv` | `cnv_status_newcnv.csv` |
| re-annotation | `cell_annotation_cnv.csv` | `cell_annotation_newcnv.csv` |
| figures / tables | `05_1_infercnv/` | `05_1_infercnv/newcnv/` |

## Why

Three things were measured on the existing call and each is fixed by one change below.

**1. The verdict was taken per cell; Shiao et al. take it per epithelial subcluster.** Their
STAR Methods compute an i6 HMM state per gene per cell, then the median % of altered genes per
*subcluster*, then a z-score across subclusters. The HMM is how they produce the per-cell
statistic; the robustness comes from the median, not from the HMM. `newcnv` therefore keeps
the two continuous Puram/Neftel axes and adds the aggregation level, using
`optscib_epi_leiden` from `04_epi/shiao_epi.h5ad` - an epithelial clustering computed from
expression, frozen on disk, and never exposed to an inferCNV output, so using it here is not
circular. Measured cost of not doing this: **11,271 epithelial cells** sit in clusters that are
aneuploid at the median but fall below the per-cell threshold on one axis.

**2. Seven cohorts had no detectable tumour and were called anyway.** Comparing each cohort's
epithelial crossing rate against its own *stromal* crossing rate splits the dataset at a clean
gap: 4.0-6.2% for seven cohorts - what the null produces on itself - against 11.6-95% for the
rest. Six are dropped downstream by `MIN_CELLS_PER_COHORT`; Patient52 is not, and enters the
`tum` subset with 9,243 epithelial cells, 385 of them called malignant. `newcnv` adds a
per-cohort gate and those cells become `not_tested`.

**3. The stromal null was used as one block.** It holds fibroblasts, endothelium and
perivascular cells, and which of the three dominates changes the call in five cohorts (worst:
Patient10, 3.6% malignant off the fibroblast null against 48.8% off the endothelial one). The
operating cut still uses the pooled block - that is the highest-power estimate - but the three
are now computed separately, written to `stromal_subnull_sensitivity_newcnv.csv`, and cohorts
whose call moves by more than 15 points are flagged.

## What deliberately does not change

Per-patient runs, the immune reference (T/NK + myeloid), B/plasma excluded from every run, the
two continuous axes, no HMM.

The **reference composition** is the one axis kept further from Shiao than it could be. They
use every non-epithelial cell type as reference; here fibroblasts, endothelium and
perivascular cells stay out of it, as observations, because a control folded into the
reference can no longer detect anything - and it is that control which produced all three
findings above. B/plasma stay out of both: adding them would mean re-running 33 cohorts to
gain a null that is knowingly bad, since a clonally skewed immunoglobulin transcriptome puts
structure on chr2/chr14/chr22.

## The four values of `cnv_status`

| value | which cells | meaning |
|---|---|---|
| `malignant` | in a cluster aneuploid at the median, in a cohort that passed the gate | aneuploid |
| `borderline` | the two criteria disagree | evidence both ways |
| `non_malignant` | below the null on both criteria | diploid as far as this test can tell |
| `not_tested` | B/plasma, uncovered cohorts, **and every epithelial cell of a gated-out cohort** | no evidence either way |

`borderline` is a disagreement between two criteria, not a third quantile chosen by hand. It
is small - 302 cells - and its size is the calibration: a `borderline` set that grew would say
`MAL_MED` / `NORM_MED` are badly placed.

Both `borderline` and `not_tested` are excluded from the tumour object by the selection
`05_2_subsetting` already applies (`cnv_status == 'malignant'`), so neither needs new
downstream code. Both still go through `recelltypist_nonmalignant.py` and come back with a
normal-breast label, for the reason that script's docstring gives - restricting its input
would change the neighbourhood structure its majority vote is computed over - and `cnv_status`
travels next to that label so 05_2 reads the status rather than the label.

## What it produces

| | `cnv_status.csv` | `cnv_status_newcnv.csv` |
|---|---|---|
| malignant | 37,014 | **46,920** |
| borderline | - | 302 |
| non_malignant | 38,424 | 15,463 |
| not_tested (epithelial only) | 0 | 12,753 |
| cohorts with >= 200 malignant | 19 | **20** |
| cells in the `tum` subset | 36,192 | **46,030** |

Patient10 enters (396 cells), Patient52 leaves. The 997 epithelial cells of Patient01,
Patient04, Patient23 and Patient30 have no cluster label - 04_1 drops cohorts missing a
treatment timepoint and 05 does not - and fall back to the per-cell rule alone; the column
`epi_cluster` is empty for them, so they are identifiable.

## What to look at before trusting it

Two things this call surfaces and does not resolve.

* **Cluster 0** is called malignant, holds 16,478 cells and spans 29 cohorts. Shared arm-level
  events in TNBC make a large cross-patient aneuploid cluster possible, but it is the one
  cluster whose per-cohort heatmaps are worth opening by eye.
* **Patient26** has immune *reference* cells crossing its own cut at 6.3%, ten times any other
  cohort, against 1.2% for its stroma. It contributes 4,876 malignant cells. The reference is
  not used by the rule, so this is a genuine specificity failure in that run, not a tautology.
  The notebook prints it and `cohort_gate_newcnv.csv` carries `ref_rate` per cohort, so it
  cannot pass unnoticed; nothing acts on it automatically.


---

# The revised architecture (current)

Everything above describes how the step was first built. Three things changed after the
measurements in `sensitivity_*`, and this section is what the code now does. The
`newcnv` call's structure is unchanged in outline - two continuous axes, a per-cohort
stromal null, a cohort gate - but what the verdict is aggregated over, and what the run
contains, are different.

```
                 PATIENT
                    |
        +-----------+-----------+
        |                       |
   T/NK + myeloid           stromal (fibro / endo / PVL, stratified)
   REFERENCE                NULL              + 500 held-out immune (specificity)
        |                       |
        +-----------+-----------+
                    v
                inferCNV   (analysis_mode="subclusters", leiden 0.005)
                    |
                    v
             epithelial cells
                    |
                    v
          patient-specific subclusters
                    |
                    v
          median nullpos per subcluster
                    |
                    v
         compared with the stromal null      -> MAL_MED / NORM_MED
                    |
                    v
             cohort gate, 10%
                    |
        +-----------+-----------+
        v                       v
   CNV-positive            CNV-negative
   epithelium              epithelium (-> not_tested)
```

## What each run contains, and the two things that changed

**The reference is unchanged**: `ref_tcell` + `ref_myeloid`, 1,000 each, never merged.
`sensitivity_reference.py` establishes that this is not a lever - see below.

**The stromal null is now stratified and larger.** Sampled uniformly over its union, the
block followed whatever mix the cohort happened to have (Patient30 came out 81% fibroblast,
Patient64 60% endothelial), so the cut was a quantile of a different population in every
run - and the three do not sit at the same level. It is now drawn up to `cap/3` from each of
fibroblast / endothelium / perivascular, with the remainder redistributed when a stratum is
short, and the cap is 2,000 rather than 1,000. Patient53, which used to yield ~77%
fibroblast, now yields 668 / 666 / 666.

**`immune_heldout` is new, and it is the specificity control.** 500 immune cells drawn
BEFORE the reference and removed from its pool, carried as observations. The reference
cannot serve this purpose: `cnv_score` is a residual against the mean of the reference
cells, so those sit at ~0 as a matter of arithmetic and their crossing rate proves nothing.
These cells set no threshold and centre nothing. Across every sensitivity run - eight
reference configurations and sixteen resolution runs, spanning cohorts from 0% to 95%
tumour content - they were called **0.00%** of the time.

A fourth, smaller change: cells whose raw CellTypist prediction is epithelial are excluded
from the reference and from the held-out group (1.6-4.5% of immune cells, depending on the
type). The filter is deliberately narrow - requiring the voted and raw labels to agree would
discard 60% of some cohorts' immune cells for within-immune churn that has nothing to do
with tumour contamination.

## The unit the verdict is taken on

The call is aggregated over subclusters rather than taken per cell, and those subclusters
are computed by inferCNV itself, inside each per-patient run. Two things motivate this. The
first is robustness: the per-cell call has no gap to find - 40% of the epithelium sits
within 1.5x of its own cohort's cut - so where the line falls decides the fate of thousands
of cells, whereas a median over the cells of a subcluster is stable. This is also how Shiao
et al. take the verdict; the i6 HMM in their methods is only how they produce the per-cell
statistic being averaged, so the aggregation is adopted here without the HMM. The second is
that the previous aggregation used the epithelial leiden of phase 04, computed on the
**integrated** object: its cluster 0 held 16,478 cells drawn from 29 of the 33 cohorts. A
copy-number profile is private to a patient, so a cross-patient cluster cannot be a clone,
and "this cluster is aneuploid at the median" could not mean what it appeared to mean. A
partition computed inside a single run cannot have that defect.

**The circularity, stated rather than hidden.** Phase 04's clustering was chosen precisely
because nothing in 04 has ever seen an inferCNV output, which made it non-circular by
construction. The per-run subclusters are the opposite: they are computed *on the residual
matrix whose values are then aggregated over them*. A partition will therefore exist, and
will separate high from low, even in a patient with no aneuploidy at all - leiden partitions
noise as readily as signal. What keeps the verdict from being circular is where the
threshold comes from: the cut is a quantile of the stromal null, and the subclustering never
sees the stroma's relationship to the epithelium, nor the reference. **The subcluster
decides which cells are judged together; the null decides what counts as high.**

Two consequences follow, and neither is optional. The **cohort gate remains necessary** - it
is the only thing standing between "leiden found a highest subcluster" and "this patient has
a tumour". Patient52 is the case it exists for and the demonstration that the concern is
real: 9,243 epithelial cells, epithelium indistinguishable from its own stromal null
(AUC 0.50), and leiden still produced **146 subclusters**, one of which crossed. Without the
gate that one subcluster would have become a clone. With it, 1.35% of the epithelium is far
below `GATE_FRAC` and the patient is correctly CNV-negative. And subcluster ids are **per
run**: `epi_s1` exists in every cohort and denotes different cells in each, so the key is
always `(cohort, subcluster)`.

`MIN_SUB = 20` covers the rest: a subcluster below it gets no verdict of its own, because a
median over a handful of cells is not the aggregate this design rests on. Its cells fall
back to the per-cell rule and `subcluster_call` records that this happened. At the operating
resolution that is ~2% of the epithelium.

## Sensitivity 1: does the call depend on the clustering granularity?

`sensitivity_resolution.sh` sweeps 0.005 / 0.01 / 0.02 / 0.05 on four cohorts chosen to span
the range: Patient16 (tumour-rich), Patient43 (intermediate), Patient64 (thin stromal null),
Patient52 (no detectable tumour). Sixteen runs.

| cohort | subclusters, 0.005 → 0.05 | fraction of epithelium CNV-positive | verdict |
|---|---|---|---|
| Patient16 | 113 → 745 | 0.957 - 0.958 (spread **0.001**) | CNV-positive throughout |
| Patient43 | 53 → 377 | 0.513 - 0.548 (spread 0.035) | CNV-positive throughout |
| Patient64 | 3 → 16 | 0.902 - 0.902 (spread **0**) | CNV-positive throughout |
| Patient52 | 146 → 974 | 0.013 - 0.021 (spread 0.007) | **CNV-negative throughout** |

Over a tenfold range of resolution the number of subclusters grows about sevenfold and the
per-patient verdict does not change once. The reason is structural: the final aggregation
counts **cells**, not subclusters, so splitting the same cells more finely leaves each piece
on the same side of the threshold.

What does degrade is the aggregation itself. The percentage of epithelial cells sitting in
subclusters of at least 20 cells falls from 98-100% at 0.005 to 23-57% at 0.05, where the
"median per subcluster" has quietly become a per-cell call again. **That, and not the
verdict, is why the operating resolution is 0.005** - the coarse end of a range over which
the answer is invariant.

## Sensitivity 2: does the call depend on the baseline population?

`sensitivity_reference.sh` runs four reference configurations on the same observations and
the same gene set, so only the baseline differs: the current T/NK + myeloid, T + B, five to
six separate immune subtypes, and a single pooled immune group. This matters because
inferCNV takes the residual against the **bounds** of the per-group means, so anything an
observed cell has inside that band becomes exactly zero - the number of reference groups is
a conservativeness dial, and the four span it from no protection to the most the data
supports.

On Patient16 the AUC separating epithelium from the stromal null moves by **0.003** across
all four (0.972 - 0.975); on Patient52 it stays at chance for every one of them
(0.452 - 0.525). Splitting finer is not free - five groups blank 1.76% of genes outright
against 0.84% for two, and push the epithelium sitting within 1.5x of its cut from 12.5% to
23.1% - and T + B is worse on every metric while putting the band on the immunoglobulin
loci. The conclusion to quote is not that the current reference is optimal but that **the
reference composition is not a lever on this call**; it is kept at two groups because that
retains a protection which costs nothing.

## What has to be re-run

The subclusters live only in the object inferCNV returns, so `summary/*_cnv.csv` from before
this change has no `subcluster` column and `call_malignant_newcnv.ipynb` will refuse it. The
previous call is preserved as `summary_v1/` and `cohort_census_v1.csv`.

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
./infercnv_all.sh --force --threads 20      # --force, or prepare skips the existing inputs
#   then call_malignant_newcnv.ipynb
python3 recelltypist_nonmalignant.py --tag newcnv
```

`--force` is not decoration: without it `prepare_infercnv_input.py` reports `[have]` for
every cohort that already has an `annotations.tsv` and the runs proceed on the old inputs,
which have no `immune_heldout` and an unstratified null.
