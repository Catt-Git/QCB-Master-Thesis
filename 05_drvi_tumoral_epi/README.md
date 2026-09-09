# 05_drvi_tumoral_epi

Sixth phase of the thesis, and the third inside **Part 2 (biological interpretation)**.

Phase 04 read DRVI's latent dimensions on the epithelial compartment and said, in its own
README, exactly what it could not claim: `cell_type` comes from a CellTypist model trained on
the **normal** adult breast, that model has no malignant class, and so every "epithelial state"
it reports is a state of a mixture of normal and tumour epithelium. This phase removes that
limitation instead of restating it. It calls copy-number variation per cell, separates the
aneuploid epithelium, re-annotates what is left with the same CellTypist model on a population
the model was actually trained for, and then applies 04's procedure to the tumour.

**Phase 04 is not modified and is not re-run.** 04 and 05 are two branches off `shiao.h5ad`
that share no output file: 04 keeps its results and its honest caveat, 05 carries its own. The
step numbering of the two runs in parallel, offset by one, because 05 has an extra step at the
front:

| 05 | 04 | |
|---|---|---|
| `05_1_infercnv/` | - | inferCNV → `malignant` / `non_malignant`, CellTypist re-run |
| `05_2_subsetting/` | `04_1_subsetting/` | from `shiao.h5ad` to the malignant object |
| `05_3_drvi_run/` | `04_2_drvi_run/` | DRVI on that object |
| `05_4_signatures/` | `04_3_signatures/` | the lists of one collection → `.gmt`, coverage, Jaccard |
| `05_5_cytotrace2/` | `04_4_cytotrace2/` | per-patient potency |
| `05_6_cell_first/` | `04_5_cell_first/` | Route A |
| `05_7_factor_first/` | `04_6_factor_first/` | Route B |
| `05_8_convergence/` | `04_7_convergence/` | Route C |
| `05_9_cycle_confound/` | `04_8_cycle_confound/` | how much of "stemness" is the cycle |

Input: `shiao.h5ad` (619,693 × 30,869) from phase 01, read-only. Integration
`batch_key = 'cohort'`.

## The cell set: `tum` by default, `epi` as a control

`CELL_SET` decides what the phase is about. It is one flag, not two copies of the chain, and
the two sets write different file prefixes so both can be run in the same directory:

| `CELL_SET` | cells | prefix | role |
|---|---|---|---|
| `tum` (default) | `cnv_status == 'malignant'` | `shiao_tum_*` | the phase |
| `epi` | all epithelium under the post-CNV labels | `shiao_epicnv_*` | a control |

A second flag, `HVG_SET`, is orthogonal to it and picks **which 2,000 genes** the DRVI input is
made of: unset is the panel `05_2/reduce_data_tum.py` selected, `nomt` the one
`05_2/hvg_no_mt.py` rebuilt with the 11 mitochondrial (`MT-`) genes replaced by the next 11 of
the same batch-aware ranking. It tags the files and the run id the same way `CELL_SET` does
(`drvi_epicnv_64_nomt` beside `drvi_epicnv_64`), so nothing is overwritten and the two are read
against each other. See `05_2_subsetting/README.md`.

**Why `tum` is the primary line.** 05_6 and 05_7 read the DRVI latent dimensions and ask which
gene programmes load on them. If DRVI is trained on all epithelium, its dominant dimensions
encode malignant-versus-normal — a contrast that is nearly **constant** inside the malignant
subset — so those two steps would be interpreting what is left over after the biggest axis has
been spent on a question they are not asking. This is the argument 04 makes for existing at all
against 03 (*"its HVGs, its PCA and its latent dimensions describe fibroblasts and endothelium
at least as much as epithelium"*), one compartment further down.

**What `epi` is for.** It answers "how much did the wrong labels cost phase 04?" — a methods
question, worth a paragraph, and not the biology. Its results are a comparison against 04, not
a second main line.

## What the malignant subset costs

19 cohorts, **36,192 cells**, against 29 cohorts and 74,441 cells in 04. Eleven of 04's cohorts
do not reach 200 malignant cells — a cohort can have 1,649 epithelial cells and 102 aneuploid
ones.

The `MIN_CELLS_PER_COHORT = 200` rule is inherited from 04_1 and here it is nearly free,
because the per-cohort malignant counts are bimodal — a cohort has thousands of malignant cells
or a handful, with almost nothing between:

| threshold | cohorts | cells |
|---|---|---|
| ≥ 50 | 27 | 36,901 |
| ≥ 100 | 21 | 36,433 |
| ≥ 200 | **19** | **36,192** |
| ≥ 300 | 18 | 35,991 |

Going from 100 to 200 costs 241 cells out of 36,433 (0.7%) and buys back per-cohort batch
parameters DRVI can actually estimate. 04_1's other filter — dropping cohorts missing a
treatment timepoint — is **off** here: it existed for a treatment phase that is no longer
planned, and on this subset several cohorts clear it with two cells in a timepoint while
Patient63 would lose 1,673 malignant cells for having none.

> **A limit to state, not to hide.** Keeping the cohorts with the most tumour is a selection,
> and "most tumour" could correlate with something biological. Against treatment response the
> eleven dropped cohorts split 4 NR / 4 R1 / 3 R2 against 6 / 4 / 8 kept — no visible skew, but
> with 8-11 cohorts per group a moderate one would not be visible either.

## Three things that do not carry over from 04

Duplicating 04's scripts unchanged would break these silently. Each is a decision, and each is
written down where it is taken.

**1. The leiden resolution has no NMI target.** 04_1 picks the resolution by maximising NMI
against `cell_type`. On the malignant subset the post-CNV `cell_type` is the constant
`malignant`, and NMI against a constant is zero everywhere. `05_2/clustering_tum.py` falls back
to `cell_type_01_4`, the pre-CNV CellTypist label, which is **not** constant inside the tumour —
18,943 `Lumsec-prol`, 12,028 `Lumsec-basal`, 3,975 `LummHR-SCGB`, 1,174 `Lumsec-KIT`. It is a
borrowed label used to pick one number, never propagated as a biological claim.

**2. Route A must not standardise within `cell_type`.** 04_5 uses
`GROUPBY = ["cohort", "cell_type"]`, and there `cell_type` was a **lineage** — luminal versus
basal are different cells, and standardising within them asks "inside a cell type, does this
dimension track stemness?". Inside one malignant compartment those groups are **states** of the
same tumour, and state is the quantity this phase measures. 05_6 therefore uses
`GROUPBY = ["cohort"]` only, and reports the pre-CNV label as a covariate rather than
regressing it out. The concrete danger of not doing this: `Lumsec-prol` means *proliferating*,
and the cell cycle is a named risk of the `scie` collection — standardising within that label
would remove the proliferation axis **by construction** instead of measuring it, which is what
05_9 exists to do.

**3. The caveat text is new.** 04's `CAVEAT` says no CNV inference has been run, which stays
true of 04's tables and figures. Phase 05 needs its own wording, describing its own limits
(19 cohorts, a threshold-based call, `not_tested` cells excluded), not a patch to 04's.

## What EMT can and cannot be asked here

The `emt` collection's target region is *epithelial-high × mesenchymal-high*, i.e.
co-expression, with the HYBRID lists validating the call rather than making it. Running on the
malignant subset **narrows a named risk of that collection**: 04 lists "fibroblast ambient RNA /
doublets" as a way a cell could score mesenchymal without being in a hybrid state, and a
high-mesenchymal cell here has to be aneuploid to be in the object at all. Ambient RNA survives
— it is contamination, not identity — but an actual fibroblast does not.

What it does not buy is discreteness. Partial EMT in carcinoma is a continuum; expect to place
each tumour cell on an E–M axis and identify the co-expressing region, not to find three leiden
clusters that name themselves. A sparse mesenchymal tail is the usual result and is a result.

## Which run 05_4 - 05_8 read

Those five steps train nothing: they read what 05_3 wrote, and what they need from it is its
**name**. That name has three segments, each an environment variable resolved by
`05_2_subsetting/cell_set.py` and nowhere else, so 05_3 and the interpretation chain cannot
spell the same run differently:

| variable | default | what it changes |
|---|---|---|
| `CELL_SET` | `tum` | `tum` the malignant subset, `epi` the all-epithelium control → `drvi_tum_…` / `drvi_epicnv_…` |
| `N_LATENT` | `32` | the DRVI latent size 05_3 was run at |
| `HVG_SET` | unset | `nomt` selects the gene panel without the MT- genes → the `_nomt` suffix |

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets

python3 build_signatures_tum.py --collection gavish                    # drvi_tum_32
N_LATENT=64 python3 build_signatures_tum.py --collection gavish        # drvi_tum_64
CELL_SET=epi HVG_SET=nomt N_LATENT=64 \
    python3 build_signatures_tum.py --collection gavish                # drvi_epicnv_64_nomt
```

05_3 keeps its `--n-latent` flag and the flag wins; its default is `$N_LATENT`, so one export
moves the whole phase. These five steps deliberately have **no flag of their own** for any of
the three: a flag would be a fourth place the run id is spelled, and the failure that invites is
silent — a step reading last week's embedding and writing a table that says nothing about it.

Nothing is shared between two runs. Every table and figure carries the run id in its filename
**and in a folder of its own**, so several runs of one collection sit side by side without
either the reader or the code having to filter a directory by suffix:

```
tables/gavish/
    drvi_tum_32/          coverage_gavish_drvi_tum_32.csv, convergence_..., 12 files
    drvi_tum_64/          the same twelve, for that run
    drvi_epicnv_64_nomt/  and for that one
    signatures_gavish_tum.gmt      <- above the run level, see below
    signatures_gavish_epicnv.gmt

figures/05_8_convergence/gavish/drvi_tum_64/routes_side_by_side_gavish_drvi_tum_64.png
```

The run id stays in the filename as well, and that is not redundancy: a file dragged out of
its folder has to remain identifiable, and the folder is there to make the directory readable
rather than to make the name unique.

The one exception is the `.gmt`, which depends on the **object** rather than on the run and is
therefore named `signatures_<collection>_<compartment>.gmt` and sits one level up, beside the
run folders rather than inside one — `shiao_tum.h5ad` has 24,779 genes
and `shiao_epicnv.h5ad` 26,379, so a single name would have let a `CELL_SET=epi` run overwrite
the malignant set's mapping in place and 05_7 would then have done its ORA against the wrong
gene universe without complaining. `N_LATENT` and `HVG_SET` are **not** in that name, and must
not be: neither changes which genes the object has.

## The third collection: naming the dimensions from outside

`scie` and `emt` each ask whether a **named state exists**, and each defines a region of the
cell-first plane to call cells in. A third collection, `gavish`, asks the question the other way
round — *given a latent dimension, what is it?* — and therefore defines **no** target region.
`Collection.has_target` is False for it, and each step drops exactly what depends on a region:
05_6 writes the per-cell scores, the confounder table and the dimensions × metaprograms
correlations but skips A5, the consensus quadrant and the two figures that draw it; 05_8 runs
unchanged except for the target-axis test at the end. 05_4 and 05_7 are untouched.

The lists are the pan-cancer metaprograms of Gavish et al. 2023 (Nature 618:598-606), in
`$DATA_DIR/signatures/GAVISH_metaprograms/`, derived by NMF over ~1,000 tumours of 24 cancer
types with no knowledge of this project. Two things they buy:

* **a name from outside.** A dimension on which Route A and Route B independently land on the
  same metaprogram is a dimension named by prior knowledge that was not chosen for this project.
  On `drvi_tum_32` that is 28 of the 64 dimension-directions, led by `DR 9- = MP5_STRESS`
  (ρ 0.813, FDR 7e-40), `DR 5- = MP2_CELL_CYCLE_G1_S` and `DR 17+ = MP6_HYPOXIA`.
* **an independent check on `emt`.** MP12–MP16 are four EMT metaprograms plus the glioma
  mesenchymal one, curated by other people from other tumours. An EMT axis visible on the
  collaborator's lists *and* on those is an axis that does not depend on whose EMT list was used.

Eleven of the 40 describe lineages that cannot be in this compartment (neural, skin
pigmentation, haematopoietic). They are marked `primary=False` and kept on purpose: they are the
floor the other readings are measured against, and `MP36_IG` doubles as this dataset's
ambient-immunoglobulin readout. `MP1` (Cell Cycle G2/M) is absent from `datasets/GAVISH.csv` and
so from the text files — there are 40 lists, not 41.

**Only 27 of the 40 are scored by default.** Eighteen name a lineage or a tissue a
triple-negative breast carcinoma cannot express, and each of them is a column of both heatmaps
and one more test inside the FDR correction of 05_7. `--collection gavish` therefore scores the
22 states reported in basal-like / TNBC malignant cells — cycle (MP2, MP3), chromatin (MP4),
stress and hypoxia (MP5, MP6), proteostasis (MP8–MP11), the four EMT programmes (MP12–MP15),
interferon / MHC-II (MP17, MP18), epithelial senescence (MP19), MYC (MP20), respiration (MP21),
the secreted pair (MP22, MP23) and the detoxification pair (MP38, MP39) — plus five kept only to
bound them: `MP7` (in-vitro stress, the dissociation control on the stress axis), `MP36` (IG) and
`MP33` (RBCs), the two ambient-RNA readouts that phase 06 removes with SoupX, `MP25`
(astrocytes) as the impossible lineage that shows what zero looks like, and `MP41` (unassigned)
as the residual sink. `--all-metaprograms` restores all 40. The set is one editable tuple,
`_GAVISH_TNBC_MPS` in `utils/sig_collections.py`; `MP30` and `MP40` are the first candidates to
add back, their genes being the generic secretory-epithelial ones that are luminal in breast.

The two widths are two collections on disk — slug `gavish_tnbc` and slug `gavish` — so their
tables and figures never share a folder or a filename, and the numbers quoted above, measured on
all 40, stay exactly where they were written.

## Status

`05_1` through `05_8` are written. `05_9_cycle_confound` is not, and neither is an
embedding-control step — 04's 04_9 has no counterpart here yet (`utils/signature_common.py`
registers DRVI only, and says what adding one would take). See each step's README.

```bash
conda env create -f ../environments/infercnv-r.yml     # once, for 05_1
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets

cd 05_1_infercnv && ./infercnv_all.sh --threads 12     # then call_malignant.ipynb
python3 recelltypist_nonmalignant.py

cd ../05_2_subsetting                                  # then subset_and_qc.ipynb
./subsetting_all.sh

cd ../05_3_drvi_run                                    # drvi_tum.ipynb, or headless:
python3 run_drvi_tum.py                                # n_latent 32, see that README

cd ../05_4_signatures && python3 build_signatures_tum.py   # + --collection emt|gavish [--all-metaprograms]
cd ../05_5_cytotrace2                                      # in the cytotrace2-py env
python3 cytotrace2_tum.py
cd ../05_6_cell_first    && python3 cell_first_tum.py      # Route A   + --collection emt|gavish
cd ../05_7_factor_first  && python3 factor_first_tum.py    # Route B   + --collection emt|gavish
cd ../05_8_convergence   && python3 convergence_tum.py     # Route C   + --collection emt|gavish
```

Steps 05_4 - 05_8 share `utils/signature_common.py` (paths, the caveat, the figure and table
writers) and `utils/sig_collections.py` (the three collections and everything that differs
between them), both duplicated from 04's `utils/` as every phase in this repo duplicates rather
than imports. What is *not* a copy is listed at the top of each: for `signature_common.py` it is
the object and run id, the caveat, the grouping keys that replace the constant `cell_type`, and
an embedding registry holding DRVI alone.

**A fourth thing that does not carry over from 04** joins the three listed above, and it is
05_3's: `cell_type` is the constant `malignant` in the DRVI object too, so every figure 04_2
draws by cell type is drawn against `optscib_tum_leiden` first (05_2's clustering, computed on
these cells) and `cell_type_01_4` second, as a landmark. The pre-CNV label is obsolete as an
identity — that is the point of the phase — and it survives as a grouping only because on the
malignant subset it is the last non-constant CellTypist column: the post-CNV re-annotation ran
on the non-malignant cells and leaves `celltypist_predicted_cnv` at `malignant` here. Same
substitution as the leiden NMI target, same caveat, written down in
`05_3_drvi_run/README.md`.
