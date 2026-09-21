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

## The vanished-dimension control

05_4 - 05_8 keep **every** dimension and **every** direction. The reasoning is in
`utils/signature_common.py` next to `PRUNE_VANISHED` and it is a choice, not an oversight:
dropping dimensions before the correlations and the ORA decides, ahead of the analysis, which
axes are allowed to mean something. `var['vanished']` is read and reported throughout, so any
dimension that does come out significant can be checked against its flag.

`PRUNE_VANISHED=1` is the control that measures what that choice costs. It restores DRVI's own
behaviour everywhere at once and **writes nowhere near the reported run**: `OUT_TAG` moves the
tables to `tables_pruned/`, the figures to `figures_pruned/`, and 05_7's two `.tsv` to a
`_pruned` name. The per-cell files pruning cannot change — CytoTRACE2 and
`signature_scores_<coll>_<run>.csv` — are read from their existing paths and left untouched, so
the control needs no re-scoring and 05_5 never runs twice.

```bash
PRUNE_VANISHED=1 N_LATENT=64 ./signature_interpretation_all.sh --collection scie
```

The two routes prune at **different granularities**, and that is the whole finding:

| | kept on `drvi_tum_64_nomt` | dropped |
|---|---|---|
| Route A (05_6) | 56 of 64 **dimensions** | DR 57 - DR 64, the 8 flagged `vanished` |
| Route B (05_7) | 111 of 128 **directions** | those 16, plus **DR 55+** alone |
| Route C (05_8) | 112 rows | inherits Route A's list |

Route A's numbers do not move on a dimension it keeps — the Spearman correlations are computed
per dimension and cannot see each other — so every difference below is Route B's, and all of it
comes from the BH denominator shrinking from 128 to 111 directions: the same p-values, corrected
over fewer tests, so FDRs fall and a few rows cross the bar. Nothing here is a new measurement.

| collection | rows lost to pruning | of them significant on ≥1 route | verdicts changed among the 112 rows both runs share |
|---|---|---|---|
| `scie` | 16 | 13, incl. **7 convergent** | 1 |
| `emt` | 16 | 1 | 14, all `neither`/`cell_only` → `factor_only`/`convergent` at FDR 0.0538 → 0.0482 |
| `gavish_tnbc` | 16 | 15, incl. 11 `both_routes_different_family` | 1 |

Two things to take from it.

**The vanished tail is not empty.** On `scie`, seven of the 37 convergent rows live on dimensions
DRVI flagged vanished — DR 59+ at Route B FDR 7e-5, DR 62+ at 1e-3 — and on `gavish_tnbc`
fifteen of the sixteen dropped rows were significant on at least one route, DR 63+ at FDR 4e-9.
These carry ~1e-05 of the latent variance and are exactly the axes DRVI considers unused; that
they still land on interferon and EMT metaprograms is a statement about the ORA on 200-gene
lists, not evidence of a state. But it is the reason not to prune silently: pruning removes them
without ever putting them in a table where they could be dismissed on their merits.

**The `emt` flips are a threshold artefact, not a result.** All fourteen sit at FDR 0.0538 →
0.0482 across a 0.05 bar. A verdict that moves because seventeen untested directions left the
denominator is a verdict that was never resolved; read the FDR, not the label.

**The one real loss is `DR 55+`.** DRVI flagged that *direction* vanished while keeping the
dimension, so Route A keeps DR 55 and Route B stops testing its positive side. 05_8 still writes
a `DR 55+` row, reads no Route B value for it, and reports FDR 1.0 — where the unpruned run has
**7e-4**, the collection's strongest hit on that axis. This is the asymmetry `interpretability_
scores` documents, and on this object it costs a real result. It is the single strongest
argument for the default.

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
  On `drvi_tum_64_nomt`, with the 22-list default, that is **39 of the 128 dimension-directions**,
  led by `DR 3+ = MP1_CELL_CYCLE_G2_M` (ρ 0.845, FDR 4e-46), `DR 8+ = MP5_STRESS` (ρ 0.811) and
  `DR 5+ = MP2_CELL_CYCLE_G1_S` (ρ 0.775). On `drvi_epicnv_64_nomt` it is 33 of 128. The earlier
  figure quoted here, 28 of 64 on `drvi_tum_32`, was the 27-list set on the 32-dimension run and
  has not been recomputed.

  **MP1 is why the top row changed.** Before it existed, DR 3+ had no metaprogram of its own and
  read as `MP2_CELL_CYCLE_G1_S` — DRVI had separated G2/M from G1/S into two dimensions and the
  vocabulary could only name one of them. MP1 now converges on three directions (`DR 3+`,
  `DR 25-`, `DR 37-`), and `DR 3+` is the single strongest convergent row of the run.
* **an independent check on `emt`.** MP12–MP16 are four EMT metaprograms plus the glioma
  mesenchymal one, curated by other people from other tumours. An EMT axis visible on the
  collaborator's lists *and* on those is an axis that does not depend on whose EMT list was used.

Eleven of the 41 describe lineages that cannot be in this compartment (neural, skin
pigmentation, haematopoietic). In the **full** collection they are marked `primary=False` and
kept on purpose: they are the floor the other readings are measured against, and `MP36_IG`
doubles as this dataset's ambient-immunoglobulin readout.

`MP1` (Cell Cycle - G2/M) is **not** in `datasets/GAVISH.csv` — the export has 40 columns and
not 41, and the MSigDB release it came from has no MP1 either. It is now on disk anyway, taken
from the authors' own `MP_list.RDS` (github.com/tiroshlab/3ca, `$Cancer[[1]]`, 50 genes in the
paper's order) and written by `utils/gavish_extraction.py` next to the forty from the CSV. The
two sources are the same list: that object's MP2 is identical gene-for-gene to
`MP2_CELL_CYCLE_G1_S.txt` apart from two symbols MSigDB updated (`HIST1H4C`→`H4C3`,
`KIAA0101`→`PCLAF`). MP1 keeps the paper's spelling, which is why all 50 of its genes map here:
this dataset is on an older reference and carries `HIST1H4C` but not `H4C3`.

**Only 22 of the 41 are scored by default.** Seventeen name a lineage or a tissue a
triple-negative breast carcinoma cannot express, and each of them is a column of both heatmaps
and one more test inside the FDR correction of 05_7. `--collection gavish` therefore scores the
21 states reported in basal-like / TNBC malignant cells — the whole cycle (MP1 G2/M, MP2 G1/S,
MP3 HMG-rich), chromatin (MP4), stress and hypoxia (MP5, MP6), proteostasis (MP8–MP10), the four
EMT programmes (MP12–MP15), interferon / MHC-II (MP17, MP18), epithelial senescence (MP19), MYC
(MP20), respiration (MP21), the secreted pair (MP22, MP23) and metal-response (MP39) — plus
`MP7` (in-vitro stress), the one list kept only to bound them: it is what separates a stress
dimension from a dissociation one.

**Two lists were measured out rather than argued out: `MP38` and `MP11`.** The null is 200
random 50-gene lists, matched bin-for-bin to the expression profile of the real metaprograms,
scored with the Route A settings and correlated against the same 64 dimension-directions of
`drvi_tum_64_nomt`; its max |ρ| has p50 0.225, p95 0.325, p99 0.372 and never exceeded 0.432 in
200 draws. `MP38` (Glutathione) reaches **0.165**, the 25th percentile of that null — a random
list beats it three times in four — and its genes (`ACSM2A/B`, `AGXT2`, `CUBN`, `AMN`, `CLTRN`,
`FOLR1`, `AQP1`) are proximal-tubule renal, so it belongs with MP30 / MP31 / MP40. `MP11`
(Translation initiation) reaches **0.216**, the 42nd percentile, while firing 18 significant
pairs on Route B: `EIF2`/`EIF3` genes land in the top decoder genes of several dimensions
without the per-cell score tracking any of them, which is what a library-complexity artefact
looks like from both sides. The null is a Methods number, not a threshold applied anywhere in
the pipeline, and these two are the only decision taken with it.

What the null does **not** license is removing a list for having few genes inside the 2,000-HVG
panel: that count predicts nothing about whether a programme can be found. `MP21` has 2 of them
and fires 23 times on Route B; `MP2` has 28 and fires 6. `MIN_SIGNATURE_GENES` is a floor on the
MAPPED count, never on the HVG one, and the HVG count is reported as a warning for that reason.


**No lineage control is scored by default any more,** and that is a decision taken by reading
the full 41 rather than the earlier default, which kept five of them. What it costs is worth
knowing: there is no longer a row saying what a correlation of nothing looks like on this data,
the two ambient-RNA readouts (`MP36` IG, `MP33` RBCs) are gone with it, and `MP41` (unassigned)
is no longer there as the sink for a dimension that matches nothing. The first is what
`--all-metaprograms` is for — the same tables under the other slug, once per latent space. The
second is the least costly here: 05's own diagnosis put the ambient contribution about two
orders of magnitude below the dimensions it could have explained, and SoupX in phase 06 checks
the same thing without going through a metaprogram. `primary=False` and the
`lineage_control_claim_bounds_the_rest` flag therefore never fire on the default variant; they
still fire on `gavish`.

The set is one editable tuple, `_GAVISH_TNBC_MPS` in `utils/sig_collections.py`; `MP30` and
`MP40` are the first candidates to add back, their genes being the generic secretory-epithelial
ones that are luminal in breast, and `MP16` / `MP24` after them.

The two widths are two collections on disk — slug `gavish_tnbc` and slug `gavish` — so their
tables and figures never share a folder or a filename.

**What has been re-run on the 22-list set.** `05_4` → `05_8` on `drvi_tum_64_nomt` and
`drvi_epicnv_64_nomt`. After that run, `tables/gavish_tnbc/` and `figures/*/gavish_tnbc/` hold
those two runs and nothing else, so everything currently under `tables/` and `figures/` for this
collection is the 22-list set. What is superseded is in `tables_v1/` and `figures_v1/`: the
27-list `gavish_tnbc` on `drvi_tum_32`, `drvi_tum_64` and `drvi_epicnv_64_nomt`, and the whole
`gavish` slug, which was the 40-list registry from before MP1 existed. Nothing there has been
recomputed; `--all-metaprograms` on the 41 has not been run since MP1 was added.

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

# 05_4 - 05_8 in one command, resuming; see ./signature_interpretation_all.sh --help
cd .. && ./signature_interpretation_all.sh --collection scie      # or emt | gavish
PYTHON=~/miniconda3/envs/cytotrace2-py/bin/python \
    ./signature_interpretation_all.sh cytotrace                   # the one step in its own env

# or by hand, which is what the driver does:
cd 05_4_signatures       && python3 build_signatures_tum.py   # + --collection emt|gavish [--all-metaprograms]
cd ../05_5_cytotrace2                                         # in the cytotrace2-py env
python3 cytotrace2_tum.py
cd ../05_6_cell_first    && python3 cell_first_tum.py         # Route A   + --collection emt|gavish
cd ../05_7_factor_first  && python3 factor_first_tum.py       # Route B   + --collection emt|gavish
cd ../05_8_convergence   && python3 convergence_tum.py        # Route C   + --collection emt|gavish
```

`signature_interpretation_all.sh` is the phase-level driver, the counterpart of 04's. It walks
the six steps in the one order they work in — 05_7 reads the row order off the table 05_6
writes, and 05_8 reads five tables from the two — resumes by output existence, and derives the
run id and the output tag from `utils/resolve_run.py` rather than re-spelling either rule in
bash. Everything else it does is pass the environment through.

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
