# 05_3_drvi_run

DRVI on the malignant subset. The counterpart of
[04_2_drvi_run](../../04_drvi_epithelial/04_2_drvi_run/), on the object 05_2 wrote instead of
04_1's, with the same architecture and the same seed — so the only thing that differs between
the two latent spaces is the cells they were trained on.

Input: `$DATA_DIR/05_tum/shiao_tum_hvg_2k.h5ad` — **36,192 cells × 2,000 HVGs**, 19 cohorts,
the HVGs re-selected inside the tumour by `05_2/reduce_data_tum.py`. `batch_key = 'cohort'`,
`.layers['counts']` are the raw counts DRVI trains on.

Nothing is inherited from 04 — not the model, not the HVGs, not the latent size — and the run
ids do not collide (`drvi_tum_*` here against `drvi_epi_*` there, `drvi_nonimm_*` in 03_2 and
`drvi_unscaled_*` in 02_2), so every phase keeps its own runs.

## The files

| File | What it is | Where |
|---|---|---|
| `drvi_tum.ipynb` | The step. Trains, builds the embedding, draws every figure, screens the factors with SMI, then names them by their genes (Enrichr / g:Profiler / decoupler). This is where the latent size is **judged**. | local (notebook) |
| `run_drvi_tum.py` | The same computation without a kernel: same input, same parameters, same outputs, same figures. Not the SMI and enrichment sections, which are read rather than run. For a terminal or the cluster. | local / cluster |
| `submit_drvi_tum.slurm` | One SLURM job wrapping the runner. CPU, `long` partition, 16 cpus / 32G, no `--time`. | cluster |

The two halves are deliberately redundant: train wherever it is convenient, then open the
notebook with `OVERWRITE = False` and it reads the model and the embedding off disk instead of
recomputing them.

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
python3 run_drvi_tum.py                  # n_latent 32, the run of this step
python3 run_drvi_tum.py --n-latent 64    # another size, beside it, nothing overwritten
python3 run_drvi_tum.py --overwrite      # retrain and rewrite everything
CELL_SET=epi python3 run_drvi_tum.py     # the same run on the control set
```

`CELL_SET` works here as in 05_2 and through the same `05_2_subsetting/cell_set.py`, imported
rather than duplicated: `tum` (default) is the malignant subset, `epi` the control set of all
epithelium under the post-CNV labels. They produce different run ids — `drvi_tum_32` against
`drvi_epicnv_32` — and therefore different models, embeddings, figures and tables.

## Why 32 latent dimensions

`N_LATENT` is not tuned against a score; it is read off the **vanished count**, the number of
dimensions the KL pressure emptied out because the model did not need them. Both the notebook
and the runner print it as `n vanished / n_latent`, and it is the one number to look at before
deciding anything else:

- **several vanished** → the size was already generous, keep it;
- **none vanished** → the space is too tight to say what DRVI did not need, and the run is
  worth repeating at twice the size.

32 is where this step starts, and the argument is 04's own read backwards. 04_2 tried 32 first
— on 74,441 cells and 10 CellTypist labels — got **0 / 32** vanished, and settled at 64. This
subset is half of that (36,192 cells) and, unlike 04's, is a single compartment whose
annotation is one constant value, so there is less for the model to encode, not more. If 32
still leaves nothing vanished the answer is the same one 04 arrived at: re-run at 64 and
compare. The run id carries the size, so that costs a training and no bookkeeping.

Everything else is 04_2 verbatim: encoder `[256, 128]` / decoder `[128, 256]`,
`dispersion='gene-batch'` (the parameter that matters on this TNBC dataset, and the reason
05_2 kept the 200-cell cohort floor — a dispersion column per gene *and* per batch cannot be
fitted on a cohort of 14 malignant cells), 400 full epochs, no early stopping, `SEED = 123`.

> **Why no early stopping.** scvi-tools raises the KL weight linearly over
> `n_epochs_kl_warmup`, 400 by default, so a run that stops earlier never trains at full KL
> weight — and it is that pressure that makes the unused dimensions vanish, i.e. produces the
> statistic the latent size is chosen on. `--early-stopping` exists and only makes sense
> together with a shorter warmup.

## The one thing that does not carry over from 04_2: the label

`cell_type` in this object is **`malignant` for all 36,192 cells**. It was written by 05_2 from
the CNV call, not by CellTypist, and it is what makes phase 05 different from phase 04. What it
costs is every figure 04_2 draws *by cell type*: the per-label UMAP panels, the
latent-dimension heatmaps, the SMI screening. Drawn against one category, all three say
nothing.

Two groupings replace it, **in this order**, and they are not interchangeable:

| # | grouping | what it is | how a hit reads |
|---|---|---|---|
| 1 | `optscib_tum_leiden` | 05_2's clustering, computed on these very cells | a **tumour state** on the subset's own terms |
| 2 | `cell_type_01_4` | the CellTypist label these cells carried **before** the CNV call, counted in `shiao_tum_hvg_2k.h5ad` itself: 18,617 `Lumsec-prol`, 11,771 `Lumsec-basal`, 3,893 `LummHR-SCGB`, 1,167 `Lumsec-KIT`, 594 `LummHR-major`, 148 across four more | a **landmark**. The tumour state has a normal-breast counterpart the model can put a name to — never evidence that the cell *is* that normal type |

The order is the argument. The first key is the one that leads: it gets the extra grouped
heatmap, it comes first in the UMAP panels, and it is the first SMI target. If the leiden
column cannot be attached, `cell_type_01_4` takes that slot by default rather than by choice.

### Isn't `cell_type_01_4` obsolete?

**As an identity, yes — establishing that is what phase 05 is for.** It is kept as a
*grouping* because on this subset it is the only non-constant CellTypist column left:

```
cell_type                 malignant  36,192                     (constant)
celltypist_predicted_cnv  malignant  36,192                     (constant: 05_1 re-ran CellTypist
                                                                 on the NON-malignant cells only)
celltypist_predicted      Lumsec-prol 9,877 … Fibro-matrix 1,871, pericytes 658, …
                                                                (raw, pre-voting: calls fibroblasts
                                                                 and pericytes on aneuploid cells)
cell_type_01_4            Lumsec-prol 18,617, Lumsec-basal 11,771, …
```

So the alternatives are a constant, a noisier version of the same wrong model, or nothing.
This is the same substitution `05_2/clustering_tum.py` makes for its NMI target, with the same
caveat stated in the same place it is used: a *borrowed grouping used to read the latent
space*, not a biological claim about these cells. Whatever 05_4 onwards writes about "states"
has to keep that distinction and that order — the borrowed label is there to *name* a dimension
the leiden partition and the gene programs have already established, never to establish one.

**Neither grouping is clean, and they fail differently.** The CellTypist label comes from a
model with no malignant class. The leiden partition comes from the **unintegrated** PCA of
05_2, so it carries part of the cohort structure DRVI is meant to correct — which is why every
heatmap is drawn against `cohort` as well — and it is not fully independent of the label above
either: `clustering_tum.py` picked its **resolution** by maximising NMI against
`cell_type_01_4`, which is recorded in `uns['optscib_tum_leiden_label_key']`. The clusters come
from the expression graph; the borrowed label chose how many of them there are, not which cell
goes in which. A dimension that answers to both groupings and not to `cohort` is the safe
reading; one that answers to leiden alone is a state the normal-breast vocabulary has no word
for, which is what this phase is looking for.

### None of this is an input to DRVI

Worth stating once, because it is what makes an imperfect grouping affordable at all:

```python
DRVI.setup_anndata(adata, layer="counts", batch_key=BATCH_KEY)   # and nothing else
```

`setup_anndata` takes a `labels_key` and it is left `None`. The counts layer and `cohort` are
the entire input, so the latent space is built without knowing any of these labels and **none
of them can bias it**. They are read back afterwards, as captions on dimensions that already
exist — a wrong caption can mislead a reader, it cannot corrupt the model. That is the argument
for keeping a borrowed grouping rather than dropping it: dropping it removes the caption, not a
source of error, and leaves 32 anonymous axes.

The caption that owes nothing to any annotation is the third one, and this step writes it too:
the **OOD / IND interpretability scores**, which name a dimension by its genes. That is the
route 05_6 (cell-first) and 05_7 (factor-first) take. The two groupings here are the
cross-check on it, never the evidence.

The leiden column is not in `shiao_tum_hvg_2k.h5ad` — `reduce_data_tum.py` wrote that file
before `clustering_tum.py` ran — so both the notebook and the runner fetch it from
`shiao_tum.h5ad`, reading only `obs` (a few MB out of ~330) and realigning by cell name. If
that object is not around, the leiden figures are skipped and nothing else changes.

Constant columns are dropped with a line in the log saying which and why, rather than drawn:
`compartment`, `fraction`, `cnv_status` and, under `tum`, `cell_type`. Under `CELL_SET=epi`
`cell_type` is *not* constant and comes back as the interesting key — one flag, one code path.

## Two panels 04_2 has no equivalent for

`cnv_score` and `cnv_corr`, the two quantities 05_1 called the subset with, are plotted on the
latent UMAP. They are the check that this space is about tumour **states** and not about how
aneuploid a cell is: a gradient running along a dominant dimension would mean the model spent
its capacity on CNV burden, and that would have to be read into everything downstream. Their
unintegrated counterparts are `figures/05_2_reduce_data/umap_tum_cnv_score.png` and
`umap_tum_cnv_corr.png` — same cells, same keys, only the space changes.

## What comes out

`$DATA_DIR/05_tum/` — never the repo: the embedding alone is a few hundred MB, past GitHub's
per-file limit, and `datasets/*` is gitignored.

```
model_<run_id>.pt              the trained model, one flat file per run
embed_<run_id>.h5ad            latent space + per-dimension stats + OOD/IND gene scores
shiao_tum_<run_id>.h5ad        the 05_2 object (all genes) + obsm['X_drvi']
```

with `<run_id>` = `drvi_tum_32`. Figures go to `../figures/05_3_<run_id>/` and tables to
`../tables/05_3_<run_id>/`.

- **`embed_<run_id>.h5ad` is the file 05_4 reads.** It needs neither the model nor a GPU.
- `shiao_tum_<run_id>.h5ad` is for the steps that need genes and latent coordinates in the same
  object (05_5 cytotrace2, 05_6 cell-first, 05_9 cycle-confound). It is the only output that
  needs `shiao_tum.h5ad` present; without it the run says so up front, trains, and writes the
  other two.
- Nothing is written for the scib benchmark: 02_3 and 02_4 do not run on this compartment.

Resuming is the default, as everywhere in this phase: a step whose output is already on disk is
reported as `[have]` and reused, so a crash after training does not cost the training.
`--overwrite` (or `OVERWRITE = True`) recomputes everything.

## On the cluster

```bash
cd 05_drvi_tumoral_epi/05_3_drvi_run && mkdir -p logs
export DATA_DIR=/users/genomics/albertoc/Tesi/hopes_and_dreams/datasets
sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_drvi_tum.slurm                # n_latent 32
sbatch --export=ALL,DATA_DIR=$DATA_DIR submit_drvi_tum.slurm --n-latent 64  # another size
```

Everything after the script path reaches `run_drvi_tum.py` unchanged. Three things have to be
up there first: `shiao_tum_hvg_2k.h5ad`, `shiao_tum.h5ad` (optional, see above) and the sibling
`05_2_subsetting/cell_set.py`, which the runner imports — the job checks for it in its first
second rather than after loading the model. `long` rather than `normal` because the cluster has
no GPU; 32G because this is the lightest DRVI job of the thesis by a wide margin.

## Naming the factors: two sections, two vocabularies

After the run itself the notebook asks what the dimensions are *about*, twice:

**SMI**, against both groupings, as described above — "which factor behaves like a leiden
cluster / like a normal-breast label".

**Enrichment** (the DRVI tutorial's *Identifying biological processes of DRVI factors*): the
same per-gene interpretability scores handed to three tools — Enrichr and g:Profiler for
over-representation, decoupler for TF activity — with the 2,000-HVG panel as the declared
background. It matters more here than in 04_2: this phase has no non-constant biological
annotation, so naming a dimension **by its genes** is the only reading that owes nothing to a
borrowed vocabulary, and the SMI matches are the cross-check on it rather than the evidence.

It writes `../tables/05_3_<run_id>/factor_processes_<run_id>_top<N>.csv` — one row per factor
direction, its SMI match against *each* target, the top hit of each tool, and its top genes.
That file is the notebook-level counterpart of what `05_7_factor_first` does against the
collaborator's `.gmt`: same scores, same background, a different question.

Needs network (Enrichr's library, g:Profiler's service, OmniPath) and `gseapy`,
`gprofiler-official`, `decoupler` — all three already in `benchmark-py-r`.

> **What to expect on this subset.** The cell cycle will enrich somewhere: `Lumsec-prol` is 51%
> of these cells. A mitotic dimension is a real fact about the axis and also the confound 05_6
> recomputes its target region inside G1 to measure. Read the two together.

## What this step does not do

04_2's notebook also carries the rare-state scan (*Finding rare (un-annotated) cell types with
DRVI*), which asks which dimensions describe a population the annotation does not have. On this
subset the annotation is one constant value, so "a sub-population the labels do not resolve"
would have to be defined against the leiden partition instead — a real question, and a different
one from the one 04_2 asks. It belongs with 05_6 (cell-first) rather than here.
