# 05_9_embedding_control — is a signature *present* in a coordinate system?

Counterpart of [04_9_embedding_control](../../04_drvi_epithelial/04_9_embedding_control/), and
deliberately **not** a port of it. 04_9 compares two spaces through the *target region* and
through `max |rho|` per dimension; this step asks the prior question, in a form that does not
require a target region and never reads a single axis — so it runs on `gavish`, the collection
that defines no region and is the vocabulary phase 05 actually reports.

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
cd 05_drvi_tumoral_epi/05_9_embedding_control

N_LATENT=64 conda run -n benchmark-py-r ./embedding_control_all.sh       # all three collections
N_LATENT=64 python3 signature_presence_tum.py --collection gavish        # one of them
N_LATENT=64 python3 signature_presence_tum.py --spaces harmony drvi      # without the PCA arm
N_LATENT=64 python3 signature_presence_tum.py --n-random 0               # permutation null only
```

`N_LATENT=64` is not optional in practice: 05_3b only ever ran at 64, and at any other value
the Harmony embedding this step needs is not on disk. `CELL_SET` and `HVG_SET` behave exactly
as in 05_4 – 05_8.

## The question, and why it is not 05_6 pointed at another space

[05_3b's README](../05_3b_harmony_run/README.md) calls its Harmony run a deliberate dead end,
for one reason: 05_7 and 05_8 read an **additive decoder**, Harmony has none, and its
dimensions are corrected principal components whose loadings belong to the PCA *before* the
correction. That argument is about Route B and it is correct.

It says nothing about the cells. **A signature score is a number per cell**, computed by 05_6
from the gene expression matrix, and it has never seen a latent space — which is exactly why
`signature_common.scores_csv()` names it after the *object* and not after the run. So the
question can be asked without asking Harmony for anything it does not offer:

> take the per-cell score. Ask whether cells that are **neighbours in the space** have similar
> scores.

That is the whole instrument. It never reads a dimension's gene list, never names an axis and
never assumes the programme lives on one coordinate. **A method whose axes are not individually
interpretable can still carry a programme perfectly, spread over all of them**, and this step is
built so such a method scores well rather than being punished for its architecture.

## What is measured

Four numbers per (space, readout), all on the **within-cohort z score** — absolute scores are
not comparable across patients, which is 05_6's argument and is inherited here unchanged:

| column | the question it is |
|---|---|
| `morans_i` | spatial autocorrelation of the score on the *k*-NN graph of the space: do neighbours agree? One number, no axis, no direction |
| `knn_r2_cell` | out-of-sample R² of predicting a cell's score from the mean of its *k* nearest **training** neighbours, 5 folds. Does the space know the score at all? |
| `knn_r2_cohort` | the same with folds **grouped by cohort**, so the neighbours are never from the held-out patients. Does it transfer to a patient the fold never saw? The harder one, and the honest one |
| `top_decile_fold` | of the 15 neighbours of a top-decile cell, how many are top-decile themselves, over the 0.10 baseline. Is the high end a **region**? |

plus three read against each other:

| column | |
|---|---|
| `ridge_r2_cohort` | the **linear** share: ridge on the coordinates, same cohort folds |
| `nonlinear_gap` | `knn − ridge`. **Positive**: structure no linear function of the coordinates reproduces. **Negative**: the programme is a linear gradient across the whole space and averaging 15 neighbours is simply the worse way to read it — on a corrected PCA that is the *expected* sign and a defect of nothing |
| `space_r2_cohort` | `max(knn, ridge)`. **The column to quote for "is it present".** The question is the space, not which of two readers was used on it, and a space that holds a programme as a clean gradient must not score low for being easy |

and exactly one that is DRVI's own claim and is labelled as such:

| column | |
|---|---|
| `max_abs_rho`, `best_dim` | the strongest single-dimension Spearman. **Descriptive**, and the one column that is unfair to Harmony by construction — which is why it is one column and not the table. Harmony is the floor here, never the thing being ranked; [04_9's header](../../04_drvi_epithelial/04_9_embedding_control/compare_embeddings_epi.py) states that argument in full |

## The nulls, which are the point

Every number above is positive for almost any vector on almost any graph, so **none of them
means anything alone**. Three reference levels are computed on the same graph, the same folds
and the same code path.

**1. Permutation within cohort** (`--n-perm`, default 20). The score shuffled among the cells
of its own patient: keeps the per-cohort mean, destroys the cell-level structure. This is the
floor the graph itself induces, cohort composition included. In practice it lands on
Moran's I ≈ 0, which is what makes it a sanity check rather than a result.

**2. Size- and expression-matched random gene sets** (`--n-random`, default 5 per signature).
For each signature, random sets of the **same size** drawn from the **same expression bins**
`sc.tl.score_genes` uses, scored by the same call and standardised the same way. A signature of
highly expressed genes scores differently from one of rare genes whatever the biology, so
matching bin by bin is what makes this the level of *"a gene set like this one, with no
meaning"* — and it is the number "present" has to beat. Cached in
`$DATA_DIR/05_tum/random_signature_scores_<collection>_<run_id>.csv`; `--n-random 0` skips it
and `--overwrite` redraws it.

**3. The confounder arm.** `n_genes_by_counts` and `pct_counts_mt`, z-scored within cohort,
carried through the identical code path as readouts named `__depth` and `__mito`. This one is
not optional here: **05_3b found that Harmony *concentrates* depth rather than removing it**
(ρ(n_genes) on H2 grows from +0.636 to +0.674 under the correction), so a signature that does
not clear the `__depth` row *in the Harmony space* is a depth readout in the Harmony space,
whatever its name is. The figures draw those two rows hollow and hatched so they cannot be
read as results.

## Three spaces, one set of cells

| arm | what it is | what it answers |
|---|---|---|
| `harmony` | 05_3b's 64 corrected components | the thing the step was asked about |
| `drvi` | 05_3's latent space, the phase's own | the reference level |
| `pca` | the **uncorrected** PCA Harmony started from (`pca_<harmony_run>.npy`) | did the correction cost the programme anything |

The cells and the scores are **identical** in all three. The spaces are realigned on cell
**name**, never on position — the same guard 05_3b's own DRVI comparison uses — and the step
aborts rather than proceed if any arm is missing a scored cell. The `pca` arm is the one whose
array carries no cell names of its own; it takes the Harmony embedding's index and is
**dropped**, not lined up on a guess, if the two disagree in length.

The *k*-NN graph is built the way 05_3 and 05_3b build theirs — `k = 15`, the whole space as
the representation, **no rescaling of the coordinates** — so the only thing that differs
between the arms is the space, and "neighbours in the space" means in this table what it means
in those two steps' UMAPs. The graph Moran's I runs on is binary and symmetric rather than
UMAP-weighted, deliberately: a weighting scheme is a second thing that could differ between the
arms and there is nothing to gain from it here.

## What this step does not do

It does not name a Harmony dimension, does not rank Harmony against DRVI as integration
methods — that is phase 02, on what they both promise — and does not write anything into the
per-cell tables 05_6 owns. It **re-scores nothing**: the real signatures come from 05_6's
cached csv, and `sc.tl.score_genes` is called only for the random null.

## What comes out

`../tables/<collection>/<run_id>/`:

```
signature_presence_<collection>_<run_id>.csv   one row per (space, readout), every column above
space_summary_<collection>_<run_id>.csv        one row per space: dimensions, effective rank
random_null_<collection>_<run_id>.csv          which random draw was matched to which signature
```

The run id in those names is the **DRVI** one, and that is not a mistake: it names the *cells*
and the *scores*, which 05_6 owns and which are the same in every arm. The `space` column says
which coordinate system each row was measured in, and the file header says so too — this is the
one table of the phase that `signature_common.write_table` does not stamp, precisely because
that helper writes "*\<Space\> run \<run_id\>*", which is right for a table measured in one
space and wrong for this one.

Figures in `../figures/05_9_embedding_control/<collection>/<run_id>/`:

| figure | |
|---|---|
| `morans_i_by_space` | the headline: do neighbours agree, per space, with the matched-random level ticked on each bar |
| `knn_r2_cohort_by_space` | the cross-patient version of the same |
| `presence_over_null` | observed Moran's I minus the matched-random mean, **in units of the random sd** — the figure that answers the question as asked, with guides at 0 and 3 sd |
| `knn_vs_ridge` | how much of it is linear, one panel per space, the diagonal drawn |
| `umap_top_signatures` | the eyeball check: the signatures the Harmony space predicts best, on each space's own UMAP |

## Reading it

Three traps, in the order they catch people.

**A high Moran's I is not evidence on its own.** Read it against the matched-random level on
the same bar and against `__depth` in the same space. `presence_over_null` is the figure that
does the first of those for you.

**`nonlinear_gap` is not a quality score.** On a corrected PCA a negative gap is the expected
sign and says the programme is a clean gradient; on DRVI a positive one says the programme sits
somewhere a linear read of the coordinates misses. Neither is better. `space_r2_cohort` is the
column that does not care which of the two it was.

**`max_abs_rho` is the only column that ranks the methods, and it ranks them on DRVI's own
promise.** Harmony has never claimed axis-level interpretability — it claims batch correction
with biological signal preserved, and phase 02 is where it was benchmarked on that. A low
`max_abs_rho` for Harmony next to a high `space_r2_cohort` is the useful statement of this
step, and it is a statement about **where** the programme sits, never about whether Harmony is
the worse method.
