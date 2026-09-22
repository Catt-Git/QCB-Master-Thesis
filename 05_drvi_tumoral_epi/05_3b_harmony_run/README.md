# 05_3b_harmony_run

Harmony at the same step as `05_3_drvi_run`, on the same cells and the same 2,000 genes.

**This is a control on 05_3, not a second main line, and it is a deliberate dead end.** It
writes an embedding, figures and two diagnostic tables, and stops. 05_4 – 05_8 are not run on
it, and `--downstream` (the object carrying `obsm['X_harmony']`) is off by default, because
what those steps do — read a latent dimension and ask which gene programme loads on it — is a
question about a DRVI *decoder*. Harmony has no decoder: it rotates and shifts a PCA until the
cohorts overlap, so its dimensions are corrected principal components. They have loadings, but
the loadings belong to the PCA *before* the correction, and the correction is a different
affine map per cohort. Reading a gene programme off one of them means reading the loadings of a
basis the coordinates are no longer expressed in.

## What it is for

The question phase 05 keeps running into: **is the structure DRVI finds in these cells a
property of the cells, or a property of DRVI?** A linear method with a different objective, on
the same input, either recovers the same axes or does not — and either answer is worth having.

Three things come out:

| | |
|---|---|
| `figures/05_3b_<run_id>/umap_*.png` | the same keys, palettes and seeded permutation as 05_3, so the only thing that differs from `figures/05_3_drvi_tum_64_nomt/` is the space |
| `figures/05_3b_<run_id>/umap_vs_drvi_*.png` | the two spaces side by side on `cohort` (did the batches mix), `optscib_tum_leiden` and `cell_type_01_4` (did the states survive it), `cnv_score` (is either space just drawing aneuploidy) |
| `tables/05_3b_<run_id>/*.tsv` | the two diagnostic tables below — the point of running this at all |

### `harmony_dim_stats.tsv`

One row per corrected dimension: its variance, the **η²(cohort)** before and after the
correction, and the **Spearman correlation with sequencing depth, mito/ribo fraction and CNV
burden**, again before and after. This is the 05 diagnosis — the one that found DR 1 of the
DRVI space to be depth and DR 10 to be a real axis — run on Harmony's dimensions instead. The
before/after pairing is what makes it readable: a correlation present in both columns is
something Harmony left alone, one that grows is something the correction *introduced*.

### `harmony_vs_drvi_corr.tsv` and `harmony_vs_drvi_best_match.tsv`

The 64 × 64 Spearman matrix between these dimensions and the DRVI ones (realigned on cells **by
name**, never by position), and for each DRVI dimension the Harmony dimension matching it best.
A DRVI axis with a strong linear counterpart is an axis a completely different method also
found. The converse does **not** follow: a DRVI dimension with no match may be non-linear
structure Harmony cannot represent, or it may be nothing. The table narrows the question, it
does not close it. `vanished` is carried over from the DRVI embedding so the dimensions DRVI
itself did not use can be dropped from the reading.

## The two parameters that are choices, not defaults

**64 components, not scanpy's 50.** Harmony has no latent size of its own — it corrects a PCA
and returns a matrix of the same width. 64 equals the DRVI run it is compared against
(`drvi_tum_64_nomt`), so the two spaces have the same number of coordinates and the matrix
above is square. It is *not* what `02_2_integration` uses: that runs
`scib.integration.harmony`, which takes scanpy's default 50. **So this run is comparable with
05_3 and not with the phase-02 benchmark.** `--n-comps 50` gets the other one.

**The genes are z-scored before the PCA.** `.X` of the 05_2 object is scran log-normalised
counts on 2,000 HVGs. Harmony's own recipe — Seurat's, and every harmony-pytorch example —
scales before the PCA so a few high-variance genes do not own the first components;
`scib.integration.harmony` does not, because scib controls scaling as a separate variant of the
input. This script scales by default because it runs the method the way the method is meant to
be run. `--no-scale` turns it off. **The flag is in the log and nowhere in the run id**, so
running both needs `--fig-dir`, or the second overwrites the first.

## Run

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
cd 05_drvi_tumoral_epi/05_3b_harmony_run

N_LATENT=64 conda run -n benchmark-py-r --no-capture-output python3 run_harmony_tum.py
```

`benchmark-py-r` is the environment: it is the only one holding both `harmony` (harmony-pytorch,
the implementation scib wraps) and `drvi`, and the second is needed only to read the DRVI
embedding's `var` table for the comparison. ~4 minutes on 42,096 cells, GPU if one is there.

`CELL_SET` and `HVG_SET` are honoured exactly as in 05_3, through `05_2/cell_set.py`: unset is
the malignant subset on the panel without the mitochondrial genes, `CELL_SET=epi` the
epithelial control set. The run id is built by the same `cell_set.run_id()`, with
`method='harmony'` — `harmony_tum_64_nomt` beside `drvi_tum_64_nomt`, so neither can land on
the other's files. Resuming is the default; `--overwrite` recomputes.

## What the first run said

`harmony_tum_64_nomt`, 42,096 malignant cells, 19 cohorts, 64 components, scaled.

**The batch correction works, and it is not close.** Mean η²(cohort) over the dimensions goes
from **0.184 on the uncorrected PCA to 0.008** after Harmony; dimensions where the cohort still
explains more than 5% of the variance go from **38/64 to 1/64**, and the worst one from 0.847
to 0.079. Whatever else is below, Harmony did the job it is for.

**Sequencing depth is not batch, and Harmony does not touch it — it concentrates it.** The
correlation with `n_genes_by_counts` *grows* under the correction on the dimensions that carry
it:

| dim | ρ(n_genes) before | after |
|---|---|---|
| H2 | +0.636 | **+0.674** |
| H4 | −0.403 | **−0.621** |
| H0 | +0.281 | **+0.443** |
| H6 | +0.126 | **+0.407** |

This is the expected direction and worth stating plainly: Harmony corrects **cohort**, and
depth varies inside a cohort as much as between cohorts. Removing the between-cohort variance
leaves the within-cohort depth gradient a larger share of what is left.

**DRVI spreads the same confounder; Harmony concentrates it.** In the DRVI space the depth
signal sits on many dimensions at moderate strength — DR 2 (ρ = +0.52), DR 1 (−0.39), DR 63
(−0.39), DR 62 (+0.31), DR 8 (+0.31), DR 33 (−0.31) — where Harmony puts most of it on two
dimensions at ρ ≈ 0.62–0.67. Neither is obviously better, and they fail differently: Harmony's
depth is *droppable* (ignore H2 and H4), DRVI's is diffuse and cannot be excised by dropping
dimensions.

**The two spaces agree much less than they disagree.** Of the 56 non-vanished DRVI dimensions,
only **4 have a Harmony counterpart at |ρ| > 0.5** and 28 at > 0.3. The strongest pairings are
DR 3 ↔ H1 (+0.73), DR 5 ↔ H1 (+0.62), DR 8 ↔ H2 (+0.61), DR 2 ↔ H4 (−0.57) — and two of those
four (DR 8, DR 2) are the depth dimensions, i.e. the axes the two methods agree on most are
partly the confounder they share.

**DR 10 has no linear counterpart.** Its best match is H14 at |ρ| = 0.30. DR 1's is H8 at 0.38.
Under the reading above this is a narrowing and not a verdict: DR 10 being invisible to a
linear method is consistent with it being non-linear structure only a decoder can see, and
equally consistent with it being nothing. What it does rule out is the easy version of the
claim — DR 10 is not simply a rotation of a principal component.

## Files

```
05_3b_harmony_run/
  run_harmony_tum.py                  the whole step
  README.md                           this file
  logs/harmony_tum_64_nomt.log        the run above

$DATA_DIR/05_tum/
  embed_harmony_<run_id>.h5ad         corrected embedding, one var per dimension, UMAP
  pca_<run_id>.npy                    the UNCORRECTED PCA, kept only for the before/after
                                      columns of the stats table; safe to delete, the run
                                      recomputes it when it is missing

05_drvi_tumoral_epi/
  figures/05_3b_<run_id>/             the UMAPs, the dimension grid, the heatmaps, the
                                      confound panel and the four vs-DRVI comparisons
  tables/05_3b_<run_id>/              harmony_dim_stats.tsv
                                      harmony_vs_drvi_corr.tsv
                                      harmony_vs_drvi_best_match.tsv
```

`cell_set.run_id()` gained a `method` argument for this step; it defaults to `'drvi'`, so every
run id already on disk and every caller in 05_4 – 05_8 keeps the name it had.
