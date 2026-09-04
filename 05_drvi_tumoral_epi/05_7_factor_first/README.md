# 05_7_factor_first — Route B

The dimension defines itself; prior knowledge only names it. Counterpart of
[04_6_factor_first](../../04_drvi_epithelial/04_6_factor_first/), on this phase's run.

**The unit of analysis is the gene.** The question is what gene programme a dimension encodes,
irrespective of what anyone was looking for. This is the only **discovery** route in the stage:
it can name dimensions nobody asked about, and it is the corrective for the confirmation bias
built into Route A, which by construction can only find states brought in from outside.

It is also the route that cashes in the property DRVI was chosen for over the higher-scoring
methods of the phase-02 benchmark: the additive decoder gives every dimension a directly
readable gene-level footprint. Without this stage any integration method would have done.

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
python3 factor_first_tum.py                      # scie, the default
python3 factor_first_tum.py --collection emt     # the EMT lists
python3 factor_first_tum.py --n-top-genes 500    # a deeper list
python3 factor_first_tum.py --no-hallmark        # custom signatures only, fully offline
```

Needs `05_4` (the `.gmt`), `05_3` (the embedding) and `05_6` (the row order). No GPU, no model:
the genes × dimension-direction table is rebuilt from `embed.varm` alone.

## It carries more weight here than in 04

This phase has **no non-constant biological annotation of its own**. `cell_type` is `malignant`
for every cell; everything else is either a normal-breast landmark (`cell_type_01_4`) or a
clustering (`optscib_tum_leiden`). Naming a dimension **by its genes** is therefore not one of
two equivalent ways to read the latent space — it is the one that owes nothing to a borrowed
vocabulary. The groupings used in 05_3 and 05_6 are the cross-check on this route, not the
evidence.

## The method, and the three things that are not defaults

- **Both directions of each dimension are tested separately.** DRVI can encode two distinct
  concepts on the two sides of one axis, and pooling them cancels the programmes against each
  other.
- **Every test is offline and hypergeometric against a declared background** (`gp.enrich`), never
  against Enrichr's implicit all-human-genes universe.
- **Benjamini–Hochberg is applied once across every dimension-direction × gene-set pair**, not
  per query as Enrichr does — with dozens of directions a per-query FDR is far too permissive.
  Pairs that returned no overlap are p = 1 and cannot become significant, but they are added back
  into the denominator. The correction is computed **inside one collection**, so adding the EMT
  lists cannot move a SCIE p-value.

The correction is *conservative in count and still not a set of independent tests*: the block
structure of 05_4's Jaccard matrix is exactly the dependence BH assumes away. Both facts are
printed by the script.

## The background is the HVG panel, and it is not 04's

DRVI trained on the 2,000 batch-aware HVGs of **05_2**, re-selected inside the tumour, so a gene
outside that set could never have entered a top-gene list and the ORA background is the
**training feature set**, not the transcriptome. Each signature is effectively tested in its
HVG-restricted form.

Which part of a list that is has **moved against 04** — `ESC_WONG` 48 → 64 genes, `LIM_STEM`
160 → 143 — so this route's power per signature is not the same as in 04 in either direction.
05_4's `n_in_hvg_background` column is the one to read, from this phase's table.

## Half the denominator of 04

64 dimension-directions here (2 × 32) against 128 there (2 × 64). The same p-value therefore goes
further. That is a property of the latent size, not evidence, and it is one more reason 05_8
refuses to promote anything on this route alone.

Nothing is pruned: all 2 × `n_latent` directions are tested and the vanished flag is reported
rather than acted on.

## What comes out

```
../tables/<collection>/dim_geneset_signed_significance_*.csv   the matrix 05_8 joins
../tables/<collection>/factor_first_significant_*.csv          every significant pair
../tables/<collection>/factor_first_hallmark_significant_*.csv the sanity-check library
$DATA_DIR/05_tum/factor_first_top<N>_<collection>_<run_id>.tsv every pair tested, unfiltered
$DATA_DIR/05_tum/factor_first_top<N>_genes_<run_id>.tsv        the top-gene lists themselves
```

The signed matrix is on the **same row order as Route A**, so the two heatmaps can be read side
by side in 05_8. The sign is the direction of the axis carrying the enrichment.

Hallmark 2020 is kept alongside the collection as a sanity-check library — a dimension that
enriches for nothing in the custom lists can still be named. It is fetched once and cached, so a
re-run needs no network; `--no-hallmark` skips it entirely.

**How this route fails:** gene-level enrichment says nothing about cells. A dimension can be
strongly enriched for a stemness list while being driven by a handful of cells, by a
patient-specific effect that shares part of the gene set, or by a dissociation-stress response
overlapping the same genes. The top-gene list is also a truncation, so long diluted signatures
systematically under-enrich relative to short sharp ones. That is why 05_8 exists.
