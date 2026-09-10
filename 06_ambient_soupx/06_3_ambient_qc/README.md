# 06_3 - the ambient call

One notebook, `ambient_qc.ipynb`, and one output: `$DATA_DIR/06_amb/ambient_keep.csv`, a
boolean per cell that `06_4` applies.

This is the decision step of the phase. `06_2` produced a continuous quantity; where to cut
it is chosen by looking at a distribution, which is why it is a notebook and not a script -
exactly as `05_1/call_malignant.ipynb` is to `run_infercnv.R`.

## The axis

SoupX estimates one contamination fraction `rho` per **channel** - it is a property of the
emulsion - so `soupx_rho` is constant within a sample and cannot separate cells. What varies
cell by cell is how much of that channel's soup actually landed in a given droplet, and the
correction measures it directly:

```
soupx_frac_removed = 1 - umis_after / umis_before
```

A droplet whose profile was mostly the channel average loses most of its counts; a cell with
a strong specific transcriptome loses a few percent.

## What the notebook shows, in order

1. **rho per channel** - the level, the spread, and which channels were never actually
   estimated (`method = 'fallback_prior'`, drawn apart).
2. **`soupx_frac_removed` per cell** - the distribution, the same split by cell type, and the
   same against library size. The third panel is the control question: if the tail is only
   small cells, the threshold is a library-size filter under another name and `01_2` already
   has one.
3. **The phase-01 UMAP coloured by ambient burden** - which parts of the map phase 01 drew
   were being drawn by soup. A landmark; no claim rests on it.
4. **Cross-lineage markers, before and after** - the only panel that is evidence the
   correction *worked* rather than merely subtracted. *PTPRC*/*CD3D* in epithelial cells and
   *EPCAM*/*KRT18* in immune cells should fall, while each lineage's own markers stay put.
   A correction that also flattens *EPCAM* in epithelium is removing signal, not soup.
5. **The threshold**, then what it costs by cell type and by cohort. The per-cohort panel
   draws `05_2`'s `MIN_CELLS_PER_COHORT = 200` line, because a cohort losing cells here can
   fall below it later.

## The cell to edit

```python
MAX_FRAC_REMOVED = 0.50   # a cell that lost more than half its counts to the soup
MIN_UMIS_AFTER = 200      # ...or has fewer than this many left
MIN_GENES_AFTER = 100     # ...or expresses fewer than this many genes
DROP_FALLBACK_CHANNELS = False
```

The two floors exist because a cell can pass the fraction test and still be a 40-UMI droplet
afterwards: phase 01 filtered on the **uncorrected** totals, and this is the first chance to
filter on the real ones. `MIN_GENES_AFTER` is deliberately looser than `05_2`'s own
`MIN_GENES = 100` applied to its subset, because dropping a cell here also removes it from
`05_1`'s inferCNV run, which is a stronger action.

Setting `MAX_FRAC_REMOVED = 1.0` with both floors at 0 keeps every cell: the counts are still
corrected and the phase becomes a pure decontamination with no cell filter. If the tail turns
out to be thin, that is the correct outcome and should be recorded as one.
