# 05_5_cytotrace2

CytoTRACE2 potency per cell, one patient at a time. Counterpart of
[04_4_cytotrace2](../../04_drvi_epithelial/04_4_cytotrace2/), on the malignant subset.

CytoTRACE2 is a per-cell predictor of differentiation potency, not a gene set. It belongs to
Route A only and has **no Route B counterpart**. Its value here is that it is independent of the
collaborator's lists, so it is the only non-circular evidence on the stemness axis — which
matters, because the stemness consensus was deliberately not built (it captured proliferation
only) and no single stemness signature is primary.

```bash
conda activate cytotrace2-py          # NOT benchmark-py-r, see below
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
python3 cytotrace2_tum.py --dry-run   # export the per-patient matrices, do not score
python3 cytotrace2_tum.py             # export + score + concatenate
python3 cytotrace2_tum.py --cores 4   # cap the cores
```

## Its own environment, and that is not a preference

`cytotrace2-py` 1.1.0.4 declares `numpy<2.0.0` as a hard pin. Installing it into
`benchmark-py-r` does not add a package, it rolls the whole stack back (numpy 2.4 → 1.26,
pandas 3.0 → 2.3, scanpy 1.12 → 1.11, anndata 0.13 → 0.12) and leaves two packages with
unsatisfiable requirements. That was done once, in 04, and undone. Use
`environments/cytotrace2-py.yml`. The script checks for the package and stops with that
instruction rather than failing halfway.

Everything else in this stage runs without it: `cell_first_tum.py` simply drops the CytoTRACE2
quadrant definition if the csv is missing, and says so — at the cost of leaving the stemness
axis resting entirely on the lists.

## Why per patient

The score is computed **within the object it is given**, so scoring all 19 cohorts at once would
let the ranking absorb batch: a patient sequenced deeper would come out more potent than one
sequenced shallow, and "stemness" would be a batch axis. Each cohort is scored on its own and
the results concatenated, so the output cannot be confounded with `cohort` by construction.

The cost is that scores are comparable **within** a patient and only ordinally across patients —
which is exactly how Route A uses them, since 05_6 standardises every readout within `cohort`
anyway.

## What changes against 04

**What it buys.** In 04 a high potency score could be a normal progenitor sitting in the
epithelial compartment, and nothing in that phase could tell. Here every cell is aneuploid, so
that reading is gone.

**What it does not buy.** A high score can still be a *cycling* cell, and that risk is larger
here, not smaller: `Lumsec-prol` is **51%** of this subset against 25% of 04's epithelium.

**The summary tables therefore changed shape.** 04 prints mean potency per `cell_type`; that
column is the constant `malignant` here, so the script summarises per grouping instead — the
leiden partition of 05_2, then the pre-CNV CellTypist label — and adds `phase`, because the
question this table is asked to answer is a cycle question and not a lineage one.

> **How to read it.** Potency should not rank the groups in the same order the cycle does. If
> the top of the leiden table is exactly the S/G2M clusters and the top of the label table is
> `Lumsec-prol`, then potency is tracking proliferation, and the stemness axis of Route A rests
> on a cycle readout. 05_6 recomputes its whole target region inside G1 for the same reason;
> this is the earlier and much cheaper warning.

## What comes out

```
$DATA_DIR/05_tum/cytotrace2_<run_id>.csv          per cell: score, potency, relative, preKNN
$DATA_DIR/05_tum/cytotrace2_work_<run_id>/        per-patient matrices and raw results
../tables/scie/cytotrace2_by_<grouping>_scie_<run_id>.csv
```

The per-cell csv is **not** collection-scoped: it is a measurement. The summary tables go under
`scie`, the collection that declares CytoTRACE2 as a readout. `--write-obs` also adds the
columns to `shiao_tum_<run_id>.h5ad`; off by default, because that object is 05_3's output.

The per-patient count matrices are deleted afterwards unless `--keep-work` — they are plain-text
gene x cell tables and they are large.
