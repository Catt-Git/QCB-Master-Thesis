# 05_2_subsetting

From `shiao.h5ad` to the object the rest of phase 05 runs on. Same shape as
[04_1_subsetting](../../04_drvi_epithelial/04_1_subsetting/), with one difference that runs
through everything here: **04 defined its subset by a label, this one defines it by a
measurement.** The epithelial compartment of 04 is whatever CellTypist called epithelial; the
malignant compartment here is whatever inferCNV found aneuploid, and the CellTypist labels come
along as metadata rather than as the selection criterion.

A full re-run of the phase-01 pre-processing on the subset, for the reason that applies at
every subsetting step in this repo: scran size factors, HVGs, PCA, neighbours, UMAP and leiden
were all computed on populations of which this one is a minority (5.8% of `shiao.h5ad`, 48% of
04's epithelium), and none of them describe these cells. Only the raw counts and the metadata
carry over.

## The cell set

`CELL_SET` picks what gets subsetted, and with it the file prefix, so both can live in
`$DATA_DIR/05_tum/` without overwriting each other:

| `CELL_SET` | cells kept | prefix |
|---|---|---|
| `tum` (default) | `cnv_status == 'malignant'` | `shiao_tum_*` |
| `epi` | every epithelial cell under the post-CNV labels, malignant included, `not_tested` dropped | `shiao_epicnv_*` |

The mapping lives in `cell_set.py` and nowhere else — four scripts write into one directory,
and a prefix computed independently in each of them is a silent-overwrite bug waiting to
happen. `cell_set.py` also holds every threshold, with the argument for each in a comment.

`cell_set.py` is also where the other two knobs of the phase live, for the same reason: `HVG_SET`
(below) and `N_LATENT`, which is not a property of the cells at all but is the third segment of
the run id — and `run_id()` there is the only place that string is built. 05_3 takes `N_LATENT`
as the default of `--n-latent`; 05_4 - 05_8 read it and have no flag of their own. See
[the phase README](../README.md#which-run-05_4---05_8-read).

Resuming is per cell set, because the outputs are named differently: running `epi` after `tum`
re-runs all four steps, as it should.

## Execution order

| # | File | What it does | Where |
|---|------|--------------|-------|
| 1 | `subset_and_qc.ipynb` | Join `cnv_status.csv` and `cell_annotation_cnv.csv` onto `shiao.h5ad`; subset; write `compartment`, the post-CNV `cell_type` and the pre-CNV `cell_type_01_4`; restore `.X` = raw counts; drop every phase-01 derived slot and the stale palettes; recompute QC and re-apply the cell thresholds; gene filter `min_cells=3`; cohort census and drops. | local (notebook) |
| 2 | `subsetting_all.sh` | Driver for steps 3-6: runs them in sequence, resumes from the last completed one, logs to `logs/`. | local |
| 3 | `scran_norm_tum.py` | scran re-estimated on the subset (`quickCluster → computeSumFactors → logNormCounts`); adds `size_factors`, clears `.raw`, casts `.X` to float32. | local |
| 4 | `cell_cycle_score_tum.py` | Tirosh/Regev re-scoring → `S_score`, `G2M_score`, `phase`. | local |
| 5 | `reduce_data_tum.py` | Batch-aware HVG (`cohort`, 2000) + PCA(50) + neighbours + UMAP; writes the HVG list **and** the 2,000-gene DRVI input, plus the overlap with the 01_5, 03_1 and 04_1 HVG lists. | local |
| 6 | `clustering_tum.py` | Leiden optimal-resolution sweep (0.1-2.0) → `shiao_tum.h5ad`. | local |
| 7 | `visualization_tum.ipynb` | Figures on the final object. Read-only. | local (notebook) |
| — | `hvg_no_mt.py` | **Off-chain, on demand.** Rebuilds step 5's panel without the 11 `MT-` genes, refilling to 2,000 with the next 11 of the same ranking → `<prefix>_hvg_2k_nomt_list.csv` + `<prefix>_hvg_2k_nomt.h5ad`. See *The MT-free panel* below. | local |

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   run subset_and_qc.ipynb first
./subsetting_all.sh                 # every step, resuming
CELL_SET=epi ./subsetting_all.sh    # the same chain on the control set
./subsetting_all.sh --dry-run       # print what would run
./subsetting_all.sh reduce cluster  # only the named step(s)
```

> **Note on step 2.** The scripts are duplicated from 04_1 (which duplicated 03_1, which
> duplicated 01_3/01_4/01_5) rather than called with different paths, so this phase reads as a
> self-contained Materials & Methods section and the earlier scripts stay frozen. The
> parameters are identical except where stated below.

> **Note on where.** ~36k cells, half of 04: the whole chain runs locally in sequence and there
> is no SLURM wrapper.

## The three deliberate differences from 04_1

**The NMI target of the leiden sweep.** 04_1 maximises NMI against `cell_type`; on this subset
that column is the constant `malignant`, and NMI against a constant is zero at every
resolution. `clustering_tum.py` detects that and falls back to `cell_type_01_4`, the pre-CNV
CellTypist label, which is not constant inside the tumour. The full argument, including why the
alternative (an internal criterion like modularity) was not taken, is in that script's
docstring. The column that produced the profile is recorded in
`uns['optscib_tum_leiden_label_key']`, so a reader never has to guess which one was used.

**No treatment-completeness filter.** 04_1 dropped cohorts missing any of BASE/PD1/RTPD1
because a treatment phase was planned after it. That phase is not planned any more, and on the
malignant subset the filter is actively harmful: several cohorts clear it with two cells in a
timepoint, while Patient63 would lose 1,673 malignant cells for having none.
`DROP_INCOMPLETE_COHORTS = False`.

**`not_tested` cells are excluded, not assumed normal.** Under `CELL_SET=epi` this is an actual
filter: a cell no inferCNV run covered has no evidence either way, and calling it normal
epithelium is the mistake phase 05 exists to stop making. Under `tum` it happens by
construction. Either way the post-condition asserts no `not_tested` cell survived.

## What comes out

19 cohorts, ~36,192 cells before the QC re-application. `$DATA_DIR/05_tum/`:

```
<prefix>_raw.h5ad                     the notebook's output, raw counts
<prefix>_norm.h5ad                    + scran
<prefix>_norm_cc.h5ad                 + cell cycle
<prefix>_reduced.h5ad                 + HVG / PCA / neighbours / UMAP
<prefix>_hvg_2k_list.csv              the 2,000 selected genes
<prefix>_hvg_2k.h5ad                  the DRVI input for 05_3
<prefix>_hvg_2k_nomt_list.csv         the same panel without the MT- genes  (hvg_no_mt.py)
<prefix>_hvg_2k_nomt.h5ad             the DRVI input for 05_3 under HVG_SET=nomt
<prefix>.h5ad                         + leiden - the definitive object
<prefix>_leiden_resolution_profile.csv
```

Figures go to `../figures/05_2_*/`.


## The MT-free panel (`HVG_SET=nomt`)

Step 5 selects 2,000 batch-aware HVGs and, on the `epi` object, **11 of them are
mitochondrial**: `MT-ND1`, `MT-ND2`, `MT-CO1`, `MT-CO2`, `MT-ATP8`, `MT-CO3`, `MT-ND3`,
`MT-ND4L`, `MT-ND4`, `MT-ND6`, `MT-CYB`. What they vary with is the mitochondrial fraction of a
cell's RNA - dissociation stress, membrane damage, ambient lysate - which is the quantity
`MAX_PCT_MT` already filters cells on. That is QC, not epithelial state, and a DRVI latent
dimension answering to it is a dimension not spent on biology.

`hvg_no_mt.py` drops them and refills the panel to 2,000 with **the next 11 genes of the same
ranking**: `PTPRB`, `SVOPL`, `TMEM47`, `TGFB3`, `FLRT2`, `KCNA1`, `AC091057.3`, `HIST2H2AB`,
`RHBDL1`, `CILP`, `AC009041.2`, all highly variable in 6 cohorts with normalized dispersion
0.502-0.509 against the 0.510 of the last gene inside the cut. The other **1,989 genes are gene
for gene** the ones already in `<prefix>_hvg_2k_list.csv` - the script reproduces scib's
selection over the *unchanged* `sc.pp.highly_variable_genes(n_top_genes=2000)` call and asserts
the reproduction against the file on disk before touching anything, because
`hvg_batch(target_genes=2011)` would feed 2011 into scanpy and reshuffle the ranking rather than
extend it.

The metallothioneins `MT1A/E/F/G/H/M/X` and `MT2A`, and `MTRNR2L12`, share the prefix and
nothing else: nuclear genes, real epithelial biology, kept. The match is the exact prefix `MT-`.

Nothing is overwritten. The outputs carry the `_nomt` tag `cell_set.py` derives from `$HVG_SET`,
and so does the 05_3 run id built on them (`drvi_epicnv_64_nomt` against `drvi_epicnv_64`), so
the two runs sit side by side and the difference between them is readable. `<prefix>.h5ad` and
its `optscib_tum_leiden` clustering are **not** rebuilt - leiden stays the partition both runs
are read against.

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
CELL_SET=epi python3 hvg_no_mt.py          # ~10 min: the HVG scan, then PCA/neighbours/UMAP
CELL_SET=epi python3 hvg_no_mt.py --force  # rebuild over an existing panel
```
