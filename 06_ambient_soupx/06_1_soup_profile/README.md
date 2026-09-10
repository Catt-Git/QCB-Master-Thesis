# 06_1 - the ambient profile of each channel

The only step of phase 06 that runs on the cluster, and the reason is one file:
`raw_feature_bc_matrix.h5`.

SoupX corrects a cell against the *composition of the soup* - which gene contributes what
share of the cell-free mRNA in that emulsion. That quantity cannot be read off the cells; it
has to be measured in the droplets that contain no cell. `shiao.h5ad` was built from
`filtered_feature_bc_matrix.h5` (`00_6_build_h5ad/build_combined_h5ad.py`), so those droplets
were discarded at the very first step of the thesis and only exist in the Cell Ranger output
that is still on the cluster next to the FASTQ.

What this step is, therefore, is a **reduction**: several million droplets per channel go in,
one table of ~33,000 numbers comes out, and only the tables travel.

```
$CELLRANGER_DIR/<sample>/outs/raw_feature_bc_matrix.h5    0.5-2 GB, on the cluster
                          |
                          v
$DATA_DIR/06_amb/soup/<sample>.csv.gz                     ~300 KB, copied down
```

## Which droplets count as empty

Those with **1 to 100 UMIs**, which is SoupX's own `soupRange = c(0, 100)`. Reproducing the
package's rule rather than inventing one keeps the profile consistent with what `autoEstCont`
assumes about it downstream. The upper bound is the method's assumption about where cells
stop, and it is the parameter of this phase most worth varying if a channel looks odd
(`--max-umi`).

The empty droplets are deliberately **not** intersected with "barcodes Cell Ranger did not
call": any called cell has far more than 100 UMIs, so the ceiling excludes them by
construction, and the rule then works whether or not the filtered matrix is at hand.

## Genes

`var_names_make_unique()` on the same GRCh38-3.0.0 reference `build_combined_h5ad.py` used,
so the symbols here and the symbols in `shiao.h5ad` are the same vocabulary. They are not the
same *length* - phase 01 filtered genes down to 30,869 - and the alignment happens later, in
`06_2/prepare_soupx_input.py`, where the object's gene order is known and `est` is
renormalised over the genes that survive.

## Running it

```bash
cd 06_1_soup_profile && mkdir -p logs
export DATA_DIR=/users/genomics/albertoc/Tesi/hopes_and_dreams/datasets
export CELLRANGER_DIR=/users/genomics/albertoc/Tesi/hopes_and_dreams/cellranger_out
sbatch --export=ALL,DATA_DIR=$DATA_DIR,CELLRANGER_DIR=$CELLRANGER_DIR submit_soup_profile.slurm
```

**Set `CELLRANGER_DIR`.** The default is `$DATA_DIR/cellranger_out`, which is what
`00_5/run_cellranger_batch.sh` documents as its `OUTDIR` - but on the cluster this thesis ran
on, the batch landed one level up, as a **sibling** of `datasets/` and not inside it. The
script stops immediately with that path in the message rather than failing 100 channels in.

`normal` partition, 16 GB, one channel at a time, no `--time` - the repo's SLURM convention.
The memory is smaller than the step sounds: a `raw_feature_bc_matrix.h5` here is ~6 MB on
disk and a few hundred MB expanded, because the 737k-barcode whitelist is almost entirely
empty. Resuming is the default: a channel whose `.csv.gz` exists is skipped, `--force`
recomputes.

**Check the directory names before launching.** The profiles are joined to the channels **by
name**, in `06_2/prepare_soupx_input.py`, so a `--id` that differs from `obs['sample']`
produces a profile nothing can use - and it surfaces locally, after the job:

```bash
ls -d $CELLRANGER_DIR/*/outs/raw_feature_bc_matrix.h5 | awk -F/ '{print $(NF-2)}' | sort -u > /tmp/cr.txt
tail -n +2 ../../00_raw_data_processing/00_6_build_h5ad/sample_metadata_final.tsv | cut -f1 | sort -u > /tmp/obj.txt
comm -13 /tmp/cr.txt /tmp/obj.txt
```

That table has 150 sample ids against the object's 148: `P30_C_P` and `P57_A_P` never made it
into `shiao.h5ad`, so seeing those two - and nothing else - is the expected output. Extra
profiles are harmless either way; `prepare_soupx_input.py` iterates the object's channels, not
the ones it finds on disk.

Then bring back the two things 06_2 needs:

```bash
rsync -av <cluster>:$DATA_DIR/06_amb/soup/           $LOCAL/datasets/06_amb/soup/
rsync -av <cluster>:$DATA_DIR/06_amb/soup_census.csv $LOCAL/datasets/06_amb/
```

`soup_census.csv` is worth reading before moving on: it records, per channel, how many
droplets were in the soup range, how many UMIs they carried, and which gene dominates. A
channel missing from it - a Cell Ranger run that never finished, a matrix that was cleaned up -
makes `06_2` stop rather than silently correct 147 channels out of 148.
