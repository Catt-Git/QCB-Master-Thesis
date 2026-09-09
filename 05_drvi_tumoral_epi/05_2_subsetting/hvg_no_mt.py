"""
05_2 (variant): the same 2,000-gene panel with the mitochondrial genes taken out.

`reduce_data_tum.py` selects 2,000 batch-aware HVGs and 11 of them, on the CELL_SET=epi
object, are transcripts of the mitochondrial genome (`MT-ND1`, `MT-ND2`, `MT-CO1`, `MT-CO2`,
`MT-ATP8`, `MT-CO3`, `MT-ND3`, `MT-ND4L`, `MT-ND4`, `MT-ND6`, `MT-CYB`). What those genes
vary with is the fraction of a cell's RNA that came from mitochondria - dissociation stress,
membrane damage, ambient lysate - and not a state of the epithelium. They are already what
the `pct_counts_mt <= MAX_PCT_MT` filter is computed on, so the cells that survive QC are the
ones where the signal is *residual*, and a latent dimension spent on it is a dimension not
spent on biology. This script rebuilds the panel without them.

The eight metallothioneins (`MT1A/E/F/G/H/M/X`, `MT2A`) and `MTRNR2L12` share the prefix and
nothing else: they are nuclear genes, they are real epithelial biology (zinc/copper handling,
a well-known stress-response programme), and they STAY. The rule here is the exact prefix
`MT-`, which on GENCODE is the mitochondrial genome and only it.

## How the eleven replacements are chosen

The list has to stay 2,000 genes, so the eleven are refilled from the same ranking that
produced the panel in the first place - the next eleven, nothing hand-picked.

Doing that by calling `scib.preprocessing.hvg_batch(target_genes=2011)` would NOT be it.
`hvg_batch` passes `target_genes` straight into `sc.pp.highly_variable_genes(n_top_genes=...)`,
so a larger target changes which genes are flagged *within each batch*, hence
`highly_variable_nbatches` for every gene, hence the ranking itself: the 2,011-gene answer is
not the 2,000-gene answer plus eleven. So `hvg_batch`'s selection loop is reproduced here
(and checked against the panel on disk, gene for gene) over the *unchanged*
`sc.pp.highly_variable_genes(n_top_genes=2000, batch_key='cohort')` call, then continued past
2,000 until eleven non-`MT-` genes have been added. The 1,989 survivors are therefore
bit-for-bit the ones already in `<prefix>_hvg_2k_list.csv`.

Input : $DATA_DIR/05_tum/<prefix>_norm_cc.h5ad    (what reduce_data_tum.py selected on)
        $DATA_DIR/05_tum/<prefix>_reduced.h5ad    (for obs/uns; the cells and their
                                                   annotation are untouched by this step)
        $DATA_DIR/05_tum/<prefix>_hvg_2k_list.csv (the panel being repaired, as a check)
Output: $DATA_DIR/05_tum/<prefix>_hvg_2k_nomt_list.csv
        $DATA_DIR/05_tum/<prefix>_hvg_2k_nomt.h5ad   <- what 05_3 trains on under HVG_SET=nomt

Nothing already on disk is overwritten: both outputs carry the `_nomt` tag `cell_set.py`
derives from `$HVG_SET`, and so does the 05_3 run id built on top of them
(`drvi_epicnv_64_nomt` against `drvi_epicnv_64`). The definitive object `<prefix>.h5ad` and
its `optscib_tum_leiden` clustering are NOT rebuilt: leiden comes off the unintegrated PCA of
the original panel and stays exactly as it was, which is what keeps the two DRVI runs
comparable against the same partition.

PCA / neighbours / UMAP inside the output are recomputed on the new panel, with the
parameters `scib.preprocessing.reduce_data` uses (50 comps, arpack), so no array in the file
is left describing the old gene set.

Usage:
  export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
  CELL_SET=epi python3 hvg_no_mt.py          # the control set, the one 05_3 retrains on
  CELL_SET=tum python3 hvg_no_mt.py          # the same repair on the malignant subset
  CELL_SET=epi python3 hvg_no_mt.py --force  # rebuild even if the outputs exist

$HVG_SET is set by the script itself: it is what *writes* the `nomt` panel, so reading the
flag would be circular.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

# scib's hvg_batch uses np.in1d, removed in numpy>=2.4 in favour of np.isin. Same shim as
# reduce_data_tum.py, kept because the selection below is scanpy's, called the way scib
# calls it.
if not hasattr(np, "in1d"):
    np.in1d = np.isin

import matplotlib
matplotlib.use("Agg")
import scanpy as sc  # noqa: E402

import cell_set as C  # noqa: E402

# The prefix of the mitochondrial genome on GENCODE. Exact, and deliberately not a regex:
# MT1A / MT2A / MTRNR2L12 are nuclear and must not match.
MT_PREFIX = "MT-"


def ranked_hvgs(adata, target: int) -> list[str]:
    """`scib.preprocessing.hvg_batch`'s selection, in order, extendable past N_HVGS.

    The scanpy call underneath is the one reduce_data_tum.py made - `n_top_genes=C.N_HVGS`,
    never `target` - so the ranking is fixed and `target` only says how far down it to read.
    Genes come out in selection order: the full-intersect set by descending normalized
    dispersion first, then the n_batches-1 set, and so on, exactly as scib fills up.
    """
    sc.pp.highly_variable_genes(
        adata,
        flavor="cell_ranger",
        n_top_genes=C.N_HVGS,
        n_bins=20,
        batch_key=C.BATCH_KEY,
    )

    n_batches = len(adata.obs[C.BATCH_KEY].cat.categories)
    disp = adata.var["dispersions_norm"]
    nbatches = adata.var["highly_variable_nbatches"]

    selected: list[str] = []
    for k in range(n_batches, 0, -1):
        tier = disp[nbatches == k].sort_values(ascending=False)
        selected.extend(tier.index[: target - len(selected)])
        if len(selected) >= target:
            break

    if len(selected) < target:
        raise SystemExit(
            f"only {len(selected)} genes are highly variable in at least one cohort; "
            f"cannot reach {target}"
        )
    return selected


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--force", action="store_true",
                    help="rebuild even if the outputs are already on disk")
    args = ap.parse_args()

    C.banner("05_2 HVG panel without the mitochondrial genes")

    in_path = C.path("_norm_cc.h5ad")
    reduced_path = C.path("_reduced.h5ad")
    base_csv = C.path("_hvg_2k_list.csv")

    # The outputs, named by cell_set.py rather than by hand: this script IS the producer of
    # the `nomt` panel, so it sets the flag instead of reading it.
    os.environ["HVG_SET"] = "nomt"
    out_csv = C.hvg_path("_list.csv")
    out_h5ad = C.hvg_path(".h5ad")
    assert C.hvg_tag() == "_nomt", "cell_set.hvg_tag() did not pick up HVG_SET=nomt"

    for p in (in_path, reduced_path, base_csv):
        if not p.exists():
            raise SystemExit(f"missing input: {p}\nRun subsetting_all.sh for this CELL_SET first.")

    if out_h5ad.exists() and out_csv.exists() and not args.force:
        print(f"[have] {out_csv}", flush=True)
        print(f"[have] {out_h5ad}", flush=True)
        print("Nothing to do; --force rebuilds.", flush=True)
        return

    # ---------------------------------------------------------------- the ranking
    print(f"Reading {in_path}", flush=True)
    adata = sc.read_h5ad(in_path)
    print(f"Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes", flush=True)
    assert C.BATCH_KEY in adata.obs, f"missing batch key {C.BATCH_KEY!r}"

    base = ranked_hvgs(adata, C.N_HVGS)

    # The panel on disk was written in .var order, this one comes out in selection order, so
    # they are compared as sets. A mismatch means the ranking here is not the one that
    # produced the file - in which case the eleven replacements would not be "the next
    # eleven" of anything and there is no point continuing.
    on_disk = set(pd.read_csv(base_csv, header=None)[0].astype(str))
    assert len(on_disk) == C.N_HVGS, f"{base_csv} has {len(on_disk)} genes, expected {C.N_HVGS}"
    missing = on_disk - set(base)
    if missing:
        raise SystemExit(
            f"reproduced panel differs from {base_csv.name} in {len(missing)} genes "
            f"(e.g. {sorted(missing)[:5]}). The scanpy/scib version this runs under is not "
            "the one reduce_data_tum.py ran under; the replacements would not be comparable."
        )
    print(f"Reproduced the {C.N_HVGS}-gene panel of {base_csv.name} exactly", flush=True)

    # ---------------------------------------------------------------- the swap
    mito = [g for g in base if g.startswith(MT_PREFIX)]
    print(f"\nMitochondrial genes in the panel: {len(mito)}", flush=True)
    print("  " + ", ".join(mito), flush=True)
    if not mito:
        raise SystemExit("no MT- genes in the panel: nothing to repair, and no file written.")

    kept = [g for g in base if not g.startswith(MT_PREFIX)]
    target = C.N_HVGS + len(mito)
    extended = ranked_hvgs(adata, target) if target > C.N_HVGS else base
    # ranked_hvgs is deterministic and its first C.N_HVGS entries are `base`, so the extras
    # are literally the next ones down the same list.
    assert extended[: C.N_HVGS] == base, "the extended ranking is not an extension of the panel"

    added: list[str] = []
    for gene in extended[C.N_HVGS:]:
        if gene.startswith(MT_PREFIX):
            continue          # a mitochondrial gene just below the cut, dropped for the same reason
        added.append(gene)
        if len(added) == len(mito):
            break
    if len(added) < len(mito):
        raise SystemExit(
            f"only {len(added)} non-mitochondrial genes available below rank {target}"
        )

    print(f"\nReplacements, ranks {C.N_HVGS + 1}-{target} of the same ranking:", flush=True)
    for rank, gene in enumerate(added, start=C.N_HVGS + 1):
        nb = int(adata.var.loc[gene, "highly_variable_nbatches"])
        print(f"  {gene:<12} highly variable in {nb} cohorts, "
              f"dispersion_norm {adata.var.loc[gene, 'dispersions_norm']:.3f}", flush=True)

    panel = kept + added
    assert len(panel) == C.N_HVGS, f"panel is {len(panel)} genes, expected {C.N_HVGS}"
    assert len(set(panel)) == C.N_HVGS, "duplicate genes in the panel"
    assert not any(g.startswith(MT_PREFIX) for g in panel), "an MT- gene survived"

    del adata

    # ---------------------------------------------------------------- the DRVI input
    # Built from the reduced object rather than from _norm_cc, so the categorical palettes in
    # .uns (which the 05_3 notebook reuses to keep a label the same colour across figures)
    # come along. Everything derived from the OLD panel is then either recomputed (PCA,
    # neighbours, UMAP) or dropped (the per-gene HVG statistics).
    print(f"\nReading {reduced_path}", flush=True)
    adata = sc.read_h5ad(reduced_path)
    print(f"Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes", flush=True)

    order = pd.Index(panel)
    unknown = order.difference(adata.var_names)
    assert unknown.empty, f"genes absent from {reduced_path.name}: {list(unknown)[:5]}"

    hvg = adata[:, order].copy()
    del adata
    assert hvg.n_vars == C.N_HVGS, f"got {hvg.n_vars} genes after subsetting"

    # The old selection's bookkeeping describes a panel this object no longer is.
    for col in ("highly_variable", "highly_variable_nbatches", "highly_variable_intersection",
                "means", "dispersions", "dispersions_norm", "highly_variable_rank"):
        if col in hvg.var:
            del hvg.var[col]
    hvg.var["highly_variable"] = True
    hvg.uns["hvg_panel"] = "nomt"
    hvg.uns["hvg_panel_dropped"] = mito
    hvg.uns["hvg_panel_added"] = added

    # Same guards reduce_data_tum.py puts on the file it hands to DRVI: the model reads
    # .layers['counts'] and ignores .X entirely, so a swap upstream has to fail here rather
    # than after hours of training.
    counts = hvg.layers["counts"]
    counts_values = counts.data if hasattr(counts, "data") else np.asarray(counts).ravel()
    assert np.allclose(counts_values, np.floor(counts_values)), (
        "layers['counts'] contains non-integer values: it is no longer raw counts."
    )
    assert counts_values.min() >= 0, "layers['counts'] contains negative values"
    x_values = hvg.X.data if hasattr(hvg.X, "data") else np.asarray(hvg.X).ravel()
    assert not np.allclose(x_values, np.floor(x_values)), (
        ".X looks like integer counts, but it should be scran log-normalized expression."
    )

    # PCA / neighbours / UMAP on the new panel, with reduce_data's parameters. The arrays
    # inherited from _reduced.h5ad describe the old gene set and are replaced, not kept.
    print("\nRecomputing PCA / neighbours / UMAP on the new panel...", flush=True)
    for key in ("X_pca", "X_umap"):
        hvg.obsm.pop(key, None)
    # mask_var=None, not use_highly_variable: the object is already the panel and nothing
    # else, so every gene is in - and the flag's meaning here is "this is the panel", not a
    # selection PCA still has to make.
    sc.pp.pca(hvg, n_comps=50, svd_solver="arpack", mask_var=None)
    sc.pp.neighbors(hvg, use_rep="X_pca")
    sc.tl.umap(hvg)

    hvg.write_h5ad(out_h5ad, compression="gzip")
    print(f"\nWrote {out_h5ad} "
          f"({hvg.n_obs:,} x {hvg.n_vars:,}, "
          f"{os.path.getsize(out_h5ad) / 1024 ** 3:.2f} GB on disk)", flush=True)

    # The list last: it is the cheap file, and writing it only after the h5ad means a crash
    # never leaves a list without the object it describes.
    pd.Series(hvg.var_names).to_csv(out_csv, index=False, header=False)
    print(f"Wrote {len(hvg.var_names)} genes to {out_csv}", flush=True)

    print(f"\n05_3 reads this with:  CELL_SET={C.cell_set()} HVG_SET=nomt "
          f"python3 run_drvi_tum.py --n-latent 64", flush=True)


if __name__ == "__main__":
    sys.exit(main())
