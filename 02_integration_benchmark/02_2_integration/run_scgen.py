#!/usr/bin/env python
"""
02_2 integration: scGen, in its own environment.

scGen has a separate script rather than a branch in run_integration.py because it
runs in the `scgen-py` environment (scvi-tools 0.19, torch 2.0), which cannot
coexist with the main stack. That environment does not have scib, so this calls
SCGEN directly instead of scib.integration.scgen, reproducing the same call
sequence scib uses (and the one validated in environments/check_scgen_env.py):
setup_anndata -> SCGEN -> train -> batch_removal.

Output type: full. scGen returns a corrected expression matrix in .X; the batch
correction is a full reconstructed feature space, not an embedding.

The input .h5ad was written by 02_1 with h5ad_compat, so anndata 0.10 (this
environment) can read it. Cell order is preserved: batch_removal returns the
cells in input order, and it is checked here regardless.

Usage:
  conda activate scgen-py
  export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
  python run_scgen.py -i $DATA_DIR/shiao_hvg_2k.h5ad \
                      -o $DATA_DIR/02_integration/scgen_unscaled.h5ad
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import anndata as ad
from scgen import SCGEN


def parse_args():
    p = argparse.ArgumentParser(description="Run scGen integration")
    p.add_argument("-i", "--input", required=True, help="prepared .h5ad input")
    p.add_argument("-o", "--output", required=True, help="integrated .h5ad output")
    p.add_argument("-b", "--batch-key", default="cohort")
    p.add_argument("-l", "--label-key", default="cell_type")
    p.add_argument("--epochs", type=int, default=100,
                   help="training epochs [default %(default)s, scib's value]")
    return p.parse_args()


def main():
    args = parse_args()

    print(f"[read] {args.input}", flush=True)
    adata = ad.read_h5ad(args.input)
    print(f"[read] {adata.n_obs:,} cells x {adata.n_vars:,} genes", flush=True)

    for key in (args.batch_key, args.label_key):
        assert key in adata.obs, f"missing obs column {key!r}"

    obs_names_before = adata.obs_names.to_numpy().copy()

    # Same parametrisation as scib.integration.scgen (from the scGen tutorial).
    print(f"[run] scGen (batch_key={args.batch_key}, labels_key={args.label_key}, "
          f"epochs={args.epochs})", flush=True)
    SCGEN.setup_anndata(adata, batch_key=args.batch_key, labels_key=args.label_key)
    model = SCGEN(adata)
    model.train(
        max_epochs=args.epochs,
        batch_size=32,
        early_stopping=True,
        early_stopping_patience=25,
    )
    corrected = model.batch_removal()

    # batch_removal returns a new object; make sure it matches cell for cell.
    assert corrected.n_obs == adata.n_obs, "scGen changed the cell count"
    assert np.array_equal(corrected.obs_names.to_numpy(), obs_names_before), (
        "scGen changed the cell order; the metrics would be computed on "
        "mismatched cells"
    )
    assert np.isfinite(corrected.X).all(), "scGen produced non-finite values"

    # Keep the output aligned with the Python dispatcher: full output in .X, no
    # raw counts (the metrics do not read them). scGen leaves its own latent space
    # in obsm['latent']; drop it to avoid it being mistaken for an X_emb.
    if "counts" in corrected.layers:
        del corrected.layers["counts"]
    # scGen leaves its own latent space in obsm; the key is 'corrected_latent' in
    # this scvi-tools version (older ones used 'latent'). Drop whichever is present
    # so it cannot be mistaken for an X_emb.
    for latent_key in ("latent", "corrected_latent"):
        corrected.obsm.pop(latent_key, None)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    print(f"[write] {args.output}", flush=True)
    corrected.write_h5ad(args.output, compression="gzip")

    size_gb = os.path.getsize(args.output) / 1024 ** 3
    print(f"[done] {corrected.n_obs:,} x {corrected.n_vars:,}, {size_gb:.2f} GB", flush=True)


if __name__ == "__main__":
    main()
