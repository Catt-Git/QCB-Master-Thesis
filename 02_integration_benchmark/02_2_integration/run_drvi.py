#!/usr/bin/env python
"""
02_2 integration: DRVI, headless (the compute half of shiao_drvi_128.ipynb).

DRVI is the one method of the grid that run_all.sh leaves out: the latent size is
chosen by eye, from how many dimensions vanish, so the run lives in a notebook.
That choice is a decision, not a computation - and the computation underneath is
the longest of the phase. This script is that computation without a kernel, so a
chosen latent size can be trained on the cluster (submit_drvi.slurm) or in a
terminal, and the notebook re-opened afterwards with OVERWRITE = False to look at
the result instead of recomputing it.

Nothing is parameterized differently from the notebook: same input, same
architecture (encoder [256, 128] / decoder [128, 256], dispersion='gene-batch'),
same batch_key, same seed, the same 400-epoch cap with early stopping after 50
epochs without improvement, and the same four outputs under DATA_DIR, all named
from the run id `drvi_unscaled_<n_latent>`:

    02_drvi/model_drvi_<N>.pt        the trained model, one flat file per run
    02_drvi/embed_<run_id>.h5ad      latent space + dimension stats + OOD/IND scores
    02_integration/<run_id>.h5ad     the benchmark output, obsm['X_emb'] (what 02_4 scores)
    02_embeddings/<run_id>.npy       the embedding alone, for 02_3

The last two mirror what run_integration.py writes for the other embed methods
(harmony, scVI, scANVI), so 02_3 and 02_4 read them as they are.

The notebook's figures are redrawn here too, into the same figures/<run_id>/
folder, so a cluster run leaves nothing to redo locally except looking at them. A
failure while plotting is reported but does not fail the run: by then every
artifact above is already on disk, and the training behind it is hours long.

Resuming is the default, as in the notebook: a step whose output already exists
is reported as [have] and reused, so a crash after training does not cost the
training. --overwrite recomputes everything instead.

Usage:
  export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
  python run_drvi.py                   # n_latent 128, the run in the grid
  python run_drvi.py --n-latent 64     # another size, side by side with the first
  python run_drvi.py --overwrite       # retrain and rewrite everything

On the cluster the same script is submitted by submit_drvi.slurm, which trains on
CPU (the cluster has no GPU) and takes its arguments unchanged.
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless: figures are written to disk, never shown
import matplotlib.pyplot as plt  # noqa: E402

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import scanpy as sc  # noqa: E402
import scvi  # noqa: E402
import drvi  # noqa: E402
from drvi.model import DRVI  # noqa: E402

BATCH_KEY = "cohort"       # the batch to correct, as everywhere else in the benchmark
LABEL_KEY = "cell_type"

# The keys of the notebook's UMAP panels, in the same order and with the same
# names, so a cluster run and a local run produce the same file names.
UMAP_QC_KEYS = {
    "cell_type": "cell_type",
    "fraction": "fraction",
    "cohort": "cohort",
    "treatment": "treatment",
    "response": "response",
    "phase": "phase",
    "n_genes_by_counts": "n_genes_by_counts",
    "total_counts": "total_counts",
    "mito": "pct_counts_mt",
    "ribo": "pct_counts_ribo",
    "size_factors": "size_factors",
}
UMAP_COMBINED_KEYS = ["cell_type", "fraction", "cohort", "treatment", "response", "phase"]
UMAP_SEED = 0


def parse_args():
    p = argparse.ArgumentParser(description="Run DRVI on the prepared Shiao object")
    p.add_argument("-n", "--n-latent", type=int, default=128,
                   help="latent dimensions; the run id follows it [default: 128]")
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--epochs", type=int, default=400,
                   help="max epochs, the maximum suggested in the literature; "
                        "early stopping usually fires first [default: 400]")
    p.add_argument("--data-dir", default=os.environ.get("DATA_DIR"),
                   help="directory holding the prepared inputs [default: $DATA_DIR]")
    p.add_argument("--fig-dir", default=None,
                   help="where the figures go [default: the repo's "
                        "02_integration_benchmark/figures/<run_id>]")
    p.add_argument("--overwrite", action="store_true",
                   help="retrain and rewrite everything instead of reusing what is on disk")
    return p.parse_args()


def save_model(model, path):
    """Save a DRVI model as the single file `path`.

    model.save() writes `<dir>/model.pt` and gives no say over the file name, so it
    writes into a scratch directory *inside* the destination folder (same
    filesystem, so the move below is a rename and not a copy) and the one file it
    produced is then renamed.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=path.parent) as scratch:
        model.save(scratch, overwrite=True)
        Path(scratch, "model.pt").replace(path)
    return path


def load_model(path, adata):
    """Load a DRVI model from the single file `path`; the counterpart of save_model.

    DRVI.load() insists on a directory containing `model.pt`, so it gets a scratch
    one holding a symlink to the real file: nothing is copied, and the scratch is
    gone by the time this returns.
    """
    path = Path(path).resolve()
    with tempfile.TemporaryDirectory(dir=path.parent) as scratch:
        os.symlink(path, Path(scratch, "model.pt"))
        return DRVI.load(scratch, adata)


def savefig(fig_dir, run_id, name, fig=None, dpi=300):
    """Save a matplotlib figure (the current one by default) into fig_dir.

    For the DRVI and seaborn plots, which draw on the pyplot state instead of
    going through sc.pl and so ignore sc.settings.figdir. The run id is appended
    to the name, so figures from different latent sizes stay distinguishable even
    once they are pulled out of their folder. Headless, the figure is closed
    rather than shown: hundreds of them are drawn in one process.
    """
    fig = plt.gcf() if fig is None else fig
    path = Path(fig_dir) / f"{name}_{run_id}.png"
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {path}", flush=True)
    return path


def draw_figures(model, adata, embed, fig_dir, run_id):
    """The plotting section of the notebook, in order and with the same names.

    The five standard QC panels that compare DRVI against the other methods
    (cohort / cell_type, integrated vs unintegrated) are deliberately NOT drawn
    here: 02_3_plot_method_umap draws them for every method with the same palette,
    the same seeded point order and the same UMAP parameters, and that uniformity
    is the whole point of the comparison.
    """
    sc.settings.figdir = fig_dir

    # Same seeded permutation as the notebook, so overplotting hides the same
    # cells in both.
    order = np.random.default_rng(UMAP_SEED).permutation(embed.n_obs)
    embed_plot = ad.AnnData(
        obs=embed.obs.iloc[order].copy(),
        obsm={"X_umap": embed.obsm["X_umap"][order]},
        # Palettes carried over from the input object, so a cell type keeps the same
        # colour here, in the 01_6 unintegrated UMAPs and in 03_1.
        uns={k: v for k, v in adata.uns.items() if k.endswith("_colors")},
    )

    for label, col in UMAP_QC_KEYS.items():
        if col not in embed_plot.obs:
            print(f"[skip] umap {label}: no obs column {col!r}", flush=True)
            continue
        sc.pl.umap(embed_plot, color=col, show=False, save=f"_{label}_{run_id}.png")
        plt.close("all")

    combined = [k for k in UMAP_COMBINED_KEYS if k in embed_plot.obs]
    with plt.rc_context({"figure.figsize": (7, 7)}):
        sc.pl.umap(embed_plot, color=combined, ncols=2, wspace=0.8,
                   hspace=0.25, show=False, save=f"_combined_{run_id}.png")
        plt.close("all")

    # One panel per cell type: with ~50 CellTypist labels the single `cell_type`
    # UMAP is unreadable, and where each population sits is the point.
    cell_types = [c for c in embed_plot.obs["cell_type"].cat.categories
                  if (embed_plot.obs["cell_type"] == c).any()]
    ncols = 4
    nrows = int(np.ceil(len(cell_types) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.4))
    axes = np.atleast_1d(axes).flatten()
    for ax, ct in zip(axes, cell_types):
        sc.pl.umap(embed_plot, color="cell_type", groups=[ct], title=ct, ax=ax,
                   show=False, legend_loc="none", na_in_legend=False, size=4)
        ax.set_xlabel("")
        ax.set_ylabel("")
    for ax in axes[len(cell_types):]:
        ax.set_visible(False)
    plt.tight_layout()
    savefig(fig_dir, run_id, "umap_per_cell_type", fig=fig)
    del embed_plot

    # Latent dimensions: per-dimension statistics, with and without the vanished
    # ones. The second is the plot behind the choice of n_latent.
    drvi.utils.pl.plot_latent_dimension_stats(embed, ncols=2, show=False)
    savefig(fig_dir, run_id, "latent_dimension_stats")
    drvi.utils.pl.plot_latent_dimension_stats(embed, ncols=2, remove_vanished=True, show=False)
    savefig(fig_dir, run_id, "latent_dimension_stats_rmVanished")

    # Every non-vanished dimension on the UMAP (remove_vanished defaults to True).
    drvi.utils.pl.plot_latent_dims_in_umap(embed, show=False)
    savefig(fig_dir, run_id, "latent_dims_in_umap")

    # Which dimension responds to which cell type. Sorted by dimension first, then
    # grouped by label: the second ordering makes it obvious when one label owns a
    # whole dimension.
    drvi.utils.pl.plot_latent_dims_in_heatmap(embed, "cell_type", title_col="title", show=False)
    savefig(fig_dir, run_id, "latent_dims_in_heatmap_cell_type")
    drvi.utils.pl.plot_latent_dims_in_heatmap(embed, "cell_type", title_col="title",
                                              sort_by_categorical=True, show=False)
    savefig(fig_dir, run_id, "latent_dims_in_heatmap_cell_type_sorted")

    # The same heatmap against the batch, then treatment, then phase.
    for key in ("cohort", "treatment", "phase"):
        drvi.utils.pl.plot_latent_dims_in_heatmap(embed, key, title_col="title", show=False)
        savefig(fig_dir, run_id, f"latent_dims_in_heatmap_{key}")

    # Interpretability. OOD comes from the decoder reconstructions (fast, favours
    # the genes *specific* to a dimension); OOD_min/max_possible are its two
    # halves; IND averages each factor's effect over all cells (broader).
    model.plot_interpretability_scores(embed, adata, show=False)
    savefig(fig_dir, run_id, "ood_interpretability_scores")
    model.plot_interpretability_scores(embed, adata, key="OOD_max_possible", show=False)
    savefig(fig_dir, run_id, "ood_max_interpretability_scores")
    model.plot_interpretability_scores(embed, adata, key="OOD_min_possible", show=False)
    savefig(fig_dir, run_id, "ood_min_interpretability_scores")
    model.plot_interpretability_scores(embed, adata, key="IND_linear_weighted_mean", show=False)
    savefig(fig_dir, run_id, "ind_linear_weighted_mean")


def main():
    args = parse_args()

    if not args.data_dir:
        sys.exit("set DATA_DIR (or pass --data-dir) to the directory holding the prepared inputs")

    # The run id follows the grid's `drvi_unscaled` with the latent size appended,
    # so the 32 / 64 / 128 runs sit side by side instead of overwriting each other.
    run_id = f"drvi_unscaled_{args.n_latent}"

    script_dir = Path(__file__).resolve().parent
    phase_dir = script_dir.parent                     # 02_integration_benchmark/
    data_dir = Path(args.data_dir)

    # Input: the same prepared object every other method of the grid integrates.
    input_h5ad = data_dir / "shiao_hvg_2k.h5ad"

    # Outputs. All of them live under DATA_DIR, never in the repo: the latent .h5ad
    # alone is ~730 MB, well past GitHub's per-file limit, and /datasets/* is
    # gitignored.
    model_dir = data_dir / "02_drvi"                              # models and latent spaces, flat
    model_path = model_dir / f"model_drvi_{args.n_latent}.pt"     # one file (see save_model)
    latent_h5ad = model_dir / f"embed_{run_id}.h5ad"              # latent space + interpretability
    integrated_h5ad = data_dir / "02_integration" / f"{run_id}.h5ad"   # scored object, obsm['X_emb']
    emb_npy = data_dir / "02_embeddings" / f"{run_id}.npy"             # small durable copy

    fig_dir = Path(args.fig_dir) if args.fig_dir else phase_dir / "figures" / run_id

    for d in (model_dir, integrated_h5ad.parent, emb_npy.parent, fig_dir):
        d.mkdir(parents=True, exist_ok=True)

    print(f"run id      {run_id}", flush=True)
    print(f"input       {input_h5ad}")
    print(f"model       {model_path}")
    print(f"latent      {latent_h5ad}")
    print(f"integrated  {integrated_h5ad}")
    print(f"embedding   {emb_npy}")
    print(f"figures     {fig_dir}")
    print(f"scvi-tools {scvi.__version__}, DRVI {drvi.__version__}", flush=True)

    if not input_h5ad.exists():
        sys.exit(f"missing {input_h5ad}: run 02_1_prepare first")

    sc.settings.set_figure_params(dpi=300, facecolor="white")
    sc.settings.figdir = fig_dir

    adata = sc.read_h5ad(input_h5ad)
    print(f"{adata.n_obs:,} cells x {adata.n_vars:,} HVGs", flush=True)
    assert "counts" in adata.layers, "raw counts layer missing, DRVI has nothing to train on"

    # Setup the anndata for scvi-tools
    DRVI.setup_anndata(
        adata,
        layer="counts",       # raw unnormalised data (used by default)
        batch_key=BATCH_KEY,  # the batch to fix
    )

    # Setting seed
    scvi.settings.seed = args.seed

    # Setup the model
    model = DRVI(
        adata,
        n_latent=args.n_latent,
        encoder_dims=[256, 128],
        decoder_dims=[128, 256],
        dispersion="gene-batch",   # really important parameter for TNBC
    )

    # Train only if the model is not on disk yet, or if --overwrite asks for a
    # retrain. The training uses the GPU if one is available, otherwise the CPU.
    if args.overwrite or not model_path.exists():
        print(f">>> training {run_id}: up to {args.epochs} epochs, "
              f"early stopping after 50 without improvement", flush=True)
        model.train(
            max_epochs=args.epochs,
            early_stopping=True,             # activate early stopping
            early_stopping_patience=50,      # stop after 50 epochs without improvement
        )
        save_model(model, model_path)
        print(f"[model] saved to {model_path}", flush=True)
    else:
        print(f"[have] {model_path} already exists, not retraining (--overwrite to force)",
              flush=True)

    # Load the model back, so a resumed run works from exactly the same state a
    # fresh one does.
    model = load_model(model_path, adata)

    # The DRVI-specific view of the run: the latent space as .X, one var per latent
    # dimension carrying its stats and interpretability scores.
    if args.overwrite or not latent_h5ad.exists():
        embed = ad.AnnData(model.get_latent_representation(), obs=adata.obs)

        # We set latent dimension statistics
        print("Setting latent dimension stats ...", flush=True)
        model.set_latent_dimension_stats(embed, vanished_threshold=0.5)

        # We immediately calculate the interpretability gene scores with different approaches
        print("Calculating gene scores per factor ...", flush=True)
        # out-of-distribution (OOD) approach uses decoder reconstructions to calculate gene scores (faster)
        model.calculate_interpretability_scores(embed, "OOD")
        # within-distribution (IND) approach iterates over all cells and calculates gene scores
        model.calculate_interpretability_scores(embed, "IND")

        print("Dimension reduction ...", flush=True)
        sc.pp.neighbors(embed, n_neighbors=15, use_rep="X", n_pcs=embed.X.shape[1])
        sc.tl.umap(embed)
        sc.pp.pca(embed)

        print(f"[write] {latent_h5ad}", flush=True)
        embed.write_h5ad(latent_h5ad)
    else:
        print(f"[have] {latent_h5ad}, reading it back", flush=True)
        embed = sc.read_h5ad(latent_h5ad)

    n_vanished = int(embed.var["vanished"].sum())
    print(f"{n_vanished} vanished / {embed.n_vars} latent dimensions "
          f"({embed.n_vars - n_vanished} effectively used)", flush=True)

    # The benchmark output: the same shape run_integration.py writes for the other
    # embed methods, so 02_3 and 02_4 read it without knowing DRVI produced it.
    if args.overwrite or not integrated_h5ad.exists():
        # Cell order is load-bearing for the metrics, so the latent space is realigned
        # on adata's cells rather than assumed to be in the same order.
        assert set(embed.obs_names) == set(adata.obs_names), "embed and adata disagree on cells"

        integrated = adata.copy()
        integrated.obsm["X_emb"] = np.asarray(embed[adata.obs_names].X, dtype=np.float32)
        assert np.isfinite(integrated.obsm["X_emb"]).all()

        # The metrics never read raw counts; dropping them keeps the file light.
        if "counts" in integrated.layers:
            del integrated.layers["counts"]

        print(f"[write] {integrated_h5ad}", flush=True)
        integrated.write_h5ad(integrated_h5ad, compression="gzip")

        np.save(emb_npy, integrated.obsm["X_emb"])
        print(f"[write] {emb_npy}  {integrated.obsm['X_emb'].shape}", flush=True)

        del integrated
    else:
        print(f"[have] {integrated_h5ad} (--overwrite to rewrite)", flush=True)

    print(">>> figures", flush=True)
    try:
        draw_figures(model, adata, embed, fig_dir, run_id)
    except Exception:
        # Every artifact is already written and the training behind it is hours
        # long: report the failure, keep the exit status, and let the notebook
        # redraw whichever plot broke.
        print("[warn] the figures section failed; the artifacts above are complete",
              flush=True)
        traceback.print_exc()

    print(f"[ok] {run_id} -> {integrated_h5ad}", flush=True)


if __name__ == "__main__":
    main()
