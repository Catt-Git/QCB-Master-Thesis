#!/usr/bin/env python
"""
03_2 DRVI run: the headless half of drvi_nonimm.ipynb.

The notebook is where the latent size is chosen, by eye, from how many dimensions
vanish; the computation underneath it is deterministic and long, and this script
is that computation without a kernel. It exists so a chosen size can be trained
on the cluster (submit_drvi_nonimm.slurm) or in a terminal and the notebook
re-opened afterwards with OVERWRITE = False, which reads the model and the
embedding from disk instead of recomputing them.

Nothing is parameterized differently from the notebook: same input, same
architecture (encoder [256, 128] / decoder [128, 256], dispersion='gene-batch'),
`batch_key='cohort'`, `SEED=123`, the same 400-epoch cap with early stopping
after 50 epochs without improvement, and the same three outputs, named from the
run id `drvi_nonimm_<n_latent>`:

    03_nonimm/model_<run_id>.pt              the trained model, one flat file per run
    03_nonimm/embed_<run_id>.h5ad            latent space + dimension stats + OOD/IND scores
    03_nonimm/shiao_nonimm_<run_id>.h5ad     the 03_1 object (all genes) + obsm['X_drvi']

Those three are everything this phase produces: 02_3 and 02_4 do not run on this
compartment, so nothing is written for the scib benchmark. 03_3 reads the
embedding directly (it needs neither the model nor a GPU); the third file is the
compartment itself in the DRVI space, for whatever downstream step needs the
genes alongside the latent coordinates. It is the only output that needs
`shiao_nonimm.h5ad` to be around: if that object is missing the run says so up
front, writes the other two and skips it.

As in the notebook, this is a full re-run on the compartment and inherits nothing
from the 02_2 whole-dataset run: `drvi_nonimm_*` and `drvi_unscaled_*` share the
recipe and nothing else.

The notebook's figures are redrawn here too, into the same figures/03_2_<run_id>/
folder, so a cluster run leaves nothing to redo locally except looking at them. A
failure while plotting is reported but does not fail the run: by then every
artifact above is on disk and the training behind it is hours long.

Resuming is the default, as everywhere in this phase: a step whose output already
exists is reported as [have] and reused, so a crash after training does not cost
the training. --overwrite recomputes everything instead.

Usage:
  export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
  python run_drvi_nonimm.py                 # n_latent 64, the run of this phase
  python run_drvi_nonimm.py --n-latent 32   # another size, side by side with it
  python run_drvi_nonimm.py --overwrite     # retrain and rewrite everything

On the cluster the same script is submitted by submit_drvi_nonimm.slurm, which
trains on CPU (the cluster has no GPU) and takes its arguments unchanged.
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

BATCH_KEY = "cohort"       # the batch to correct, as in phase 02 and in 03_1
LABEL_KEY = "cell_type"    # CellTypist, 18 labels observed in this compartment

# The keys of the notebook's UMAP panels, in the same order and with the same
# names. `fraction` is not among them: it is constant (`non_imm`) after the 03_1
# subset and carries no information, and `dataset_origin` (the technical CD45
# sort) is left out as everywhere in this phase.
UMAP_QC_KEYS = {
    "cell_type": "cell_type",
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
UMAP_COMBINED_KEYS = ["cell_type", "cohort", "treatment", "response", "phase"]
UMAP_SEED = 0


def parse_args():
    p = argparse.ArgumentParser(description="Run DRVI on the non-immune compartment")
    p.add_argument("-n", "--n-latent", type=int, default=64,
                   help="latent dimensions; the run id follows it [default: 64]")
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--epochs", type=int, default=400,
                   help="max epochs, the maximum suggested in the literature; "
                        "early stopping usually fires first [default: 400]")
    p.add_argument("--data-dir", default=os.environ.get("DATA_DIR"),
                   help="directory holding the datasets [default: $DATA_DIR]")
    p.add_argument("--fig-dir", default=None,
                   help="where the figures go [default: the repo's "
                        "03_drvi_non_immune/figures/03_2_<run_id>]")
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

    Each UMAP panel has its direct counterpart in
    figures/03_1_visualization/umap_*_nonimm.png: the same cells in the
    *unintegrated* PCA space, with the same keys and the same palettes, so only
    the space changes between the two.
    """
    sc.settings.figdir = fig_dir

    # Same seeded permutation as the notebook, so overplotting hides the same
    # cells in both.
    order = np.random.default_rng(UMAP_SEED).permutation(embed.n_obs)
    embed_plot = ad.AnnData(
        obs=embed.obs.iloc[order].copy(),
        obsm={"X_umap": embed.obsm["X_umap"][order]},
        # Palettes carried over from the input object, so a cell type keeps the same
        # colour here and in the 03_1 UMAPs of the same cells.
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

    # One panel per cell type: the three largest labels (Fibro-matrix, Lumsec-prol,
    # LummHR-major) dominate the single `cell_type` panel, and where each
    # population sits is the point.
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

    # The same heatmap against the batch, then treatment, response and phase.
    for key in ("cohort", "treatment", "response", "phase"):
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
        sys.exit("set DATA_DIR (or pass --data-dir) to the directory holding the datasets")

    # `nonimm` in the run id, so nothing here can be confused with the whole-dataset
    # `drvi_unscaled_*` runs of 02_2, and the latent size keeps the sizes side by side.
    run_id = f"drvi_nonimm_{args.n_latent}"

    script_dir = Path(__file__).resolve().parent
    phase_dir = script_dir.parent                     # 03_drvi_non_immune/
    data_dir = Path(args.data_dir)
    nonimm_dir = data_dir / "03_nonimm"               # every heavy file of this phase lives here

    # Input: the 2,000-gene object written by 03_1's reduce_data_nonimm.py. HVGs were
    # selected on the non-immune cells only, .layers['counts'] are the raw counts DRVI
    # trains on.
    input_h5ad = nonimm_dir / "shiao_nonimm_hvg_2k.h5ad"
    # The definitive object of 03_1 (all genes, scran log-norm, leiden): read once
    # at the end, only to carry the latent space over for 03_3.
    full_h5ad = nonimm_dir / "shiao_nonimm.h5ad"

    # Outputs, all under DATA_DIR and never in the repo (the embedding alone is a
    # few hundred MB, past GitHub's per-file limit, and /datasets/* is gitignored).
    model_path = nonimm_dir / f"model_{run_id}.pt"        # one file (see save_model)
    latent_h5ad = nonimm_dir / f"embed_{run_id}.h5ad"     # latent space + interpretability
    downstream_h5ad = nonimm_dir / f"shiao_nonimm_{run_id}.h5ad"   # 03_1 object + obsm['X_drvi']

    fig_dir = Path(args.fig_dir) if args.fig_dir else phase_dir / "figures" / f"03_2_{run_id}"

    nonimm_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"run id      {run_id}", flush=True)
    print(f"input       {input_h5ad}")
    print(f"model       {model_path}")
    print(f"latent      {latent_h5ad}")
    print(f"downstream  {downstream_h5ad}")
    print(f"figures     {fig_dir}")
    print(f"scvi-tools {scvi.__version__}, DRVI {drvi.__version__}", flush=True)

    if not input_h5ad.exists():
        sys.exit(f"missing {input_h5ad}: run 03_1_subsetting/subsetting_all.sh first")

    # Checked here rather than at the end, where it is used: the downstream object
    # is written after the training, and a missing input discovered there would
    # only surface hours in. Not fatal - the model and the embedding do not need
    # it - but worth saying before anything starts, especially on the cluster,
    # where this file has to be copied up as well.
    if not full_h5ad.exists():
        print(f"[warn] {full_h5ad} is missing: the run will produce the model and the "
              f"embedding,\n       but not {downstream_h5ad.name}", flush=True)

    sc.settings.set_figure_params(dpi=300, facecolor="white")
    sc.settings.figdir = fig_dir

    adata = sc.read_h5ad(input_h5ad)
    print(f"{adata.n_obs:,} cells x {adata.n_vars:,} HVGs", flush=True)

    # Sanity checks on what 03_1 produced: one compartment only, the counts layer present.
    assert adata.obs["fraction"].astype(str).nunique() == 1, "not a single-compartment object"
    assert "counts" in adata.layers, "raw counts layer missing, DRVI has nothing to train on"
    print("fraction    :", adata.obs["fraction"].astype(str).unique().tolist())
    print(f"{adata.obs[BATCH_KEY].nunique()} cohorts, {adata.obs[LABEL_KEY].nunique()} cell types",
          flush=True)

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
    # dimension carrying its stats and interpretability scores. This is the file
    # 03_3 reads.
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

    # The number the latent size is chosen on: how many of the n_latent dimensions
    # DRVI did not need. This is what to read in the log before deciding whether to
    # re-run at another size.
    n_vanished = int(embed.var["vanished"].sum())
    print(f"{n_vanished} vanished / {embed.n_vars} latent dimensions "
          f"({embed.n_vars - n_vanished} effectively used)", flush=True)

    # The object for anything downstream that needs genes AND the latent space:
    # the definitive 03_1 object (all genes, scran log-norm, leiden) with the
    # latent space added as obsm['X_drvi']. The embedding above is the model's own
    # view of the run - one var per latent dimension, its stats and its
    # interpretability scores - and by construction carries no genes; this is the
    # other half, the compartment itself in the DRVI space. 03_3 as it stands
    # works from the embedding alone and does not open this file.
    if not full_h5ad.exists():
        print(f"[warn] {full_h5ad} is missing, so {downstream_h5ad.name} was not written.\n"
              f"       Run 03_1_subsetting/subsetting_all.sh (or copy that object here) and\n"
              f"       re-run: the model and the embedding above are already done, so this\n"
              f"       second pass only reads them back.", flush=True)
    elif args.overwrite or not downstream_h5ad.exists():
        print(f"[read] {full_h5ad}", flush=True)
        full = sc.read_h5ad(full_h5ad)

        # Cell order is not assumed: the latent space is realigned on the 03_1
        # object's cells by name, so a reordering anywhere upstream cannot pair a
        # cell with another cell's coordinates.
        assert set(embed.obs_names) == set(full.obs_names), \
            "the embedding and the 03_1 object disagree on cells"
        full.obsm["X_drvi"] = np.asarray(embed[full.obs_names].X, dtype=np.float32)
        assert np.isfinite(full.obsm["X_drvi"]).all()

        print(f"Writing {downstream_h5ad.name} (gzip, a few minutes) ...", flush=True)
        full.write_h5ad(downstream_h5ad, compression="gzip")
        print(f"[write] {downstream_h5ad} "
              f"({full.n_obs:,} x {full.n_vars:,}, obsm['X_drvi'] "
              f"{full.obsm['X_drvi'].shape}, "
              f"{downstream_h5ad.stat().st_size / 1024 ** 3:.2f} GB on disk)", flush=True)
        del full
    else:
        print(f"[have] {downstream_h5ad} (--overwrite to rewrite)", flush=True)

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

    print(f"[ok] {run_id} -> {latent_h5ad}", flush=True)
    print("next: 03_3_enrichment/enrichment_nonimm.ipynb on that embedding", flush=True)


if __name__ == "__main__":
    main()
