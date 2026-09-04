#!/usr/bin/env python
"""
05_3 DRVI run: the headless half of drvi_tum.ipynb.

The notebook is where the latent size is confirmed, by eye, from how many dimensions
vanish; the computation underneath it is deterministic and long, and this script is
that computation without a kernel. It exists so a chosen size can be trained on the
cluster (submit_drvi_tum.slurm) or in a terminal and the notebook re-opened afterwards
with OVERWRITE = False, which reads the model and the embedding from disk instead of
recomputing them.

Same recipe as 04_2/run_drvi_epi.py - same architecture (encoder [256, 128] / decoder
[128, 256], dispersion='gene-batch'), `batch_key='cohort'`, `SEED=123`, 400 full epochs
with no early stopping - on the object 05_2 wrote, with the three differences below.

**1. The input is the malignant subset.** `$DATA_DIR/05_tum/shiao_tum_hvg_2k.h5ad`:
36,192 cells x 2,000 HVGs, 19 cohorts, HVGs re-selected inside the tumour by
`05_2/reduce_data_tum.py`. Nothing is inherited from 04 - not the model, not the HVGs,
not the latent size - and the run ids do not collide (`drvi_tum_*` against 04's
`drvi_epi_*`), so both phases keep their outputs.

**2. The default latent size is 32, not 64.** 04 settled at 64 after 32 left *none* of
its dimensions vanished on 74,441 cells and 10 labels. This subset is half of it, 36,192
cells of one compartment whose annotation is a single constant value, so 32 is where the
scan starts here rather than where it ended there. The number to read in the log is the
same one: `n vanished / n_latent`. Several vanished dimensions mean 32 was already
generous; none means the space is too tight and the run is worth repeating at 64
(`--n-latent 64`, side by side - the run id carries the size, so nothing is overwritten).

**3. `cell_type` is constant and cannot group anything.** After the CNV call every cell
in this object is labelled `malignant`, so every figure 04_2 draws *by cell type* - the
per-label UMAP panels, the latent-dimension heatmaps, the SMI screening of the notebook -
would be drawn against one category. Two groupings replace it, in this order:

  `optscib_tum_leiden`  05_2's clustering, computed ON these cells. Attached from the
                        definitive object (only `obs` is read). It leads every figure that
                        takes a grouping. Not fully independent of the label below, and
                        the dependency is exact: the clusters come from the expression
                        graph, the RESOLUTION was picked by maximising NMI against
                        `cell_type_01_4` (`clustering_tum.py`, recorded in
                        `uns['optscib_tum_leiden_label_key']`) - the borrowed label chose
                        how many clusters there are, not which cell goes in which.
  `cell_type_01_4`      the pre-CNV CellTypist label, not constant inside the tumour (as
                        counted in the DRVI input itself: 18,617 Lumsec-prol, 11,771
                        Lumsec-basal, 3,893 LummHR-SCGB, 1,167 Lumsec-KIT, 594
                        LummHR-major and 148 across four more). A landmark, second.

The second one is obsolete as an *identity* - that is what phase 05 exists to establish -
and it is kept because on this subset it is the only non-constant CellTypist column left:
the post-CNV re-annotation of 05_1 ran on the non-malignant cells and leaves
`celltypist_predicted_cnv` at the constant `malignant` here, while the raw pre-voting
`celltypist_predicted` puts Fibro-matrix on 1,871 aneuploid cells and pericytes on 658. It
is the same substitution `05_2/clustering_tum.py` makes for its NMI target, with the same
caveat: a *borrowed grouping used to read the latent space*, never propagated as a
biological claim about these cells.

Neither grouping is clean, and the two fail differently: the CellTypist label comes from a
model with no malignant class; the leiden partition comes from the *unintegrated* PCA of
05_2, so it carries part of the cohort structure DRVI is meant to correct, and its
resolution was picked against that same CellTypist label. A dimension that answers to both
and not to `cohort` is the safe reading; one that answers to leiden alone is a state the
normal-breast vocabulary has no word for, which is what this phase is looking for.

None of this reaches the model. `DRVI.setup_anndata` takes a `labels_key` and it is left
None here: the counts layer and `batch_key='cohort'` are the whole input, so the latent
space is built without knowing any of these labels and none of them can bias it. They are
read back afterwards, as captions on dimensions that already exist. The caption that owes
nothing to any annotation is the third one and it is written by this script too: the OOD /
IND interpretability scores, which name a dimension by its genes. That is the route 05_6
and 05_7 take; the groupings here are the cross-check on it, never the evidence.

Two keys are plotted here that 04_2 has no equivalent for: `cnv_score` and `cnv_corr`
from 05_1. If the dominant dimensions of this space tracked CNV burden, the latent space
would be describing how aneuploid a cell is rather than what state it is in - these are
the two panels that show it.

Outputs, named from the run id `drvi_<compartment>_<n_latent>`:

    05_tum/model_<run_id>.pt              the trained model, one flat file per run
    05_tum/embed_<run_id>.h5ad            latent space + dimension stats + OOD/IND scores
    05_tum/shiao_tum_<run_id>.h5ad        the 05_2 object (all genes) + obsm['X_drvi']

Those three are everything this step produces: 02_3 and 02_4 do not run on this
compartment, so nothing is written for the scib benchmark. 05_4 reads the embedding
directly (it needs neither the model nor a GPU); the third file is the subset itself in
the DRVI space, for the downstream steps that need the genes alongside the latent
coordinates (05_5 cytotrace2, 05_6 cell-first, 05_9 cycle-confound). It is the only
output that needs `shiao_tum.h5ad` to be around: if that object is missing the run says
so up front, writes the other two and skips it.

`CELL_SET` is honoured as everywhere in this phase, through `05_2/cell_set.py`: unset (or
`tum`) trains on the malignant cells, `CELL_SET=epi` on the control set of every
epithelial cell under the post-CNV labels. The two write different prefixes and different
run ids (`drvi_tum_32` against `drvi_epicnv_32`), so both can live in `05_tum/`.

The notebook's figures are redrawn here too, into the same figures/05_3_<run_id>/ folder,
so a cluster run leaves nothing to redo locally except looking at them. A failure while
plotting is reported but does not fail the run: by then every artifact above is on disk
and the training behind it is hours long.

Resuming is the default, as everywhere in this phase: a step whose output already exists
is reported as [have] and reused, so a crash after training does not cost the training.
--overwrite recomputes everything instead.

Usage:
  export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
  python3 run_drvi_tum.py                  # n_latent 32, the run of this phase
  python3 run_drvi_tum.py --n-latent 64    # another size, side by side with it
  python3 run_drvi_tum.py --overwrite      # retrain and rewrite everything
  python3 run_drvi_tum.py --early-stopping # stop early instead (see --epochs)
  CELL_SET=epi python3 run_drvi_tum.py     # the same run on the control set

On the cluster the same script is submitted by submit_drvi_tum.slurm, which trains on CPU
(the cluster has no GPU) and takes its arguments unchanged.
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

# The CELL_SET -> prefix mapping and every threshold of this phase live in one module,
# next to the scripts that wrote the input. Imported rather than duplicated for the
# reason its own docstring gives: a prefix recomputed independently in each script is a
# silent-overwrite bug waiting to happen.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "05_2_subsetting"))
import cell_set as C  # noqa: E402

BATCH_KEY = C.BATCH_KEY              # 'cohort', the batch to correct, as in 02, 03, 04 and 05_2
LABEL_KEY = C.LABEL_KEY              # 'cell_type': the POST-CNV label, constant under CELL_SET=tum
GROUP_KEY = C.PRIOR_LABEL_KEY        # 'cell_type_01_4': the pre-CNV label, what actually groups
LEIDEN_KEY = "optscib_tum_leiden"    # 05_2's clustering, as named by clustering_tum.py

# The keys of the notebook's UMAP panels, in the same order and with the same names.
# `compartment`, `fraction` and `cnv_status` are not among them: after the 05_2 subset all
# three are constant ('tum', 'non_imm', 'malignant') and carry no information, and so is
# `cell_type` under CELL_SET=tum - it is listed because it is the interesting key of the
# `epi` control set, and dropped at draw time wherever it is constant. `dataset_origin`
# (the technical CD45 sort) is left out as everywhere in Part 2. `cnv_score` and
# `cnv_corr` are new against 04_2 and are the two panels that say whether this latent
# space is organised by aneuploidy.
UMAP_QC_KEYS = {
    "cell_type": LABEL_KEY,          # constant under `tum`, skipped; the point of `epi`
    "leiden": LEIDEN_KEY,            # the primary grouping; only there once attach_leiden ran
    "cell_type_01_4": GROUP_KEY,     # the secondary one, a landmark (see PANEL_KEYS)
    "cohort": BATCH_KEY,
    "treatment": "treatment",
    "response": "response",
    "phase": "phase",
    "cnv_score": "cnv_score",
    "cnv_corr": "cnv_corr",
    "n_genes_by_counts": "n_genes_by_counts",
    "total_counts": "total_counts",
    "mito": "pct_counts_mt",
    "ribo": "pct_counts_ribo",
    "size_factors": "size_factors",
}
UMAP_COMBINED_KEYS = [LABEL_KEY, LEIDEN_KEY, GROUP_KEY, BATCH_KEY, "treatment", "response",
                      "phase"]

# The categorical keys the latent dimensions are read against, in order; the first one also
# gets the extra grouped ordering. LEIDEN_KEY is only there once attach_leiden has run, and
# LABEL_KEY only says anything under CELL_SET=epi - `groupable` drops whatever is constant
# or missing rather than drawing a one-row heatmap.
#
# The order is the argument of this step. `optscib_tum_leiden` comes first because it is the
# only grouping of these cells computed ON these cells; `cell_type_01_4` follows as a
# landmark, and if the leiden column could not be attached it takes the first slot by
# default rather than by choice. Neither is clean: the CellTypist label comes from a
# normal-breast model with no malignant class, and the leiden partition comes from the
# UNINTEGRATED PCA of 05_2, so it carries some of the cohort structure DRVI is meant to
# correct. A dimension that answers to both, and not to `cohort`, is the safe reading.
HEATMAP_KEYS = [LEIDEN_KEY, GROUP_KEY, BATCH_KEY, "treatment", "response", "phase", LABEL_KEY]

# The groupings that get one UMAP panel per level, in order. Same argument as HEATMAP_KEYS.
PANEL_KEYS = [LEIDEN_KEY, GROUP_KEY]

UMAP_SEED = 0


def parse_args():
    p = argparse.ArgumentParser(description="Run DRVI on the malignant epithelial subset")
    p.add_argument("-n", "--n-latent", type=int, default=32,
                   help="latent dimensions; the run id follows it [default: 32]")
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--epochs", type=int, default=400,
                   help="epochs to train; scvi-tools warms the KL term up over exactly "
                        "this many by default [default: 400]")
    p.add_argument("--early-stopping", action="store_true",
                   help="stop after 50 epochs without improvement. Off by default: the KL "
                        "warmup lasts --epochs epochs, so a run that stops earlier never "
                        "trains at full KL weight, and it is that KL pressure that makes "
                        "the unused latent dimensions vanish. Pair it with a shorter warmup "
                        "if you turn it on")
    p.add_argument("--data-dir", default=os.environ.get("DATA_DIR"),
                   help="directory holding the datasets [default: $DATA_DIR]")
    p.add_argument("--fig-dir", default=None,
                   help="where the figures go [default: the repo's "
                        "05_drvi_tumoral_epi/figures/05_3_<run_id>]")
    p.add_argument("--overwrite", action="store_true",
                   help="retrain and rewrite everything instead of reusing what is on disk")
    return p.parse_args()


def save_model(model, path):
    """Save a DRVI model as the single file `path`.

    model.save() writes `<dir>/model.pt` and gives no say over the file name, so it
    writes into a scratch directory *inside* the destination folder (same filesystem, so
    the move below is a rename and not a copy) and the one file it produced is then
    renamed.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=path.parent) as scratch:
        model.save(scratch, overwrite=True)
        Path(scratch, "model.pt").replace(path)
    return path


def load_model(path, adata):
    """Load a DRVI model from the single file `path`; the counterpart of save_model.

    DRVI.load() insists on a directory containing `model.pt`, so it gets a scratch one
    holding a symlink to the real file: nothing is copied, and the scratch is gone by the
    time this returns.
    """
    path = Path(path).resolve()
    with tempfile.TemporaryDirectory(dir=path.parent) as scratch:
        os.symlink(path, Path(scratch, "model.pt"))
        return DRVI.load(scratch, adata)


def savefig(fig_dir, run_id, name, fig=None, dpi=300):
    """Save a matplotlib figure (the current one by default) into fig_dir.

    For the DRVI and seaborn plots, which draw on the pyplot state instead of going
    through sc.pl and so ignore sc.settings.figdir. The run id is appended to the name,
    so figures from different latent sizes stay distinguishable even once they are pulled
    out of their folder. Headless, the figure is closed rather than shown: hundreds of
    them are drawn in one process.
    """
    fig = plt.gcf() if fig is None else fig
    path = Path(fig_dir) / f"{name}_{run_id}.png"
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {path}", flush=True)
    return path


def attach_leiden(embed, full_h5ad, key=LEIDEN_KEY):
    """Add 05_2's leiden label to `embed.obs`, reading only `obs` out of the h5ad.

    The DRVI input is the 2,000-gene object, written by `reduce_data_tum.py` *before*
    `clustering_tum.py` ran, so the clustering is not in `adata.obs` and therefore not in
    the embedding either. It is worth going and getting: `cell_type` is constant here and
    `cell_type_01_4` is a CellTypist label from a normal-breast model, so leiden is the
    only grouping of these cells that was computed on these cells.

    Only the `obs` element is read (a few MB out of a ~330 MB file), the cells are
    realigned by name, and any failure is reported and swallowed - the run does not depend
    on this.
    """
    if key in embed.obs:
        return True
    full_h5ad = Path(full_h5ad)
    if not full_h5ad.exists():
        print(f"[skip] {key}: {full_h5ad.name} is not on disk", flush=True)
        return False
    try:
        import h5py
        try:
            from anndata.io import read_elem
        except ImportError:                      # anndata < 0.11
            from anndata.experimental import read_elem
        with h5py.File(full_h5ad, "r") as f:
            obs = read_elem(f["obs"])
        if key not in obs:
            print(f"[skip] {key}: not a column of {full_h5ad.name}", flush=True)
            return False
        # By name, never by position: a reordering upstream must not pair a cell with
        # another cell's cluster.
        missing = embed.obs_names.difference(obs.index)
        if len(missing):
            print(f"[skip] {key}: {len(missing):,} cells of the embedding are not in "
                  f"{full_h5ad.name}", flush=True)
            return False
        embed.obs[key] = obs.loc[embed.obs_names, key].values
        print(f"[have] {key}: {embed.obs[key].nunique()} clusters attached from "
              f"{full_h5ad.name}", flush=True)
        return True
    except Exception as exc:
        print(f"[skip] {key}: {type(exc).__name__}: {exc}", flush=True)
        return False


def groupable(embed, keys):
    """The keys of `keys` that are in `embed.obs` and have more than one level.

    A constant column is not an error here, it is the shape of this subset: `cell_type`
    is `malignant` for every cell. Passing it to a heatmap would draw one row and claim
    nothing, so it is dropped with a line in the log saying why.
    """
    out = []
    for key in keys:
        if key not in embed.obs:
            print(f"[skip] grouping {key!r}: not an obs column", flush=True)
        elif embed.obs[key].astype(str).nunique() < 2:
            print(f"[skip] grouping {key!r}: constant "
                  f"({embed.obs[key].astype(str).iloc[0]!r})", flush=True)
        else:
            out.append(key)
    return out


def draw_figures(model, adata, embed, fig_dir, run_id, heatmap_keys):
    """The plotting section of the notebook, in order and with the same names.

    Each UMAP panel has its direct counterpart in figures/05_2_reduce_data/umap_tum_*.png:
    the same cells in the *unintegrated* PCA space, with the same keys and the same
    palettes, so only the space changes between the two.
    """
    sc.settings.figdir = fig_dir

    # Same seeded permutation as the notebook, so overplotting hides the same cells in
    # both.
    order = np.random.default_rng(UMAP_SEED).permutation(embed.n_obs)
    embed_plot = ad.AnnData(
        obs=embed.obs.iloc[order].copy(),
        obsm={"X_umap": embed.obsm["X_umap"][order]},
        # Palettes carried over from the input object, so a label keeps the same colour
        # here and in the 05_2 UMAPs of the same cells.
        uns={k: v for k, v in adata.uns.items() if k.endswith("_colors")},
    )

    for label, col in UMAP_QC_KEYS.items():
        if col not in embed_plot.obs:
            print(f"[skip] umap {label}: no obs column {col!r}", flush=True)
            continue
        # A constant column paints one colour over 36,192 cells: dropped, with the reason.
        if embed_plot.obs[col].dtype.kind not in "fiu" and \
                embed_plot.obs[col].astype(str).nunique() < 2:
            print(f"[skip] umap {label}: {col!r} is constant "
                  f"({embed_plot.obs[col].astype(str).iloc[0]!r})", flush=True)
            continue
        sc.pl.umap(embed_plot, color=col, show=False, save=f"_{label}_{run_id}.png")
        plt.close("all")

    combined = [k for k in UMAP_COMBINED_KEYS
                if k in embed_plot.obs and embed_plot.obs[k].astype(str).nunique() > 1]
    with plt.rc_context({"figure.figsize": (7, 7)}):
        sc.pl.umap(embed_plot, color=combined, ncols=2, wspace=0.8,
                   hspace=0.25, show=False, save=f"_combined_{run_id}.png")
        plt.close("all")

    # One panel per level of each grouping, each highlighting only its own cells. Both
    # groupings get one, and they answer different questions: on the leiden partition,
    # where each state computed on these cells sits in the DRVI space; on
    # `cell_type_01_4`, which of DRVI's tumour states line up with what CellTypist had
    # called the cells *before* the CNV call - a landmark, not what kind of cell each one
    # is. The second is also the more crowded: Lumsec-prol and Lumsec-basal are 84% of the
    # subset between them.
    for key in PANEL_KEYS:
        if key not in embed_plot.obs or embed_plot.obs[key].astype(str).nunique() < 2:
            print(f"[skip] panels {key}: missing or constant", flush=True)
            continue
        groups = embed_plot.obs[key].astype("category")
        labels = [c for c in groups.cat.categories if (groups == c).any()]
        ncols = 4
        nrows = int(np.ceil(len(labels) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.4))
        axes = np.atleast_1d(axes).flatten()
        for ax, label in zip(axes, labels):
            n = int((groups == label).sum())
            sc.pl.umap(embed_plot, color=key, groups=[label], title=f"{label} ({n:,})",
                       ax=ax, show=False, legend_loc="none", na_in_legend=False, size=4)
            ax.set_xlabel("")
            ax.set_ylabel("")
        for ax in axes[len(labels):]:
            ax.set_visible(False)
        plt.tight_layout()
        savefig(fig_dir, run_id, f"umap_per_group_{key}", fig=fig)
    del embed_plot

    # Latent dimensions: per-dimension statistics, with and without the vanished ones.
    # The second is the plot behind the choice of n_latent.
    drvi.utils.pl.plot_latent_dimension_stats(embed, ncols=2, show=False)
    savefig(fig_dir, run_id, "latent_dimension_stats")
    drvi.utils.pl.plot_latent_dimension_stats(embed, ncols=2, remove_vanished=True, show=False)
    savefig(fig_dir, run_id, "latent_dimension_stats_rmVanished")

    # Every non-vanished dimension on the UMAP (remove_vanished defaults to True).
    drvi.utils.pl.plot_latent_dims_in_umap(embed, show=False)
    savefig(fig_dir, run_id, "latent_dims_in_umap")

    # Which dimension responds to which group. The first key gets the extra ordering
    # (grouped by category as well as by dimension): it makes it obvious when one label
    # owns a whole dimension.
    for i, key in enumerate(heatmap_keys):
        drvi.utils.pl.plot_latent_dims_in_heatmap(embed, key, title_col="title", show=False)
        savefig(fig_dir, run_id, f"latent_dims_in_heatmap_{key}")
        if i == 0:
            drvi.utils.pl.plot_latent_dims_in_heatmap(embed, key, title_col="title",
                                                      sort_by_categorical=True, show=False)
            savefig(fig_dir, run_id, f"latent_dims_in_heatmap_{key}_sorted")

    # Interpretability. OOD comes from the decoder reconstructions (fast, favours the
    # genes *specific* to a dimension); OOD_min/max_possible are its two halves; IND
    # averages each factor's effect over all cells (broader).
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
    os.environ["DATA_DIR"] = str(args.data_dir)   # cell_set reads it from the environment

    C.banner("05_3 DRVI run")

    # The compartment goes into the run id, so nothing here can be confused with 04's
    # `drvi_epi_*`, 03_2's `drvi_nonimm_*` or 02_2's whole-dataset `drvi_unscaled_*`, and
    # the latent size keeps the sizes side by side.
    run_id = f"drvi_{C.compartment()}_{args.n_latent}"

    phase_dir = Path(__file__).resolve().parent.parent   # 05_drvi_tumoral_epi/
    tum_dir = C.tum_dir()                                # $DATA_DIR/05_tum

    # Input: the 2,000-gene object written by 05_2's reduce_data_tum.py. HVGs were
    # selected on the malignant cells only, .layers['counts'] are the raw counts DRVI
    # trains on.
    input_h5ad = C.path("_hvg_2k.h5ad")
    # The definitive object of 05_2 (all genes, scran log-norm, leiden): read for the
    # leiden column and, at the end, to carry the latent space over for the steps after.
    full_h5ad = C.path(".h5ad")

    # Outputs, all under DATA_DIR and never in the repo (the embedding alone is a few
    # hundred MB, past GitHub's per-file limit, and /datasets/* is gitignored).
    model_path = tum_dir / f"model_{run_id}.pt"          # one file (see save_model)
    latent_h5ad = tum_dir / f"embed_{run_id}.h5ad"       # latent space + interpretability
    downstream_h5ad = C.path(f"_{run_id}.h5ad")          # 05_2 object + obsm['X_drvi']

    fig_dir = Path(args.fig_dir) if args.fig_dir else phase_dir / "figures" / f"05_3_{run_id}"
    fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"run id      {run_id}", flush=True)
    print(f"input       {input_h5ad}")
    print(f"model       {model_path}")
    print(f"latent      {latent_h5ad}")
    print(f"downstream  {downstream_h5ad}")
    print(f"figures     {fig_dir}")
    print(f"scvi-tools {scvi.__version__}, DRVI {drvi.__version__}", flush=True)

    if not input_h5ad.exists():
        sys.exit(f"missing {input_h5ad}: run 05_2_subsetting/subsetting_all.sh first")

    # Checked here rather than at the end, where it is used: the downstream object is
    # written after the training, and a missing input discovered there would only surface
    # hours in. Not fatal - the model and the embedding do not need it - but worth saying
    # before anything starts, especially on the cluster, where this file has to be copied
    # up as well.
    if not full_h5ad.exists():
        print(f"[warn] {full_h5ad} is missing: the run will produce the model and the "
              f"embedding,\n       but not {downstream_h5ad.name}, and the {LEIDEN_KEY} "
              f"figures will be skipped", flush=True)

    sc.settings.set_figure_params(dpi=300, facecolor="white")
    sc.settings.figdir = fig_dir

    adata = sc.read_h5ad(input_h5ad)
    print(f"{adata.n_obs:,} cells x {adata.n_vars:,} HVGs", flush=True)

    # Sanity checks on what 05_2 produced: the right cell set, the counts layer present.
    # `compartment` is what tells this object from 04's ('epi') and from 03's ('non_imm').
    expected = C.compartment()
    assert adata.obs[C.COMPARTMENT_KEY].astype(str).nunique() == 1, \
        "not a single-compartment object"
    assert adata.obs[C.COMPARTMENT_KEY].astype(str).iloc[0] == expected, \
        f"compartment is not {expected!r}: this is not the CELL_SET={C.cell_set()} object"
    assert "counts" in adata.layers, "raw counts layer missing, DRVI has nothing to train on"
    # No cell without a CNV call can be here: under `tum` by construction, under `epi`
    # because 05_2 drops them rather than assuming they are normal.
    assert (adata.obs[C.STATUS_KEY].astype(str) != "not_tested").all(), \
        "not_tested cells in the input: 05_2 should have dropped them"
    print("compartment :", adata.obs[C.COMPARTMENT_KEY].astype(str).unique().tolist())
    print("cnv_status  :", sorted(adata.obs[C.STATUS_KEY].astype(str).unique().tolist()))
    print("cell_type   :", sorted(adata.obs[LABEL_KEY].astype(str).unique().tolist()))
    print(f"{adata.obs[BATCH_KEY].nunique()} cohorts, "
          f"{adata.obs[GROUP_KEY].nunique()} pre-CNV labels", flush=True)

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

    # Train only if the model is not on disk yet, or if --overwrite asks for a retrain.
    # The training uses the GPU if one is available, otherwise the CPU.
    if args.overwrite or not model_path.exists():
        print(f">>> training {run_id}: {args.epochs} epochs, early stopping "
              f"{'on (patience 50)' if args.early_stopping else 'off'}", flush=True)
        model.train(
            max_epochs=args.epochs,
            early_stopping=args.early_stopping,  # off: the KL warmup needs all --epochs
            early_stopping_patience=50,          # only read when --early-stopping is given
        )
        save_model(model, model_path)
        print(f"[model] saved to {model_path}", flush=True)
    else:
        print(f"[have] {model_path} already exists, not retraining (--overwrite to force)",
              flush=True)

    # Load the model back, so a resumed run works from exactly the same state a fresh one
    # does.
    model = load_model(model_path, adata)

    # The DRVI-specific view of the run: the latent space as .X, one var per latent
    # dimension carrying its stats and interpretability scores. This is the file 05_4
    # reads.
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

    # The number the latent size is judged on: how many of the n_latent dimensions DRVI
    # did not need. This is what to read in the log before deciding whether to re-run at
    # another size - none vanished means the size is too tight and a larger one is worth a
    # run, which is what 04_2's first attempt at 32 reported (0 / 32) on twice these cells.
    n_vanished = int(embed.var["vanished"].sum())
    print(f"{n_vanished} vanished / {embed.n_vars} latent dimensions "
          f"({embed.n_vars - n_vanished} effectively used)", flush=True)
    if n_vanished == 0:
        print(f"[note] nothing vanished at n_latent={args.n_latent}: the space is too "
              f"tight to say what DRVI did not need. Worth a second run at "
              f"--n-latent {args.n_latent * 2}, which writes beside this one.", flush=True)

    # The clustering 05_2 computed on these very cells, for the heatmaps below. It is not
    # in the 2,000-gene input (clustering_tum.py runs after reduce_data_tum.py), so it is
    # fetched from the definitive object; failing to get it costs one figure and nothing
    # else.
    attach_leiden(embed, full_h5ad)
    heatmap_keys = groupable(embed, HEATMAP_KEYS)

    # The object for anything downstream that needs genes AND the latent space: the
    # definitive 05_2 object (all genes, scran log-norm, leiden) with the latent space
    # added as obsm['X_drvi']. The embedding above is the model's own view of the run -
    # one var per latent dimension, its stats and its interpretability scores - and by
    # construction carries no genes; this is the other half, the subset itself in the DRVI
    # space. 05_4 works from the embedding alone and does not open this file.
    if not full_h5ad.exists():
        print(f"[warn] {full_h5ad} is missing, so {downstream_h5ad.name} was not written.\n"
              f"       Run 05_2_subsetting/subsetting_all.sh (or copy that object here) and\n"
              f"       re-run: the model and the embedding above are already done, so this\n"
              f"       second pass only reads them back.", flush=True)
    elif args.overwrite or not downstream_h5ad.exists():
        print(f"[read] {full_h5ad}", flush=True)
        full = sc.read_h5ad(full_h5ad)

        # Cell order is not assumed: the latent space is realigned on the 05_2 object's
        # cells by name, so a reordering anywhere upstream cannot pair a cell with another
        # cell's coordinates.
        assert set(embed.obs_names) == set(full.obs_names), \
            "the embedding and the 05_2 object disagree on cells"
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
        draw_figures(model, adata, embed, fig_dir, run_id, heatmap_keys)
    except Exception:
        # Every artifact is already written and the training behind it is hours long:
        # report the failure, keep the exit status, and let the notebook redraw whichever
        # plot broke.
        print("[warn] the figures section failed; the artifacts above are complete",
              flush=True)
        traceback.print_exc()

    print(f"[ok] {run_id} -> {latent_h5ad}", flush=True)
    print("next: 05_4_signatures on that embedding", flush=True)


if __name__ == "__main__":
    main()
