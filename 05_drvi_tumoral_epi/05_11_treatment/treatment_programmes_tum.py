#!/usr/bin/env python3
"""05_11, idea 1: do the NAMED programmes move with the treatment?

The DRVI directions 05_8 named (the bars of `programmes_named`) are read per biopsy and
compared WITHIN the patient across the three timepoints of the trial. Nothing is retrained:
the latent space is 05_3's, the programme list is 05_8's, and the only new thing is the
grouping of cells by `treatment`.

WHAT IS MEASURED, PER BIOPSY (patient x timepoint, >= --min-cells malignant cells):

  direction   the coordinate on one named direction, oriented so that higher means more of
              the programme, z-scored over all cells - see `treatment_common`
  programme   the mean of a programme's directions: ONE readout per programme, because a
              programme named on four directions is one hypothesis, not four
  gene level  for the a priori programmes only, 05_6's raw `score_genes` value of the
              signature itself: the same question asked of the genes, with no latent space in
              between. A DR that moves while its genes do not is a property of the
              embedding, and this column is where that would show

  For each, the biopsy MEAN (prevalence and intensity together) is tested; the fraction of
  cells above 1 sd is written beside it, to read which of the two moved, and not tested.

WHAT IS TESTED. For each of the three contrasts (PD1 vs BASE, RTPD1 vs BASE, RTPD1 vs PD1),
the per-patient difference of the biopsy means, Wilcoxon signed-rank over the patients that
have both biopsies. Pairs, not pooled cells: the patient is the unit, and the patient's own
baseline is the control. Fewer than MIN_PAIRS pairs -> not tested, and said so.

TWO FAMILIES, KEPT APART. The hypotheses in `treatment_common.HYPOTHESES` were written before
the data were looked at; they are a handful of tests and get their own Benjamini-Hochberg. Every
other (readout, contrast) is exploratory and corrected as one family. A hit in the second is a
lead, never a result, whatever its FDR.

RESPONSE IS DESCRIBED, NOT TESTED AS AN INTERACTION. With ~4-8 patients per response group the
timepoint x response interaction has no power worth the name. The median delta per group is
written, with its n, and a Kruskal-Wallis across groups is reported as exploratory where each
group has >= 3 pairs. The figures draw every patient, so the reader sees the n.

Usage:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    N_LATENT=64 HVG_SET=nomt python3 treatment_programmes_tum.py
    N_LATENT=64 HVG_SET=nomt python3 treatment_programmes_tum.py --min-cells 30
    # another naming table (the defaults are PRUNE_VANISHED=1 OUT_TAG=_pruned_rho030):
    N_LATENT=64 HVG_SET=nomt PRUNE_VANISHED=1 OUT_TAG=_pruned python3 treatment_programmes_tum.py
"""

from __future__ import annotations

import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from scipy.stats import kruskal, wilcoxon

import treatment_common as T
from treatment_common import C

# Cells above this oriented z are "high" for the descriptive fraction. One sd is a plain
# reading, not a threshold anything is decided on.
HIGH_Z = 1.0

FDR_ALPHA = 0.05


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--min-cells", type=int, default=T.MIN_CELLS_PER_SAMPLE,
                   help=f"malignant cells a biopsy needs to be summarised "
                        f"(default {T.MIN_CELLS_PER_SAMPLE})")
    p.add_argument("--no-gene-level", action="store_true",
                   help="skip the gene-level check of the a priori programmes (it needs "
                        "05_6's per-cell score tables in $DATA_DIR/05_tum/)")
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


# --------------------------------------------------------------------------- #
# Pseudobulk
# --------------------------------------------------------------------------- #


def pseudobulk(R: pd.DataFrame, obs: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    """Long table: one row per (kept biopsy, readout), mean and fraction above HIGH_Z."""
    sid = T.sample_ids(obs)
    keep = sid.isin(samples.index[samples["kept"]])
    R, sid = R.loc[keep], sid.loc[keep]
    mean = R.groupby(sid.values).mean()
    frac = (R > HIGH_Z).groupby(sid.values).mean()
    long = (mean.stack().rename("mean").to_frame()
            .join(frac.stack().rename("frac_high")))
    long.index.names = ["sample", "readout"]
    long = long.reset_index().join(samples[["cohort", "treatment", "response", "n_cells"]],
                                   on="sample")
    return long


def paired(pb: pd.DataFrame, readout: str, later: str, earlier: str):
    """Per-patient delta (later - earlier) of the biopsy means, with the patient's response."""
    d = pb[pb["readout"] == readout]
    a = d[d["treatment"] == later].set_index("cohort")["mean"]
    b = d[d["treatment"] == earlier].set_index("cohort")["mean"]
    common = a.index.intersection(b.index)
    resp = d.drop_duplicates("cohort").set_index("cohort")["response"]
    return pd.DataFrame({"delta": a.loc[common] - b.loc[common],
                         "response": resp.loc[common]})


def test_readout(pb: pd.DataFrame, readout: str, kind: str) -> list[dict]:
    rows = []
    for later, earlier in T.CONTRASTS:
        pr = paired(pb, readout, later, earlier)
        n = len(pr)
        rec = {"readout": readout, "kind": kind,
               "contrast": T.contrast_name(later, earlier), "n_pairs": n,
               "median_delta": float(pr["delta"].median()) if n else np.nan,
               "mean_delta": float(pr["delta"].mean()) if n else np.nan,
               "n_up": int((pr["delta"] > 0).sum()), "p": np.nan}
        if n >= T.MIN_PAIRS and (pr["delta"] != 0).any():
            rec["p"] = float(wilcoxon(pr["delta"]).pvalue)
        groups = []
        for g in T.RESPONSES:
            x = pr.loc[pr["response"] == g, "delta"]
            rec[f"n_{g}"] = len(x)
            rec[f"median_delta_{g}"] = float(x.median()) if len(x) else np.nan
            if len(x) >= 3:
                groups.append(x.to_numpy())
        rec["p_response_kruskal"] = (float(kruskal(*groups).pvalue)
                                     if len(groups) >= 2 else np.nan)
        rows.append(rec)
    return rows


def run_tests(pb: pd.DataFrame, prog: pd.DataFrame, directions: list[str]) -> pd.DataFrame:
    rows = []
    for p in prog.index:
        rows += test_readout(pb, p, "programme")
    for d in directions:
        rows += test_readout(pb, d, "direction")
    res = pd.DataFrame(rows)

    hyp = T.HYPOTHESES.set_index(["programme", "contrast"])
    key = list(zip(res["readout"], res["contrast"]))
    res["family"] = ["a_priori" if (k in hyp.index and kind == "programme") else "exploratory"
                     for k, kind in zip(key, res["kind"])]
    res["expected_sign"] = [hyp.loc[k, "expected_sign"] if f == "a_priori" else np.nan
                            for k, f in zip(key, res["family"])]
    res["fdr"] = np.nan
    for fam, idx in res.groupby("family").groups.items():
        res.loc[idx, "fdr"] = T.bh(res.loc[idx, "p"])
    res["direction_as_expected"] = np.where(
        res["family"] == "a_priori",
        np.sign(res["median_delta"]) == res["expected_sign"], np.nan)
    missing = set(map(tuple, T.HYPOTHESES[["programme", "contrast"]].values)) - set(
        k for k, f in zip(key, res["family"]) if f == "a_priori")
    for programme, contrast in sorted(missing):
        print(f"[note] a priori {programme} / {contrast}: the programme is not among the "
              f"bars of this naming table, so the hypothesis is not tested")
    return res


# --------------------------------------------------------------------------- #
# Gene-level check of the a priori programmes
# --------------------------------------------------------------------------- #


def gene_level(obs, samples, prog) -> pd.DataFrame | None:
    """05_6's raw signature score, pseudobulked and tested like the DR readouts.

    Raw `score_` and not the within-cohort `z_`: the paired difference already removes the
    patient, and the z would rescale each patient by its own spread.
    """
    wanted = [p for p in T.HYPOTHESES["programme"].unique() if p in prog.index]
    cols = {}
    for p in wanted:
        for slug in prog.loc[p, "collections"]:
            path = C.scores_csv(T.collection_by_slug(slug))
            if not path.exists():
                print(f"[skip] gene level for {p}: missing {path.name} (05_6 --collection "
                      f"{slug})")
                continue
            header = pd.read_csv(path, comment="#", nrows=0).columns
            col = f"score_{p}"
            if col not in header:
                continue
            s = pd.read_csv(path, comment="#", index_col=0, usecols=[header[0], col])[col]
            cols[p] = s
            break
    if not cols:
        return None
    G = pd.DataFrame(cols)
    missing = obs.index.difference(G.index)
    if len(missing):
        print(f"[skip] gene level: {len(missing)} cells of the embedding have no score row")
        return None
    G = G.loc[obs.index]
    pb = pseudobulk(G, obs, samples)
    rows = []
    for p in G.columns:
        rows += test_readout(pb, p, "gene_level")
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #


def fig_counts(counts: pd.DataFrame, min_cells: int, dpi: int) -> None:
    """Malignant cells per biopsy: the denominator of every other figure of the step."""
    m = counts[list(T.TREATMENTS)].to_numpy(float)
    fig, ax = plt.subplots(figsize=(4.8, 0.28 * len(counts) + 1.6))
    from matplotlib.colors import LinearSegmentedColormap, LogNorm
    cmap = LinearSegmentedColormap.from_list("blue", ["#cde2fb", "#6da7ec", "#1c5cab"])
    norm = LogNorm(vmin=max(min_cells, 1), vmax=max(m.max(), min_cells + 1))
    ax.imshow(np.where(m >= min_cells, m, np.nan), cmap=cmap, aspect="auto", norm=norm)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            below = m[i, j] < min_cells
            if below:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, color="#f0efec", lw=0))
            ax.text(j, i, f"{int(m[i, j])}", ha="center", va="center", fontsize=7,
                    color="0.55" if below else ("white" if norm(m[i, j]) > 0.6 else "0.1"))
    ax.set_xticks(range(len(T.TREATMENTS)), T.TREATMENTS, fontsize=8)
    ax.set_yticks(range(len(counts)),
                  [f"{c}  ({r})" for c, r in zip(counts.index, counts["response"])],
                  fontsize=7)
    for t, r in zip(ax.get_yticklabels(), counts["response"]):
        t.set_color(T.RESPONSE_COLORS.get(r, "0.2"))
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(f"Malignant cells per biopsy\n(grey: below the {min_cells}-cell floor, "
                 f"not summarised)", fontsize=9, loc="left")
    fig.tight_layout()
    C.savefig("sample_cell_counts", T.STEP, T.CONSENSUS, fig=fig, dpi=dpi)
    plt.close(fig)


def draw_trajectory(ax, pb: pd.DataFrame, readout: str) -> None:
    """One line per patient across the timepoints, coloured by response; group medians bold."""
    d = pb[pb["readout"] == readout]
    x_of = {t: i for i, t in enumerate(T.TREATMENTS)}
    for cohort, g in d.groupby("cohort"):
        g = g.assign(x=g["treatment"].map(x_of)).sort_values("x")
        ax.plot(g["x"], g["mean"], color=T.RESPONSE_COLORS[g["response"].iloc[0]],
                lw=0.8, alpha=0.35, marker="o", ms=2.5, zorder=1)
    for k, resp in enumerate(T.RESPONSES):
        g = d[d["response"] == resp]
        med = g.groupby("treatment")["mean"].median().reindex(T.TREATMENTS)
        off = (k - 1) * 0.07
        ax.plot(np.arange(len(T.TREATMENTS)) + off, med.values, color=T.RESPONSE_COLORS[resp],
                lw=2.2, marker="o", ms=5.5, mec="white", mew=1.0, zorder=3)
    ax.axhline(0, color="0.85", lw=0.8, zorder=0)
    ax.set_xticks(range(len(T.TREATMENTS)), T.TREATMENTS, fontsize=7.5)
    ax.set_xlim(-0.3, len(T.TREATMENTS) - 0.7)
    ax.tick_params(axis="y", labelsize=7)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color("0.7")


def response_legend(fig, pb: pd.DataFrame) -> None:
    n = pb.drop_duplicates("cohort")["response"].value_counts()
    handles = [Line2D([], [], color=T.RESPONSE_COLORS[r], lw=2.2, marker="o", ms=5,
                      label=f"{r} (n={n.get(r, 0)} patients)") for r in T.RESPONSES]
    handles.append(Line2D([], [], color="0.5", lw=0.8, alpha=0.5, label="one patient"))
    fig.legend(handles=handles, loc="upper right", ncol=len(handles), fontsize=7.5,
               frameon=False, bbox_to_anchor=(1.0, 1.0))


def fmt_test(r) -> str:
    if pd.isna(r["p"]):
        return f"n={r['n_pairs']} pairs, not tested"
    return (f"{r['contrast'].replace('_', ' ')}: median {r['median_delta']:+.2f} sd, "
            f"p={r['p']:.2g}, FDR={r['fdr']:.2g} (n={r['n_pairs']})")


def gene_line(r) -> str:
    """The gene-level check under the DR result, when 05_6's scores were there to run it."""
    if "gene_p" not in r or pd.isna(r.get("gene_p")):
        return ""
    return (f"\ngenes (05_6 score): p={r['gene_p']:.2g}, FDR={r['gene_fdr']:.2g}, "
            f"median {'up' if r['gene_median_delta'] > 0 else 'down'}")


def fig_trajectories(pb, res, readouts, name, title, dpi, hyp=None) -> None:
    n = len(readouts)
    ncols = 3 if hyp is not None else 4
    nrows = int(np.ceil(n / ncols))
    # The a priori panels carry three lines of test text underneath, hence the extra height.
    height = 2.9 * nrows + 0.6 + (0.4 * nrows if hyp is not None else 0.0)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.6 * ncols, height), squeeze=False)
    for ax, ro in zip(axes.flat, readouts):
        draw_trajectory(ax, pb, ro)
        ax.set_title(ro.replace("_", " "), fontsize=8.5, loc="left")
        if hyp is not None:
            h = hyp[hyp["programme"] == ro].iloc[0]
            r = res[(res["readout"] == ro) & (res["contrast"] == h["contrast"])].iloc[0]
            ax.text(0.0, -0.22, f"{h['rationale']} (expected {'up' if h['expected_sign'] > 0 else 'down'})\n"
                    + "DR: " + fmt_test(r) + gene_line(r), transform=ax.transAxes,
                    fontsize=6.5, color="0.3", va="top")
        ax.set_ylabel("biopsy mean, sd", fontsize=7)
    for ax in list(axes.flat)[n:]:
        ax.set_visible(False)
    response_legend(fig, pb)
    fig.suptitle(title, x=0.01, ha="left", fontsize=11, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.95), h_pad=4.0 if hyp is not None else 1.2)
    T.footnote(fig)
    C.savefig(name, T.STEP, T.CONSENSUS, fig=fig, dpi=dpi)
    plt.close(fig)


def fig_delta_heatmap(res: pd.DataFrame, kind: str, row_labels: dict, name: str, title: str,
                      dpi: int) -> None:
    d = res[res["kind"] == kind]
    order = list(dict.fromkeys(d["readout"]))
    M = d.pivot(index="readout", columns="contrast", values="median_delta") \
        .reindex(index=order, columns=list(T.CONTRAST_NAMES))
    F = d.pivot(index="readout", columns="contrast", values="fdr") \
        .reindex(index=order, columns=list(T.CONTRAST_NAMES))
    A = d.pivot(index="readout", columns="contrast", values="family") \
        .reindex(index=order, columns=list(T.CONTRAST_NAMES))
    lim = float(np.nanmax(np.abs(M.to_numpy()))) if np.isfinite(M.to_numpy()).any() else 1.0
    fig, ax = plt.subplots(figsize=(5.4, 0.3 * len(order) + 1.8))
    im = ax.imshow(M.to_numpy(float), cmap=T.diverging_cmap(), vmin=-lim, vmax=lim,
                   aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v, f = M.iat[i, j], F.iat[i, j]
            if np.isnan(v):
                ax.text(j, i, "n.t.", ha="center", va="center", fontsize=6, color="0.55")
                continue
            star = "*" if (not np.isnan(f) and f < FDR_ALPHA) else ""
            ax.text(j, i, f"{v:+.2f}{star}", ha="center", va="center", fontsize=6.5,
                    color="white" if abs(v) > 0.6 * lim else "0.15")
            if A.iat[i, j] == "a_priori":
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="0.1",
                                       lw=1.4))
    ax.set_xticks(range(M.shape[1]), [c.replace("_", " ") for c in M.columns], fontsize=7.5)
    ax.set_yticks(range(M.shape[0]), [row_labels.get(r, r) for r in M.index], fontsize=7)
    ax.tick_params(length=0)
    ax.xaxis.tick_top()
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
    cb.set_label("median paired delta (sd)", fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    cb.outline.set_visible(False)
    ax.set_title(title, fontsize=9, loc="left", pad=28)
    fig.text(0.01, 0.005, f"* FDR < {FDR_ALPHA} within its family. Outlined: hypothesis "
             f"fixed in advance (own BH family); everything else is exploratory. "
             f"n.t.: fewer than {T.MIN_PAIRS} patient pairs.",
             fontsize=6, color="0.35", ha="left", va="bottom", wrap=True)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    T.footnote(fig)
    C.savefig(name, T.STEP, T.CONSENSUS, fig=fig, dpi=dpi)
    plt.close(fig)


# --------------------------------------------------------------------------- #


def main():
    args = parse_args()
    C.banner("05_11 - the named programmes across the treatment timepoints")
    print(f"run {C.RUN_ID}; programme table from {C.TABLE_DIR.name}/")

    adata, L = T.load_embedding()
    obs = adata.obs
    prog = T.read_programmes()

    samples = T.sample_table(obs, args.min_cells)
    T.report_samples(samples, args.min_cells)
    counts = T.count_matrix(samples)
    C.write_table(counts, "treatment_cell_counts", T.CONSENSUS)
    fig_counts(counts, args.min_cells, args.dpi)

    directions = sorted({d for ds in prog["directions"] for d in ds}, key=C.dim_sort_key)
    Z = T.oriented_directions(L, directions)
    P = T.programme_composites(Z, prog)
    R = pd.concat([P, Z], axis=1)

    pb = pseudobulk(R, obs, samples)
    pb["kind"] = np.where(pb["readout"].isin(prog.index), "programme", "direction")
    C.write_table(pb, "treatment_pseudobulk", T.CONSENSUS, index=False)

    res = run_tests(pb, prog, directions)
    if not args.no_gene_level:
        gl = gene_level(obs, samples, prog)
        if gl is not None:
            keep = ["n_pairs", "median_delta", "p"] + [f"median_delta_{g}"
                                                       for g in T.RESPONSES]
            gl = gl[["readout", "contrast"] + keep].rename(
                columns={c: f"gene_{c}" for c in keep})
            res = res.merge(gl, on=["readout", "contrast"], how="left")
            # Same family as the DR test it checks: the a priori rows, BH over those alone.
            # The deltas are in raw score_genes units, not sd - read their sign, not their size.
            ap = res["family"] == "a_priori"
            res["gene_fdr"] = np.nan
            res.loc[ap, "gene_fdr"] = T.bh(res.loc[ap, "gene_p"])
    C.write_table(res, "treatment_paired_tests", T.CONSENSUS, index=False)

    C.banner("a priori hypotheses (own BH family)")
    cols = ["readout", "contrast", "n_pairs", "median_delta", "p", "fdr",
            "direction_as_expected"] + [f"median_delta_{g}" for g in T.RESPONSES]
    cols += [c for c in ("gene_median_delta", "gene_p", "gene_fdr") if c in res.columns]
    print(res[res["family"] == "a_priori"][cols].to_string(index=False, float_format="%.3g"))
    C.banner(f"exploratory family, FDR < {FDR_ALPHA} (leads, not results)")
    hits = res[(res["family"] == "exploratory") & (res["fdr"] < FDR_ALPHA)]
    print(hits[["readout", "kind", "contrast", "n_pairs", "median_delta", "fdr"]]
          .sort_values("fdr").to_string(index=False, float_format="%.3g")
          if len(hits) else "none")

    # Figures. Row labels of the direction heatmap name the programmes each direction carries.
    carries = {}
    for p, r in prog.iterrows():
        for d in r["directions"]:
            carries.setdefault(d, []).append(p)
    dir_labels = {d: f"{d}  ({', '.join(carries[d])})" for d in directions}

    hyp = T.HYPOTHESES[T.HYPOTHESES["programme"].isin(prog.index)]
    if len(hyp):
        fig_trajectories(pb, res, list(dict.fromkeys(hyp["programme"])),
                         "a_priori_trajectories",
                         "Hypotheses fixed in advance: the programme per biopsy, by patient",
                         args.dpi, hyp=hyp)
    fig_trajectories(pb, res, list(prog.index), "programme_trajectories",
                     f"Every named programme per biopsy ({C.RUN_ID})", args.dpi)
    fig_delta_heatmap(res, "programme", {}, "paired_delta_programmes",
                      "Paired change per programme (later minus earlier biopsy)", args.dpi)
    fig_delta_heatmap(res, "direction", dir_labels, "paired_delta_directions",
                      "Paired change per named direction (exploratory)", args.dpi)

    print(f"\n[done] {C.fig_dir(T.STEP, T.CONSENSUS)}")


if __name__ == "__main__":
    main()
