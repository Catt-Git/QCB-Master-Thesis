#!/usr/bin/env python3
"""05_1 sensitivity 2 of 2: does the call depend on which population is the baseline?

The objection this answers is "the call depends on the population you chose as reference".
Four reference configurations are run on the SAME observations and the SAME gene set, so
the only thing that differs between them is the baseline:

  A_current   ref_tcell + ref_myeloid       the design in production
  B_TandB     ref_tcell + ref_bcell         B cells instead of myeloid
  C_fine      up to 6 immune subtypes       one group per cell type, the finest split
  D_pooled    one pooled immune group       no band at all - the lower bound

The four are not arbitrary. inferCNV takes the residual against the BOUNDS of the
per-group means (min and max), so anything an observed cell has inside that band becomes
exactly zero: the number of reference groups is a conservativeness dial, and these four
span it from no protection (D) to the most the data supports (C).

Two cohorts, chosen as the extremes the call has to survive:
  Patient16   tumour-rich, 8,028 epithelial cells   - can the split find the tumour?
  Patient52   no detectable tumour, 9,243 cells     - can it invent one?

Everything is written under $DATA_DIR/05_tum/sensitivity/reference/; the production
input/ and summary/ trees are never touched.

  --mode prepare   build the four input sets   (then sensitivity_reference.sh runs them)
  --mode analyze   read the results, write the tables

Result as of the run in the README: the AUC separating epithelium from the stromal null
moves by 0.003 across all four configurations on Patient16 (0.972-0.975) and stays at
chance on Patient52 (0.452-0.525). The reference composition is not a lever on this call.
What the experiment did establish is the held-out immune control: 0.00% of those cells
are called, in eight runs out of eight.
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from prepare_infercnv_input import (  # noqa: E402
    EPITHELIAL, STROMAL_STRATA, T_NK, MYELOID, B_PLASMA,
    MIN_CELLS_PER_GENE, SEED, stratified_sample,
)

COHORTS = ("Patient16", "Patient52")
FINE = ("CD8-Tem", "CD4-Treg", "CD4-naive", "NK-ILCs", "Macro-lipo", "Mono-classical")
N_REF, N_STROMAL, N_HELD, MIN_FINE = 1000, 2000, 500, 100
Q = 0.95


def auc(x: np.ndarray, y: np.ndarray) -> float:
    """P(x > y), rank-based, so it does not depend on each run's own scale."""
    r = pd.concat([pd.Series(x), pd.Series(y)]).rank().to_numpy()
    n1 = len(x)
    return (r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * len(y))


def prepare(data_dir: Path) -> None:
    import anndata as ad
    import scipy.io as sio
    import scipy.sparse as sp

    rng = np.random.default_rng(SEED)
    out_root = data_dir / "05_tum" / "sensitivity" / "reference"
    adata = ad.read_h5ad(data_dir / "shiao.h5ad", backed="r")

    go = pd.read_csv(data_dir / "05_tum" / "gene_order_hg38_gencode_v27.txt", sep="\t",
                     header=None, names=["gene", "chr", "start", "stop", "_t"],
                     usecols=[0, 1, 2, 3])
    main = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
    go = go[go["chr"].isin(main)]
    go["chr"] = pd.Categorical(go["chr"], categories=main, ordered=True)
    go = go.sort_values(["chr", "start", "stop"]).set_index("gene")
    keep = pd.Index(go.index[go.index.isin(adata.var_names.intersection(go.index))])
    gpos = pd.Series(np.arange(adata.n_vars), index=adata.var_names)[keep].to_numpy()

    for cohort in COHORTS:
        m = (adata.obs["cohort"].astype(str) == cohort).to_numpy()
        pos = np.flatnonzero(m)
        lab = adata.obs["cell_type"].astype(str)[m]
        raw = adata.obs["celltypist_predicted"].astype(str)[m]
        pure = ~raw.isin(EPITHELIAL).to_numpy()

        epi = np.flatnonzero(lab.isin(EPITHELIAL).to_numpy())
        strata = {k: np.flatnonzero(lab.isin(v).to_numpy()) for k, v in STROMAL_STRATA.items()}
        stromal = stratified_sample(strata, N_STROMAL, rng)
        imm = np.flatnonzero(lab.isin(T_NK | MYELOID).to_numpy() & pure)
        held = np.sort(rng.choice(imm, size=min(N_HELD, len(imm)), replace=False))
        avail = np.setdiff1d(imm, held)

        def take(mask_labels, cap, pool=avail):
            ix = np.intersect1d(pool, np.flatnonzero(lab.isin(mask_labels).to_numpy() & pure))
            if len(ix) > cap:
                ix = rng.choice(ix, size=cap, replace=False)
            return np.sort(ix)

        configs = {
            "A_current": {"ref_tcell": take(T_NK, N_REF), "ref_myeloid": take(MYELOID, N_REF)},
            "B_TandB": {"ref_tcell": take(T_NK, N_REF),
                        "ref_bcell": take(B_PLASMA, N_REF, pool=np.arange(m.sum()))},
            "D_pooled": {"ref_immune": np.sort(rng.choice(avail, size=min(2 * N_REF, len(avail)),
                                                          replace=False))},
        }
        fine = {f"ref_{c}": take({c}, N_REF) for c in FINE}
        configs["C_fine"] = {k: v for k, v in fine.items() if len(v) >= MIN_FINE}

        # One gene filter, computed on the OBSERVATIONS only, so all four configurations
        # are run on exactly the same features and the comparison is about the reference.
        obs_rows = np.sort(pos[np.concatenate([epi, stromal, held])])
        X = sp.csr_matrix(adata.layers["counts"][obs_rows, :])[:, gpos]
        gmask = np.asarray((X > 0).sum(axis=0)).ravel() >= MIN_CELLS_PER_GENE
        genes = keep[gmask]

        for name, ref in configs.items():
            blocks = ([("epi", epi), ("stromal", stromal), ("immune_heldout", held)]
                      + [(k, v) for k, v in ref.items()])
            loc = np.concatenate([b[1] for b in blocks])
            grp = np.concatenate([[b[0]] * len(b[1]) for b in blocks])
            rows = pos[loc]
            srt = np.argsort(rows, kind="stable")
            Xc = sp.csr_matrix(adata.layers["counts"][rows[srt], :])[:, gpos][:, gmask]
            d = out_root / cohort / name
            d.mkdir(parents=True, exist_ok=True)
            sio.mmwrite(str(d / "counts.mtx"), Xc.T.astype(np.int32), field="integer")
            pd.Series(genes).to_csv(d / "genes.tsv", index=False, header=False)
            pd.Series(adata.obs_names[rows[srt]]).to_csv(d / "barcodes.tsv", index=False,
                                                         header=False)
            pd.DataFrame({"c": adata.obs_names[rows[srt]], "g": grp[srt]}).to_csv(
                d / "annotations.tsv", sep="\t", index=False, header=False)
            pd.Series(sorted(ref)).to_csv(d / "ref_groups.txt", index=False, header=False)
            print(f"[ok] {cohort}/{name}: {Xc.shape[0]:,} cells x {Xc.shape[1]:,} genes | "
                  + ", ".join(f"{k}={len(v)}" for k, v in ref.items()), flush=True)


def analyze(data_dir: Path) -> None:
    here = Path(os.path.abspath(__file__)).parent
    tab = here.parent / "tables" / "05_1_infercnv" / "sensitivity"
    tab.mkdir(parents=True, exist_ok=True)
    rows = []
    for p in sorted(glob.glob(str(data_dir / "05_tum" / "sensitivity" / "reference"
                                  / "*" / "*_cells.csv"))):
        cohort = Path(p).parent.name
        cfg = re.sub(r"_cells\.csv$", "", Path(p).name)
        d = pd.read_csv(p)
        e = d.loc[d["group"] == "epi", "cnv_score"].to_numpy()
        s = d.loc[d["group"] == "stromal", "cnv_score"].to_numpy()
        h = d.loc[d["group"] == "immune_heldout", "cnv_score"].to_numpy()
        cs = np.quantile(s, Q)
        cc = np.quantile(d.loc[d["group"] == "stromal", "cnv_corr"], Q)
        hit = (d["cnv_score"] > cs) & (d["cnv_corr"] > cc)
        rows.append(dict(
            cohort=cohort, config=cfg,
            n_ref_groups=d.loc[d["group"].str.startswith("ref_"), "group"].nunique(),
            median_epi=round(float(np.median(e)), 6),
            median_stromal=round(float(np.median(s)), 6),
            median_heldout=round(float(np.median(h)), 6),
            contrast_epi_over_stromal=round(float(np.median(e) / np.median(s)), 2),
            auc_epi_vs_stromal=round(auc(e, s), 3),
            auc_epi_vs_heldout=round(auc(e, h), 3),
            auc_stromal_vs_heldout=round(auc(s, h), 3),
            pct_epi_called=round(100 * hit[d["group"] == "epi"].mean(), 1),
            pct_heldout_called=round(100 * hit[d["group"] == "immune_heldout"].mean(), 2)))
    r = pd.DataFrame(rows).sort_values(["cohort", "config"])
    r.to_csv(tab / "reference_configurations.csv", index=False)
    pd.set_option("display.width", 250)
    print(r.to_string(index=False))
    print()
    for cohort, d in r.groupby("cohort"):
        a = d["auc_epi_vs_stromal"]
        print(f"{cohort}: AUC epi-vs-stromal {a.min():.3f}-{a.max():.3f} "
              f"(spread {a.max()-a.min():.3f}) | held-out called "
              f"{d['pct_heldout_called'].min():.2f}-{d['pct_heldout_called'].max():.2f}%")
    print(f"\nWrote {tab / 'reference_configurations.csv'}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--mode", choices=("prepare", "analyze"), required=True)
    args = ap.parse_args()
    data_dir = Path(os.environ["DATA_DIR"]).expanduser().resolve()
    (prepare if args.mode == "prepare" else analyze)(data_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
