#!/usr/bin/env python3
"""05_1 sensitivity 1 of 2, the analysis half: read the resolution sweep, apply the FULL
decision chain at each resolution, and report whether the per-patient verdict moves.

Applying the whole chain is the point. Counting subclusters would show only that the
clustering changed, which is not in dispute; what has to be shown is that the thing the
phase actually claims - this patient's epithelium is, or is not, CNV-positive - does not
depend on the parameter. The chain, identical to call_malignant_newcnv.ipynb:

  1. nullpos per cell   min(percentile on cnv_score, percentile on cnv_corr), both taken
                        against the stromal block OF THAT COHORT
  2. median per subcluster
  3. a subcluster is CNV-positive when that median >= MAL_MED
  4. the fraction of epithelium sitting in CNV-positive subclusters
  5. the cohort gate: the patient is CNV-positive when that fraction >= GATE_FRAC

MIN_SUB is reported rather than applied: it is the diagnostic that picks the operating
resolution. A subcluster of three cells has a median, and that median is not an aggregate
of anything - so the quantity to watch as resolution rises is not the verdict (which is
stable) but how much of the epithelium still sits in subclusters big enough to average.

Input : $DATA_DIR/05_tum/sensitivity/resolution/<cohort>_res<r>.csv  (sensitivity_resolution.sh)
Output: ../tables/05_1_infercnv/sensitivity/resolution_sweep.csv
        ../tables/05_1_infercnv/sensitivity/resolution_subclusters.csv
        ../figures/05_1_infercnv/sensitivity/resolution_sweep.png
"""
from __future__ import annotations

import glob
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MAL_MED = 0.90     # a subcluster whose median cell sits above this on the null -> positive
GATE_FRAC = 0.10   # ...and the cohort needs this fraction of its epithelium in such subclusters
MIN_SUB = 20       # subclusters below this are reported on, not excluded, see the docstring


def main() -> int:
    data_dir = os.path.abspath(os.path.expanduser(os.environ["DATA_DIR"]))
    in_dir = os.path.join(data_dir, "05_tum", "sensitivity", "resolution")
    here = os.path.dirname(os.path.abspath(__file__))
    phase = os.path.dirname(here)
    tab_dir = os.path.join(phase, "tables", "05_1_infercnv", "sensitivity")
    fig_dir = os.path.join(phase, "figures", "05_1_infercnv", "sensitivity")
    os.makedirs(tab_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)

    paths = sorted(glob.glob(os.path.join(in_dir, "*_res*.csv")))
    assert paths, f"no sweep output in {in_dir}; run ./sensitivity_resolution.sh first"

    rows, sub_rows = [], []
    for p in paths:
        m = re.match(r"(\w+)_res([\d.]+)\.csv", os.path.basename(p))
        cohort, res = m.group(1), float(m.group(2))
        d = pd.read_csv(p)
        stromal = d[d["group"] == "stromal"]
        for ax in ("cnv_score", "cnv_corr"):
            d["p_" + ax] = (np.searchsorted(np.sort(stromal[ax].to_numpy()),
                                            d[ax].to_numpy()) / len(stromal))
        d["nullpos"] = np.minimum(d["p_cnv_score"], d["p_cnv_corr"])

        epi = d[(d["group"] == "epi") & d["subcluster"].notna()]
        g = epi.groupby("subcluster")
        sc = pd.DataFrame({"n": g.size(), "median_nullpos": g["nullpos"].median()})
        sc["positive"] = sc["median_nullpos"] >= MAL_MED
        frac = epi["subcluster"].isin(sc.index[sc["positive"]]).mean()
        big = sc.index[sc["n"] >= MIN_SUB]

        held = d[d["group"] == "immune_heldout"]
        rows.append(dict(
            cohort=cohort, resolution=res, n_epi=len(epi),
            n_subclusters=len(sc), n_subclusters_ge_min=int((sc["n"] >= MIN_SUB).sum()),
            median_subcluster_size=int(sc["n"].median()), min_subcluster_size=int(sc["n"].min()),
            pct_cells_in_subclusters_ge_min=round(100 * epi["subcluster"].isin(big).mean(), 1),
            n_subclusters_positive=int(sc["positive"].sum()),
            frac_epi_positive=round(frac, 4),
            verdict="CNV-positive" if frac >= GATE_FRAC else "CNV-negative",
            median_nullpos_heldout=round(held["nullpos"].median(), 4) if len(held) else np.nan,
            median_nullpos_stromal=round(
                d.loc[d["group"] == "stromal", "nullpos"].median(), 4)))
        s2 = sc.reset_index()
        s2["cohort"], s2["resolution"] = cohort, res
        sub_rows.append(s2)

    r = pd.DataFrame(rows).sort_values(["cohort", "resolution"])
    r.to_csv(os.path.join(tab_dir, "resolution_sweep.csv"), index=False)
    pd.concat(sub_rows).to_csv(os.path.join(tab_dir, "resolution_subclusters.csv"), index=False)

    pd.set_option("display.width", 250)
    print(r.to_string(index=False))
    print()
    unstable = []
    for cohort, d in r.groupby("cohort"):
        spread = d["frac_epi_positive"].max() - d["frac_epi_positive"].min()
        stable = d["verdict"].nunique() == 1
        print(f"{cohort}: {'stable' if stable else 'UNSTABLE'} "
              f"{sorted(d['verdict'].unique())} | frac_epi_positive "
              f"{d['frac_epi_positive'].min():.3f}-{d['frac_epi_positive'].max():.3f} "
              f"(spread {spread:.3f})")
        if not stable:
            unstable.append(cohort)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for cohort, d in r.groupby("cohort"):
        d = d.sort_values("resolution")
        axes[0].plot(d["resolution"], d["frac_epi_positive"], "o-", label=cohort)
        axes[1].plot(d["resolution"], d["pct_cells_in_subclusters_ge_min"], "o-", label=cohort)
    axes[0].axhline(GATE_FRAC, color="k", ls="--", lw=.8)
    axes[0].text(0.005, GATE_FRAC, " cohort gate", fontsize=7, va="bottom")
    axes[0].set_ylabel("fraction of epithelium CNV-positive")
    axes[0].set_title("The verdict does not move")
    axes[1].set_ylabel(f"% of cells in subclusters of >= {MIN_SUB}")
    axes[1].set_title("...but the aggregation degrades")
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel("leiden resolution")
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "resolution_sweep.png"), dpi=200, bbox_inches="tight")
    print(f"\nWrote {tab_dir}/resolution_sweep.csv and the figure")
    return 1 if unstable else 0


if __name__ == "__main__":
    sys.exit(main())
