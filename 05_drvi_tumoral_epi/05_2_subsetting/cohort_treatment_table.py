"""Cells per cohort x treatment, as a figure, for both cell sets.

The table is the `cohort_x_treatment` crosstab of subset_and_qc.ipynb, but read from the object
that notebook WRITES (`<prefix>_raw.h5ad`): post-QC, post gene filter and post cohort drop, so
the cohorts on the y-axis are the ones every later step of the phase actually sees.

One heatmap per cell set, counts printed in every cell. The colour scale is logarithmic because
the counts span four orders of magnitude (a cohort can have 7,100 cells in one timepoint and 2 in
another): on a linear scale every cell under ~500 is the same white and the table is unreadable.
Zeros are drawn in grey and labelled, since "no cells in this timepoint" is the thing a reader
looks for (04_1's completeness filter was about exactly that).

Only `.obs` is read (backed mode), so this runs in seconds and needs no memory.

Usage:
    DATA_DIR=/path/to/datasets python cohort_treatment_table.py
"""
from __future__ import annotations

import os
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LogNorm

import cell_set as C

PHASE_DIR = Path(__file__).resolve().parent.parent
FIG_DIR = PHASE_DIR / "figures" / "05_2_subset_and_qc"
FIG_DIR.mkdir(parents=True, exist_ok=True)

TITLES = {"tum": "malignant", "epi": "malignant + normal epithelium"}
ZERO_COLOR = "#E6E6E6"


def cohort_x_treatment(value: str) -> pd.DataFrame:
    obs = ad.read_h5ad(C.path("_raw.h5ad", value), backed="r").obs
    table = pd.crosstab(obs[C.BATCH_KEY], obs[C.TREATMENT_KEY])
    table = table.reindex(columns=list(C.REQUIRED_TREATMENTS), fill_value=0)
    return table.loc[table.sum(axis=1) > 0].sort_index()


def plot(table: pd.DataFrame, value: str) -> Path:
    n_rows = len(table)
    counts = table.to_numpy()

    fig, ax = plt.subplots(figsize=(4.6, 0.34 * n_rows + 1.6))
    sns.heatmap(
        table.where(table > 0),                     # zeros -> NaN -> drawn as the axes face
        ax=ax, cmap="Blues", norm=LogNorm(vmin=1, vmax=counts.max()),
        annot=counts, fmt=",d", annot_kws={"fontsize": 8},
        linewidths=1.5, linecolor="white", cbar_kws={"label": "cells (log scale)", "shrink": 0.6},
    )
    ax.set_facecolor(ZERO_COLOR)
    # seaborn does not annotate NaN cells, so the zeros get their label by hand.
    for i, j in zip(*np.nonzero(counts == 0)):
        ax.text(j + 0.5, i + 0.5, "0", ha="center", va="center", fontsize=8, color="#555555")
    # Per-patient total in the tick label, not as a fourth column: on the same log scale it
    # would be the darkest cell of every row and pull the eye away from the split.
    ax.set_yticklabels([f"{c}  ({t:,})" for c, t in zip(table.index, table.sum(axis=1))])

    ax.set_title(
        f"Cells per patient and treatment - {TITLES[value]}\n"
        f"{int(counts.sum()):,} cells, {n_rows} patients (post-QC)",
        fontsize=10,
    )
    ax.set_xlabel("Treatment")
    ax.set_ylabel("Patient (total cells)")
    ax.tick_params(axis="x", rotation=0, labelsize=9)
    ax.tick_params(axis="y", rotation=0, labelsize=8)
    fig.tight_layout()

    out = FIG_DIR / f"heatmap_cells_per_cohort_treatment_{C.compartment(value)}.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    for value in ("tum", "epi"):
        table = cohort_x_treatment(value)
        print(f"== {value} ({C.prefix(value)}_raw.h5ad)")
        print(table.assign(Total=table.sum(axis=1)).to_string(), "\n")
        print(f"-> {plot(table, value)}\n")


if __name__ == "__main__":
    main()
