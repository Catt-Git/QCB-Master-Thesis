#!/usr/bin/env python3
"""05_11: cells per patient in each treatment arm, one figure per arm.

Reads the three objects split_epicnv_by_treatment.py wrote and draws, for each arm, one bar
per patient (`cohort`) stacked by `cnv_status`. Patients are sorted by total count; the y axis
is shared across the three figures so the arms can be compared by eye.

Output: figures/05_11_treatment_split/epicnv_counts/cells_per_patient_epicnv_{BASE,PD1,RTPD1}.png

Usage:
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
    python plot_epicnv_counts_by_treatment.py
"""
from __future__ import annotations

import os

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

DATA_DIR = os.environ.get("DATA_DIR", os.path.expanduser("~/Desktop/QCB-Master-Thesis/datasets"))
IN_TEMPLATE = os.path.join(DATA_DIR, "05_tum", "shiao_epicnv_raw_{}.h5ad")
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PHASE_DIR = os.path.dirname(SCRIPT_DIR)
FIG_DIR = os.path.join(PHASE_DIR, "figures", "05_11_treatment_split", "epicnv_counts")
TREATMENTS = ["BASE", "PD1", "RTPD1"]

# Bottom to top. Malignant red is the `tum` colour of 05_2; non_malignant the `epi` blue.
STATUS_ORDER = ["malignant", "borderline", "non_malignant"]
STATUS_COLORS = {"malignant": "#B03A2E", "borderline": "#E0A458", "non_malignant": "#4A6FA5"}


def counts_table(treatment: str) -> pd.DataFrame:
    obs = ad.read_h5ad(IN_TEMPLATE.format(treatment), backed="r").obs
    tab = pd.crosstab(obs["cohort"].astype(str), obs["cnv_status"].astype(str))
    tab = tab.reindex(columns=STATUS_ORDER, fill_value=0)
    return tab.loc[tab.sum(axis=1).sort_values(ascending=False).index]


def plot_arm(treatment: str, tab: pd.DataFrame, ymax: float) -> str:
    fig, ax = plt.subplots(figsize=(max(6.0, 0.38 * len(tab) + 1.5), 4.2))
    x = range(len(tab))
    bottom = pd.Series(0, index=tab.index)
    for status in STATUS_ORDER:
        ax.bar(x, tab[status], bottom=bottom, width=0.75, color=STATUS_COLORS[status],
               edgecolor="white", linewidth=1, label=status.replace("_", "-"))
        bottom += tab[status]
    for i, total in enumerate(bottom):
        ax.text(i, total + ymax * 0.01, f"{total:,}", ha="center", va="bottom",
                fontsize=6.5, color="#444444", rotation=90)

    ax.set_xticks(list(x))
    ax.set_xticklabels([c.replace("Patient", "P") for c in tab.index], rotation=90, fontsize=8)
    ax.set_ylim(0, ymax * 1.12)
    ax.set_ylabel("cells")
    ax.set_title(f"{treatment}: epithelial cells per patient "
                 f"(n = {int(bottom.sum()):,}, {len(tab)} patients)", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#E5E5E5", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(title="cnv_status", frameon=False, fontsize=8, title_fontsize=8, loc="upper right")
    fig.tight_layout()

    out = os.path.join(FIG_DIR, f"cells_per_patient_epicnv_{treatment}.png")
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return out


def main() -> None:
    os.makedirs(FIG_DIR, exist_ok=True)
    tables = {t: counts_table(t) for t in TREATMENTS}
    ymax = max(tab.sum(axis=1).max() for tab in tables.values())
    for treatment, tab in tables.items():
        print(plot_arm(treatment, tab, ymax))


if __name__ == "__main__":
    main()
