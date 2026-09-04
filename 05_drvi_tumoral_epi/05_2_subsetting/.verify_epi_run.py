# AUTOGENERATO per verifica - cancellato subito dopo
import matplotlib
matplotlib.use('Agg')

# ================= cella 2 =================
print('>>> cella 2', flush=True)
# Core scverse libraries
from __future__ import annotations

# Main
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import os

# Plotting
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sns

sc.set_figure_params(dpi=300, facecolor='white')
sc.settings.verbosity = 3

# The one line that makes this notebook the control-set twin of visualization_tum.ipynb.
# setdefault, not assignment: an exported CELL_SET still wins, so this file can be pointed at
# the primary line too if you ever want a second executed copy of it.
os.environ.setdefault('CELL_SET', 'epi')
import importlib
import cell_set as C
# cell_set.py is a sibling file that gets edited while the notebook is open, and Python caches
# modules in sys.modules: without the reload an edit stays invisible until a kernel restart.
importlib.reload(C)

try:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:
    SCRIPT_DIR = os.getcwd()  # interactive fallback
PHASE_DIR = os.path.dirname(SCRIPT_DIR)  # 05_drvi_tumoral_epi/
REPO_ROOT = os.path.dirname(PHASE_DIR)

FIG_DIR = os.environ.get('FIG_DIR', os.path.join(PHASE_DIR, 'figures', '05_2_visualization'))
os.makedirs(FIG_DIR, exist_ok=True)
sc.settings.figdir = FIG_DIR

TABLE_DIR = os.environ.get('TABLE_DIR', os.path.join(PHASE_DIR, 'tables', '05_2_visualization'))
os.makedirs(TABLE_DIR, exist_ok=True)

# Same palette as subset_and_qc.ipynb, so the two notebooks of this step look like one document.
SET_COLORS = {'tum': '#B03A2E', 'epi': '#4A6FA5'}
SET_TITLES = {'tum': 'malignant', 'epi': 'epithelial (post-CNV)'}
QC_COLOR = SET_COLORS[C.cell_set()]
SET_TITLE = SET_TITLES[C.cell_set()]
# Figure filenames carry the compartment ('tum' / 'epicnv'), so the two cell sets never
# overwrite each other's output in a shared FIG_DIR.
SUFFIX = C.compartment()

from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = 'all'

print('FIG_DIR  :', FIG_DIR)
print('TABLE_DIR:', TABLE_DIR)
print('colour   :', QC_COLOR, '| title word:', SET_TITLE, '| file suffix:', SUFFIX)

# ================= cella 3 =================
print('>>> cella 3', flush=True)
DATA_DIR = os.environ.get('DATA_DIR', os.path.join(REPO_ROOT, 'datasets'))
os.environ['DATA_DIR'] = os.path.abspath(os.path.expanduser(DATA_DIR))

C.banner('05_2 visualization')

# Input: the definitive object of this cell set, after clustering_tum.py
INPUT_H5AD = str(C.path('.h5ad'))
print('Input :', INPUT_H5AD, '| exists:', os.path.exists(INPUT_H5AD))

BATCH_KEY = C.BATCH_KEY            # 'cohort'
LABEL_KEY = C.LABEL_KEY            # 'cell_type', post-CNV
PRIOR_KEY = C.PRIOR_LABEL_KEY      # 'cell_type_01_4', pre-CNV
STATUS_KEY = C.STATUS_KEY          # 'cnv_status'
CLUSTER_KEY = 'optscib_tum_leiden' # the same key for both sets, as written by clustering_tum.py

# ================= cella 4 =================
print('>>> cella 4', flush=True)
adata = ad.read_h5ad(INPUT_H5AD)
adata

# ================= cella 5 =================
print('>>> cella 5', flush=True)
# The single decision that changes how the rest of this notebook reads: which label is
# informative on THIS cell set. Same rule as reduce_data_tum.py, so the panels here and the
# diagnostic ones written by the script are colored by the same column.
LABEL_FOR_PLOTS = LABEL_KEY if adata.obs[LABEL_KEY].nunique() > 1 else PRIOR_KEY
LABEL_IS_PRIOR = LABEL_FOR_PLOTS == PRIOR_KEY
# What to put in a figure title, so a reader never has to check which column was used.
LABEL_TITLE = ('pre-CNV CellTypist label' if LABEL_IS_PRIOR else 'post-CNV cell type')
HAS_STATUS = STATUS_KEY in adata.obs and adata.obs[STATUS_KEY].nunique() > 1

print(f'{adata.n_obs:,} cells x {adata.n_vars:,} genes')
print(f'{adata.obs[BATCH_KEY].nunique()} cohorts')
print(f'{LABEL_KEY} (post-CNV): {adata.obs[LABEL_KEY].nunique()} level(s) '
      f'-> {sorted(adata.obs[LABEL_KEY].unique())[:4]}')
print(f'{PRIOR_KEY} (pre-CNV) : {adata.obs[PRIOR_KEY].nunique()} level(s)')
print(f'panels are coloured by {LABEL_FOR_PLOTS!r} ({LABEL_TITLE})')
print(f'cnv_status has both classes: {HAS_STATUS}')
print(f'HVGs: {int(adata.var["highly_variable"].sum()):,}')
print(f'leiden clusters at the selected resolution: {adata.obs[CLUSTER_KEY].nunique()}')

# ================= cella 7 =================
print('>>> cella 7', flush=True)
# What the pre-CNV annotation had called these cells. Under `tum` this is the most informative
# panel of the phase: the CellTypist model has no malignant class, so it had to put every one of
# these aneuploid cells somewhere, and this is where it put them - i.e. this is the error phase 04
# was carrying. Under `epi` the same bars split malignant from non-malignant inside each label,
# which is the same statement seen from the other side.
group_cols = [PRIOR_KEY, STATUS_KEY] if HAS_STATUS else [PRIOR_KEY]
comp = adata.obs.groupby(group_cols, observed=True).size()
comp = comp.unstack(fill_value=0) if HAS_STATUS else comp.to_frame('cells')
comp = comp.loc[comp.sum(axis=1).sort_values(ascending=False).index]

fig, ax = plt.subplots(figsize=(8, max(3, 0.38 * len(comp))))
left = np.zeros(len(comp))
palette = {'malignant': QC_COLOR, 'non_malignant': '#9AA5B1', 'cells': QC_COLOR}
for col in comp.columns:
    ax.barh(range(len(comp)), comp[col], left=left, color=palette.get(col, '#cccccc'), label=str(col))
    left += comp[col].to_numpy()
for yi, tot in enumerate(comp.sum(axis=1)):
    ax.annotate(f'{tot:,}', (tot, yi), va='center', ha='left', fontsize=8,
                xytext=(4, 0), textcoords='offset points')
ax.set_yticks(range(len(comp)))
ax.set_yticklabels(comp.index, fontsize=9)
ax.invert_yaxis()
ax.margins(x=0.12)
ax.set_xlabel('cells')
ax.set_title(f'What the pre-CNV annotation had called these cells ({SET_TITLE})')
if comp.shape[1] > 1:
    ax.legend(frameon=False, fontsize=8)
ax.grid(axis='x', linestyle='-', alpha=0.3)
ax.set_axisbelow(True)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, f'barplot_composition_by_prior_label_{SUFFIX}_unintegrated.png'),
            dpi=300, bbox_inches='tight')
plt.show()
comp

# ================= cella 8 =================
print('>>> cella 8', flush=True)
# The plane the call was actually taken in. call_malignant.ipynb thresholded `cnv_score` and
# `cnv_corr` against a per-cohort null built on the stromal cells, so the cut is NOT a single
# global line and this scatter should not be read as if it were: cohorts have different cuts.
# What it does show is whether the two axes agree, and how far from the origin the called cells sit.
rng = np.random.default_rng(0)
idx = np.sort(rng.choice(adata.n_obs, size=min(30000, adata.n_obs), replace=False))
sub_obs = adata.obs.iloc[idx]

fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

if HAS_STATUS:
    for status, color in (('non_malignant', '#9AA5B1'), ('malignant', QC_COLOR)):
        m = (sub_obs[STATUS_KEY] == status).to_numpy()
        axes[0].scatter(sub_obs['cnv_score'][m], sub_obs['cnv_corr'][m], s=2, alpha=0.25,
                        linewidth=0, color=color, label=f'{status} ({m.sum():,})')
    axes[0].legend(markerscale=6, fontsize=8, frameon=False)
else:
    axes[0].scatter(sub_obs['cnv_score'], sub_obs['cnv_corr'], s=2, alpha=0.25,
                    linewidth=0, color=QC_COLOR)
axes[0].set_xlabel('cnv_score')
axes[0].set_ylabel('cnv_corr')
axes[0].set_title('The (score, corr) plane\n(per-cohort cuts, not a global line)', fontsize=10)

for ax, col in zip(axes[1:], ('cnv_score', 'cnv_corr')):
    if HAS_STATUS:
        for status, color in (('non_malignant', '#9AA5B1'), ('malignant', QC_COLOR)):
            vals = adata.obs.loc[adata.obs[STATUS_KEY] == status, col]
            sns.kdeplot(x=vals, ax=ax, fill=True, alpha=0.4, color=color, label=status, cut=0)
        ax.legend(fontsize=8, frameon=False)
    else:
        sns.histplot(adata.obs[col], bins=80, ax=ax, color=QC_COLOR)
    ax.set_title(f'{col} ({SET_TITLE})', fontsize=10)
    ax.set_xlabel(col)

fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, f'scatter_cnv_score_vs_corr_{SUFFIX}_unintegrated.png'),
            dpi=300, bbox_inches='tight')
plt.show()

adata.obs[['cnv_score', 'cnv_corr']].describe().round(3)

# ================= cella 9 =================
print('>>> cella 9', flush=True)
# The malignant fraction per cohort. Under `epi` this is a real quantity - how much of each
# patient's epithelium is aneuploid - and its spread is what MIN_CELLS_PER_COHORT=200 interacts
# with. Under `tum` it is 100% everywhere by construction, so only the cell counts are drawn.
cells_per_cohort = adata.obs[BATCH_KEY].value_counts()
cells_per_cohort = cells_per_cohort[cells_per_cohort > 0].sort_values(ascending=False)

if HAS_STATUS:
    frac = (pd.crosstab(adata.obs[BATCH_KEY], adata.obs[STATUS_KEY], normalize='index') * 100)
    frac = frac.reindex(cells_per_cohort.index)

    fig, axes = plt.subplots(2, 1, figsize=(14, 9), sharex=True)
    cells_per_cohort.plot(kind='bar', ax=axes[0], width=0.85, color=QC_COLOR)
    axes[0].set_ylabel('cells')
    axes[0].set_title(f'Cells per cohort and malignant fraction ({SET_TITLE})')

    frac[['malignant']].plot(kind='bar', ax=axes[1], width=0.85, color=QC_COLOR, legend=False)
    axes[1].axhline(frac['malignant'].median(), color='k', ls='--', lw=1,
                    label=f"median {frac['malignant'].median():.1f}%")
    axes[1].set_ylabel('% malignant')
    axes[1].set_xlabel('Cohort')
    axes[1].legend(frameon=False, fontsize=8)
    for ax in axes:
        ax.grid(axis='y', linestyle='-', alpha=0.4)
        ax.set_axisbelow(True)
    axes[1].tick_params(axis='x', rotation=45, labelsize=9)
    axes[1].set_xticklabels(axes[1].get_xticklabels(), ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, f'barplot_malignant_fraction_per_cohort_{SUFFIX}_unintegrated.png'),
                dpi=300, bbox_inches='tight')
    plt.show()
    frac.round(1)
else:
    print(f'cnv_status is constant ({adata.obs[STATUS_KEY].iloc[0]!r}): '
          'the malignant fraction is 100% in every cohort by construction, nothing to plot here. '
          'The per-cohort counts are in the barplot section below.')

# ================= cella 10 =================
print('>>> cella 10', flush=True)
# How far the 05_1 re-annotation moved the labels, on the cells that HAVE both. Only meaningful
# under `epi`, where the non-malignant cells carry a re-run CellTypist label: under `tum` every
# cell is `malignant` post-CNV and the crosstab is a single column.
#
# Read this as CHURN, not as a measurement of contamination. Most of it is CellTypist's majority
# vote redrawing its own over-clustering (leiden at resolution 30, rebuilt from scratch on the
# reduced population and stochastic), not the tumour cells having stopped voting. On the full
# object 95.3% of the changes stayed inside the same lineage and 15.3% were CD4 <-> CD8 swaps,
# which removing epithelial tumour cells cannot cause. See the docstring of
# 05_1/recelltypist_nonmalignant.py.
if adata.obs[LABEL_KEY].nunique() > 1:
    nm = adata.obs.loc[adata.obs[LABEL_KEY] != C.MALIGNANT_LABEL]
    moved = int((nm[PRIOR_KEY].astype(str) != nm[LABEL_KEY].astype(str)).sum())
    print(f'{moved:,} of {len(nm):,} non-malignant cells carry a different label after the '
          f're-annotation ({moved / max(len(nm), 1):.1%})')

    ct = pd.crosstab(nm[PRIOR_KEY], nm[LABEL_KEY])
    ct = ct.loc[ct.sum(axis=1).sort_values(ascending=False).index]
    ct = ct[ct.sum().sort_values(ascending=False).index]
    ct_pct = ct.div(ct.sum(axis=1).replace(0, np.nan), axis=0) * 100

    fig, ax = plt.subplots(figsize=(1 + 0.55 * ct.shape[1], 1 + 0.45 * ct.shape[0]))
    sns.heatmap(ct_pct, cmap='Blues', vmin=0, vmax=100, annot=True, fmt='.0f',
                annot_kws={'fontsize': 7}, cbar_kws={'label': '% of the pre-CNV label'}, ax=ax)
    ax.set_xlabel('post-CNV label (05_1 re-annotation)')
    ax.set_ylabel('pre-CNV label (01_4)')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.set_title('Where the re-annotation moved the non-malignant cells\n'
                 '(row-normalized; churn, not a contamination measure)', fontsize=10)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, f'heatmap_reannotation_churn_{SUFFIX}_unintegrated.png'),
                dpi=300, bbox_inches='tight')
    plt.show()

    ct.to_csv(os.path.join(TABLE_DIR, f'reannotation_prior_vs_post_{SUFFIX}.csv'))
    print('wrote', os.path.join(TABLE_DIR, f'reannotation_prior_vs_post_{SUFFIX}.csv'))
    ct
else:
    print(f'{LABEL_KEY} is constant ({adata.obs[LABEL_KEY].iloc[0]!r}) on this cell set: '
          'there is no re-annotated label to compare against. Run CELL_SET=epi for this panel.')

# ================= cella 14 =================
print('>>> cella 14', flush=True)
# Cells per cohort, stacked by treatment. `compartment` (and `fraction` with it) is constant
# here by construction, so treatment takes its place as the stacking variable, exactly as in 04_1.
TREATMENT_COLORS = {'BASE': 'lightskyblue', 'PD1': 'pink', 'RTPD1': 'coral'}
treatment_order = [t for t in TREATMENT_COLORS if t in set(adata.obs['treatment'])]

cohort_treatment = pd.crosstab(adata.obs[BATCH_KEY], adata.obs['treatment'])[treatment_order]
cohort_treatment = cohort_treatment.loc[cohort_treatment.sum(axis=1) > 0]
cohort_treatment = cohort_treatment.loc[cohort_treatment.sum(axis=1).sort_values(ascending=False).index]

ax = cohort_treatment.plot(kind='bar', stacked=True, figsize=(14, 6), width=0.85,
                           color=[TREATMENT_COLORS[t] for t in treatment_order])
ax.set_title(f'Number of {SET_TITLE} cells per cohort')
ax.set_xlabel('Cohort')
ax.set_ylabel('Number of Cells')
ax.tick_params(axis='x', rotation=45, labelsize=9)
ax.set_xticklabels(ax.get_xticklabels(), ha='right')
ax.legend(title='Treatment', bbox_to_anchor=(1.02, 1), loc='upper left')
ax.grid(axis='y', linestyle='-', alpha=0.4)
ax.set_axisbelow(True)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, f'histogram_cells_per_cohort_{SUFFIX}_unintegrated.png'),
            dpi=300, bbox_inches='tight')
plt.show()

# ================= cella 15 =================
print('>>> cella 15', flush=True)
THRESHOLD = 1000  # print the exact count below this

# On a single linear axis the smallest cohorts are a couple of pixels tall and unreadable off the
# y-ticks: print the exact count for anything below THRESHOLD above its bar. The MIN_CELLS_PER_COHORT
# line is drawn so it is visible that no cohort under it survived.
fig, ax = plt.subplots(figsize=(14, 6))
cells_per_cohort.plot(kind='bar', width=0.85, color=QC_COLOR, ax=ax, legend=False)
ax.axhline(C.MIN_CELLS_PER_COHORT, color='k', ls='--', lw=1,
           label=f'MIN_CELLS_PER_COHORT = {C.MIN_CELLS_PER_COHORT}')
ax.grid(axis='y', linestyle='-', alpha=0.4)
ax.set_axisbelow(True)

for patch, value in zip(ax.patches, cells_per_cohort.values):
    if value < THRESHOLD:
        ax.annotate(f'{int(value)}', (patch.get_x() + patch.get_width() / 2, patch.get_height()),
                    ha='center', va='bottom', fontsize=8, xytext=(0, 2), textcoords='offset points')

ax.set_title(f'Number of {SET_TITLE} cells per cohort')
ax.set_xlabel('Cohort')
ax.set_ylabel('Number of Cells')
ax.legend(frameon=False, fontsize=9)
ax.tick_params(axis='x', rotation=45, labelsize=9)
ax.set_xticklabels(ax.get_xticklabels(), ha='right')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, f'histogram_cells_per_cohort_single_axis_{SUFFIX}_unintegrated.png'),
            dpi=300, bbox_inches='tight')
plt.show()

# ================= cella 16 =================
print('>>> cella 16', flush=True)
# One palette for the label used by every composition panel below, so a label keeps its colour
# across the notebook. `malignant`, when it is present (CELL_SET=epi), is forced to the phase
# colour instead of taking whatever tab20 slot its rank lands on: it is the class the figure is
# about, and it has to be findable at a glance.
label_order = adata.obs[LABEL_FOR_PLOTS].value_counts()
label_order = label_order[label_order > 0].index

_tab20 = list(plt.get_cmap('tab20').colors) + list(plt.get_cmap('tab20b').colors)
_others = [c for c in label_order if c != C.MALIGNANT_LABEL]
LABEL_COLORS = {ct: _tab20[i % len(_tab20)] for i, ct in enumerate(_others)}
if C.MALIGNANT_LABEL in set(label_order):
    LABEL_COLORS[C.MALIGNANT_LABEL] = QC_COLOR

# Cells per label, most abundant on top, with the share of the compartment annotated.
# Log x-axis: the labels span three orders of magnitude (Lumsec-prol vs Lumsec-HLA).
counts = adata.obs[LABEL_FOR_PLOTS].value_counts().loc[label_order]
total = adata.n_obs
y = np.arange(len(label_order))

fig, ax = plt.subplots(figsize=(10, 0.32 * len(label_order) + 2))
ax.barh(y, counts.values, color=[LABEL_COLORS[c] for c in label_order])
ax.set_yticks(y)
ax.set_yticklabels(label_order, fontsize=9)
ax.invert_yaxis()
for yi, val in zip(y, counts.values):
    ax.annotate(f'{val:,} ({100 * val / total:.1f}%)', (val, yi), va='center', ha='left',
                fontsize=8, xytext=(4, 0), textcoords='offset points')
ax.set_title(f'Number of cells per {LABEL_TITLE} ({SET_TITLE})')
ax.set_xlabel('Number of Cells')
ax.set_ylabel(LABEL_FOR_PLOTS)
ax.set_xscale('log')
ax.set_xlim(left=max(1, counts.min() / 2), right=counts.max() * 4)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, f'histogram_cells_per_label_{SUFFIX}_unintegrated.png'),
            dpi=300, bbox_inches='tight')
plt.show()

# ================= cella 17 =================
print('>>> cella 17', flush=True)
# Cells per cohort and label, every level with its own block so nothing is pooled into "other"
# and the stacks are directly comparable across cohorts.
cohort_label = pd.crosstab(adata.obs[BATCH_KEY], adata.obs[LABEL_FOR_PLOTS])
cohort_label = cohort_label.loc[cohort_label.sum(axis=1) > 0]
cohort_label = cohort_label.loc[cohort_label.sum(axis=1).sort_values(ascending=False).index]
cohort_label = cohort_label[[c for c in label_order if c in cohort_label.columns]]

N_PANELS = 2 if len(cohort_label) <= 20 else 3
panels = np.array_split(np.arange(len(cohort_label)), N_PANELS)

fig, axes = plt.subplots(N_PANELS, 1, figsize=(13, 3.6 * N_PANELS))
for ax, idx in zip(np.atleast_1d(axes), panels):
    block = cohort_label.iloc[idx]
    block.plot(kind='bar', stacked=True, ax=ax, width=0.85, legend=False,
               color=[LABEL_COLORS[c] for c in block.columns])
    for x, tot in enumerate(block.sum(axis=1)):
        ax.annotate(f'{tot:,}', (x, tot), ha='center', va='bottom', fontsize=7,
                    xytext=(0, 2), textcoords='offset points')
    ax.margins(y=0.12)
    ax.set_xlabel('')
    ax.set_ylabel('Number of Cells')
    ax.tick_params(axis='x', rotation=45, labelsize=9)
    ax.set_xticklabels(ax.get_xticklabels(), ha='right')
    ax.grid(axis='y', linestyle='-', alpha=0.4)
    ax.set_axisbelow(True)

np.atleast_1d(axes)[0].set_title(f'Cells per cohort and {LABEL_TITLE} ({SET_TITLE}), all levels')
np.atleast_1d(axes)[-1].set_xlabel('Cohort')
handles, legend_labels = np.atleast_1d(axes)[0].get_legend_handles_labels()
fig.legend(handles, legend_labels, title=LABEL_FOR_PLOTS, loc='center left',
           bbox_to_anchor=(0.99, 0.5), fontsize=7, ncol=1 + len(label_order) // 30)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, f'barplot_label_per_cohort_{SUFFIX}_unintegrated.png'),
            dpi=300, bbox_inches='tight')
plt.show()

# The same table on disk, in the same order as the figure, with row and column totals.
cohort_label_csv = cohort_label.copy()
cohort_label_csv.loc['total'] = cohort_label_csv.sum()
cohort_label_csv['total'] = cohort_label_csv.sum(axis=1)
out = os.path.join(TABLE_DIR, f'cells_per_cohort_and_label_{SUFFIX}.csv')
cohort_label_csv.to_csv(out)
print('wrote', out)
cohort_label

# ================= cella 18 =================
print('>>> cella 18', flush=True)
# Composition per treatment, in percent: the first look at whether the compartment shifts under
# PD1 / RT+PD1. Descriptive only - cells are not independent within a patient, and with 19 cohorts
# a per-patient model is what would be needed to say anything inferential.
comp_tr = pd.crosstab(adata.obs['treatment'], adata.obs[LABEL_FOR_PLOTS], normalize='index') * 100
comp_tr = comp_tr.reindex(index=treatment_order, columns=label_order)
comp_tr.plot(kind='bar', stacked=True, figsize=(10, 6), width=0.7,
             color=[LABEL_COLORS[c] for c in comp_tr.columns])
plt.title(f'Composition per treatment by {LABEL_TITLE} ({SET_TITLE})')
plt.xlabel('Treatment')
plt.ylabel('Percentage of Cells')
plt.xticks(rotation=0)
plt.legend(title=LABEL_FOR_PLOTS, bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, f'barplot_composition_per_treatment_{SUFFIX}_unintegrated.png'),
            dpi=300, bbox_inches='tight')
plt.show()
comp_tr.round(1)

# ================= cella 19 =================
print('>>> cella 19', flush=True)
# The same composition per response group (NR / R1 / R2). Descriptive only for the same reason,
# and one step weaker still: response is a patient-level label, so the effective n is the number
# of cohorts per group, not the number of cells.
RESPONSE_ORDER = ['NR', 'R1', 'R2']
response_order = [r for r in RESPONSE_ORDER if r in set(adata.obs['response'])]

comp_resp = pd.crosstab(adata.obs['response'], adata.obs[LABEL_FOR_PLOTS], normalize='index') * 100
comp_resp = comp_resp.reindex(index=response_order, columns=label_order)
comp_resp.plot(kind='bar', stacked=True, figsize=(10, 6), width=0.7,
               color=[LABEL_COLORS[c] for c in comp_resp.columns])
plt.title(f'Composition per response by {LABEL_TITLE} ({SET_TITLE})')
plt.xlabel('Response')
plt.ylabel('Percentage of Cells')
plt.xticks(rotation=0)
plt.legend(title=LABEL_FOR_PLOTS, bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, f'barplot_composition_per_response_{SUFFIX}_unintegrated.png'),
            dpi=300, bbox_inches='tight')
plt.show()

print('cohorts per response group (the actual n):')
adata.obs.groupby('response', observed=True)[BATCH_KEY].nunique().to_frame('cohorts')
comp_resp.round(1)

# ================= cella 21 =================
print('>>> cella 21', flush=True)
PHASE_ORDER = ['G1', 'S', 'G2M']
PHASE_COLORS = {'G1': '#377eb8', 'S': '#4daf4a', 'G2M': '#e41a1c'}

phase_counts = adata.obs['phase'].value_counts().reindex(PHASE_ORDER)
print(phase_counts.to_string())
print(f'\ncycling (S + G2M): {phase_counts[["S", "G2M"]].sum() / adata.n_obs:.1%}')

tmp = pd.crosstab(adata.obs[BATCH_KEY], adata.obs['phase'], normalize='index')
tmp = tmp.reindex(index=cells_per_cohort.index, columns=PHASE_ORDER) * 100
tmp.plot.bar(stacked=True, figsize=(13, 5), color=[PHASE_COLORS[p] for p in PHASE_ORDER])
plt.title(f'Cell cycle phase distribution by cohort ({SET_TITLE})')
plt.xlabel('Cohort')
plt.ylabel('Percentage of Cells')
plt.xticks(rotation=45, ha='right', fontsize=9)
plt.legend(title='Phase', bbox_to_anchor=(1.02, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, f'barplot_cc_vs_cohort_{SUFFIX}_unintegrated.png'),
            dpi=300, bbox_inches='tight')
plt.show()

# ================= cella 22 =================
print('>>> cella 22', flush=True)
# By label. Under `tum` this is a coherence check worth reading: `Lumsec-prol` means
# *proliferating* in the normal-breast atlas, so if the borrowed label carries any real signal at
# all its G2M share has to be the highest here. Two independent methods agreeing on one group.
tmp = (pd.crosstab(adata.obs[LABEL_FOR_PLOTS], adata.obs['phase'], normalize='index')
       .reindex(index=label_order, columns=PHASE_ORDER) * 100)
tmp.plot.bar(stacked=True, figsize=(12, 6), color=[PHASE_COLORS[p] for p in PHASE_ORDER])
plt.title(f'Cell cycle phase distribution by {LABEL_TITLE} ({SET_TITLE})')
plt.xlabel(LABEL_FOR_PLOTS)
plt.ylabel('Percentage of Cells')
plt.xticks(rotation=90, fontsize=9)
plt.legend(title='Phase', bbox_to_anchor=(1.02, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, f'barplot_cc_vs_label_{SUFFIX}_unintegrated.png'),
            dpi=300, bbox_inches='tight')
plt.show()
tmp.round(1)

# ================= cella 23 =================
print('>>> cella 23', flush=True)
# The S / G2M scores themselves, which is what `phase` is a hard call off. Worth seeing as a
# continuum: the phase assignment cuts a cloud, and how clean that cut is decides how much of
# 05_9's answer is method rather than biology.
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

for ph in PHASE_ORDER:
    m = (adata.obs['phase'] == ph).to_numpy()
    axes[0].scatter(adata.obs['S_score'][m], adata.obs['G2M_score'][m], s=2, alpha=0.3,
                    linewidth=0, color=PHASE_COLORS[ph], label=f'{ph} ({m.sum():,})')
axes[0].axhline(0, color='k', lw=0.6, ls='--')
axes[0].axvline(0, color='k', lw=0.6, ls='--')
axes[0].set_xlabel('S_score')
axes[0].set_ylabel('G2M_score')
axes[0].set_title(f'Cell cycle scores ({SET_TITLE})')
axes[0].legend(markerscale=6, fontsize=8, frameon=False)

sns.violinplot(data=adata.obs, x='phase', y='S_score', order=PHASE_ORDER, ax=axes[1],
               hue='phase', hue_order=PHASE_ORDER, palette=PHASE_COLORS, legend=False, cut=0)
axes[1].set_title('S_score by assigned phase')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, f'scatter_cc_scores_{SUFFIX}_unintegrated.png'),
            dpi=300, bbox_inches='tight')
plt.show()

# ================= cella 25 =================
print('>>> cella 25', flush=True)
fig, axes = plt.subplots(1, 2, figsize=(10, 5))

sns.histplot(adata.obs['total_counts'], bins=100, kde=False, log_scale=True, ax=axes[0],
             color=QC_COLOR)
axes[0].set_title('Library size (raw counts)')
axes[0].set_xlabel('Total raw counts per cell (log scale)')

sns.histplot(np.asarray(adata.X.sum(axis=1)).ravel(), bins=100, kde=False, ax=axes[1],
             color=QC_COLOR)
axes[1].set_title('Sum of log-normalized expression')
axes[1].set_xlabel('Sum of scran log-normalized expression per cell')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, f'histogram_total_counts_raw_vs_normalized_{SUFFIX}_unintegrated.png'),
            dpi=300, bbox_inches='tight')
plt.show()

# ================= cella 26 =================
print('>>> cella 26', flush=True)
# The same comparison as violins, on a random subsample for readability. Same caveat: the two
# panels are different units, so their widths are not comparable. Indices are drawn once and
# reused, so nothing copies the expression matrix.
rng = np.random.default_rng(0)
idx_norm = np.sort(rng.choice(adata.n_obs, size=min(20000, adata.n_obs), replace=False))

raw_log = np.log1p(np.asarray(adata.layers['counts'][idx_norm].sum(axis=1)).ravel())
norm_sum = np.asarray(adata.X[idx_norm].sum(axis=1)).ravel()

fig, ax = plt.subplots(1, 2, figsize=(12, 5))
sns.violinplot(y=raw_log, ax=ax[0], color=QC_COLOR, cut=0)
ax[0].set_title('Library size (raw counts, log1p)')
ax[0].set_ylabel('log1p(total raw counts per cell)')

sns.violinplot(y=norm_sum, ax=ax[1], color=QC_COLOR, cut=0)
ax[1].set_title('Sum of log-normalized expression')
ax[1].set_ylabel('Sum of scran log-normalized expression per cell')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, f'violin_normalization_raw_vs_scran_{SUFFIX}_unintegrated.png'),
            dpi=300, bbox_inches='tight')
plt.show()

# ================= cella 28 =================
print('>>> cella 28', flush=True)
n_bad = int((adata.obs['size_factors'] <= 0).sum())
print(f'size_factors: min={adata.obs["size_factors"].min():.4f}, '
      f'max={adata.obs["size_factors"].max():.4f}, n(<=0)={n_bad}')
assert n_bad == 0, 'a non-positive size factor survived: the normalization is not usable'

sc.pl.scatter(adata, 'size_factors', 'total_counts', color=LABEL_FOR_PLOTS,
              save=f'_size_factors_vs_total_counts_by_label_{SUFFIX}_unintegrated.png')
sc.pl.scatter(adata, 'size_factors', 'n_genes_by_counts', color=LABEL_FOR_PLOTS,
              save=f'_size_factors_vs_n_genes_by_counts_by_label_{SUFFIX}_unintegrated.png')

sns.displot(adata.obs['size_factors'], bins=50, kde=False, color=QC_COLOR)
plt.savefig(os.path.join(FIG_DIR, f'displot_size_factors_{SUFFIX}_unintegrated.png'),
            dpi=300, bbox_inches='tight')
plt.show()

adata.obs['size_factors'].describe()

# ================= cella 30 =================
print('>>> cella 30', flush=True)
hvg_genes = adata.var_names[adata.var['highly_variable']]
print(f'HVGs: {len(hvg_genes):,}')

HVG_REFS = [
    ('shiao_hvg_2k_unintegrated_list.csv', '01_5 full-object'),
    (os.path.join('03_nonimm', 'shiao_nonimm_hvg_2k_list.csv'), '03_1 non-immune'),
    (os.path.join('04_epi', 'shiao_epi_hvg_2k_list.csv'), '04_1 epithelial'),
]

rows = []
for rel, label in HVG_REFS:
    p = os.path.join(os.environ['DATA_DIR'], rel)
    if not os.path.exists(p):
        print(f'{label:20s}: list not found at {p}, skipped')
        continue
    other = set(pd.read_csv(p, header=None)[0].astype(str))
    shared = set(hvg_genes) & other
    rows.append({'reference': label, 'shared': len(shared),
                 'pct': round(100 * len(shared) / len(hvg_genes), 1)})
    print(f'{label:20s}: {len(shared):,}/{len(hvg_genes):,} ({100 * len(shared) / len(hvg_genes):.1f}%)')

# The genes this cell set selected and NO earlier phase did: the ones only visible once the
# population was narrowed this far.
all_earlier = set()
for rel, _ in HVG_REFS:
    p = os.path.join(os.environ['DATA_DIR'], rel)
    if os.path.exists(p):
        all_earlier |= set(pd.read_csv(p, header=None)[0].astype(str))
only_here = sorted(set(hvg_genes) - all_earlier)
print(f'\nnew in this selection (in none of the above): {len(only_here):,}')
print(f'first 30: {only_here[:30]}')

if rows:
    overlap = pd.DataFrame(rows)
    overlap.to_csv(os.path.join(TABLE_DIR, f'hvg_overlap_{SUFFIX}.csv'), index=False)
    overlap

# ================= cella 31 =================
print('>>> cella 31', flush=True)
# Where the HVGs sit in the mean/variance plane. scib's batch-aware hvg_batch does not keep the
# per-gene statistics in .var, so they are recomputed here from .X - cheap, and it turns a skipped
# panel into a real one.
if {'means', 'dispersions_norm'} <= set(adata.var.columns):
    means = adata.var['means'].to_numpy()
    disp = adata.var['dispersions_norm'].to_numpy()
    ylab = 'Normalized dispersion'
else:
    X = adata.X
    means = np.asarray(X.mean(axis=0)).ravel()
    sq = np.asarray(X.multiply(X).mean(axis=0)).ravel() if hasattr(X, 'multiply') \
        else np.asarray(np.square(X).mean(axis=0)).ravel()
    disp = np.sqrt(np.maximum(sq - means ** 2, 0))
    ylab = 'Standard deviation (log-normalized)'

hv = adata.var['highly_variable'].to_numpy()
fig, ax = plt.subplots(figsize=(6.5, 5))
ax.scatter(means[~hv], disp[~hv], s=2, alpha=0.2, linewidth=0, color='grey', label='other genes')
ax.scatter(means[hv], disp[hv], s=2, alpha=0.5, linewidth=0, color=QC_COLOR, label='highly variable')
ax.set_xscale('symlog', linthresh=1e-3)
ax.set_xlabel('Mean expression')
ax.set_ylabel(ylab)
ax.set_title(f'Batch-aware HVG selection ({SET_TITLE})')
ax.legend(markerscale=4, frameon=False)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, f'scatter_hvg_mean_dispersion_{SUFFIX}_unintegrated.png'),
            dpi=300, bbox_inches='tight')
plt.show()

# ================= cella 33 =================
print('>>> cella 33', flush=True)
UMAP_KEYS = {
    'label': LABEL_FOR_PLOTS,
    'cohort': BATCH_KEY,
    'treatment': 'treatment',
    'response': 'response',
    'phase': 'phase',
    'cnv_score': 'cnv_score',
    'cnv_corr': 'cnv_corr',
    'n_genes_by_counts': 'n_genes_by_counts',
    'total_counts': 'total_counts',
    'mito': 'pct_counts_mt',
    'ribo': 'pct_counts_ribo',
    'size_factors': 'size_factors',
    'leiden': CLUSTER_KEY,
}
if HAS_STATUS:
    UMAP_KEYS['cnv_status'] = STATUS_KEY

UMAP_SEED = 0

# Plot the cells in random order: otherwise the last cohort in the object is drawn on top of every
# other one and the batch panel looks far worse (or better) than it is. A lightweight AnnData with
# only obs + the UMAP coordinates, so nothing copies the expression matrix.
_order = np.random.default_rng(UMAP_SEED).permutation(adata.n_obs)
adata_plot = ad.AnnData(
    obs=adata.obs.iloc[_order].copy(),
    obsm={'X_umap': adata.obsm['X_umap'][_order]},
    uns={k: v for k, v in adata.uns.items() if k.endswith('_colors')},
)

for name, col in UMAP_KEYS.items():
    sc.pl.umap(adata_plot, color=col, save=f'_{name}_{SUFFIX}_unintegrated.png')

# ================= cella 34 =================
print('>>> cella 34', flush=True)
# The metadata panels in one grid, for the thesis figure: annotation + experimental design only,
# the QC/continuous panels stay as individual UMAPs above. No `compartment`/`fraction` panel, both
# are constant by construction. Layout reused from the 04_1 and 01_6 grids so the three can be
# shown side by side.
UMAP_COMBINED_KEYS = [LABEL_FOR_PLOTS, BATCH_KEY, 'treatment', 'response', 'phase']
if HAS_STATUS:
    UMAP_COMBINED_KEYS.insert(1, STATUS_KEY)

with plt.rc_context({'figure.figsize': (7, 7)}):
    sc.pl.umap(adata_plot, color=UMAP_COMBINED_KEYS, ncols=2, wspace=0.8, hspace=0.25,
               save=f'_combined_{SUFFIX}_unintegrated.png')

# ================= cella 35 =================
print('>>> cella 35', flush=True)
# The figure this phase exists to produce, when it has both classes: the whole epithelium with
# the aneuploid cells marked. Everything grey, malignant in the phase colour and drawn last, so
# the tumour is findable at a glance instead of being one entry in a twelve-colour legend.
# Under `tum` there is nothing to contrast against, so the panel is skipped.
if HAS_STATUS:
    xy = adata_plot.obsm['X_umap']
    is_mal = (adata_plot.obs[STATUS_KEY] == C.MALIGNANT_LABEL).to_numpy()

    fig, axes = plt.subplots(1, 2, figsize=(15, 7))

    axes[0].scatter(xy[~is_mal, 0], xy[~is_mal, 1], s=1.5, alpha=0.35, linewidth=0,
                    color='#C8CDD3', label=f'non-malignant ({(~is_mal).sum():,})')
    axes[0].scatter(xy[is_mal, 0], xy[is_mal, 1], s=1.5, alpha=0.55, linewidth=0,
                    color=QC_COLOR, label=f'malignant ({is_mal.sum():,})')
    axes[0].legend(markerscale=8, fontsize=9, frameon=False, loc='best')
    axes[0].set_title('inferCNV call')

    # The same embedding under the post-CNV labels, which is the direct comparison: where the
    # tumour sits against the states CellTypist assigns to what is left.
    for ct in label_order:
        m = (adata_plot.obs[LABEL_FOR_PLOTS] == ct).to_numpy()
        if not m.any():
            continue
        axes[1].scatter(xy[m, 0], xy[m, 1], s=1.5, alpha=0.5, linewidth=0,
                        color=LABEL_COLORS[ct], label=f'{ct} ({m.sum():,})')
    axes[1].legend(markerscale=8, fontsize=7, frameon=False, loc='center left',
                   bbox_to_anchor=(1.01, 0.5))
    axes[1].set_title(f'{LABEL_TITLE}')

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('UMAP1')
        ax.set_ylabel('UMAP2')
    fig.suptitle(f'Epithelial compartment after inferCNV ({adata.n_obs:,} cells, unintegrated)')
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, f'umap_malignant_highlight_{SUFFIX}_unintegrated.png'),
                dpi=300, bbox_inches='tight')
    plt.show()
else:
    print('cnv_status is constant on this cell set: the malignant-vs-normal panel needs '
          'CELL_SET=epi, where both classes are present.')

# ================= cella 37 =================
print('>>> cella 37', flush=True)
# Test only: t-SNE on the same unintegrated PCs the UMAP above was built on, so the two
# embeddings are directly comparable. Computed here and kept in memory -- the .h5ad on disk is not
# touched, and nothing downstream reads `X_tsne`.
# ~36k cells (tum) / ~75k (epi): expect minutes, not seconds.
TSNE_SEED = 0
N_PCS = min(50, adata.obsm['X_pca'].shape[1])

sc.tl.tsne(adata, use_rep='X_pca', n_pcs=N_PCS, perplexity=30, random_state=TSNE_SEED)

_order_tsne = np.random.default_rng(TSNE_SEED).permutation(adata.n_obs)
adata_tsne_plot = ad.AnnData(
    obs=adata.obs.iloc[_order_tsne].copy(),
    obsm={'X_tsne': adata.obsm['X_tsne'][_order_tsne]},
    uns={k: v for k, v in adata.uns.items() if k.endswith('_colors')},
)

TSNE_KEYS = [LABEL_FOR_PLOTS, BATCH_KEY] + ([STATUS_KEY] if HAS_STATUS else [])

with plt.rc_context({'figure.figsize': (7, 7)}):
    sc.pl.tsne(adata_tsne_plot, color=TSNE_KEYS, ncols=2, wspace=0.8, hspace=0.25,
               save=f'_combined_{SUFFIX}_unintegrated.png')

# ================= cella 39 =================
print('>>> cella 39', flush=True)
# The resolution/NMI profile behind the selected clustering.
profile_csv = str(C.path('_leiden_resolution_profile.csv'))
if f'{CLUSTER_KEY}_profile' in adata.uns:
    prof = pd.DataFrame({
        'resolution': adata.uns[f'{CLUSTER_KEY}_profile']['resolution'],
        'score': adata.uns[f'{CLUSTER_KEY}_profile']['nmi'],
    })
else:
    prof = pd.read_csv(profile_csv)

nmi_target = adata.uns.get(f'{CLUSTER_KEY}_label_key', '(not recorded)')
best = prof.loc[prof['score'].idxmax()]
spread = float(prof['score'].max() - prof['score'].min())

fig, ax = plt.subplots(figsize=(7.5, 4.5))
ax.plot(prof['resolution'], prof['score'], marker='o', ms=4, color=QC_COLOR)
ax.axvline(best['resolution'], color='k', ls='--', lw=1)
ax.annotate(f'best: res={best["resolution"]:g}, NMI={best["score"]:.3f}',
            (best['resolution'], best['score']), xytext=(6, -12),
            textcoords='offset points', fontsize=9)
ax.set_xlabel('Leiden resolution')
ax.set_ylabel(f'NMI vs {nmi_target}')
ax.set_title(f'Optimal-resolution sweep ({SET_TITLE})\ntarget: {nmi_target}', fontsize=10)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, f'leiden_resolution_profile_{SUFFIX}_unintegrated.png'),
            dpi=300, bbox_inches='tight')
plt.show()

print(f'NMI target      : {nmi_target}')
print(f'selected        : res={best["resolution"]:g}, NMI={best["score"]:.4f}, '
      f'{adata.obs[CLUSTER_KEY].nunique()} clusters')
print(f'profile spread  : {spread:.4f} (max {prof["score"].max():.4f} - min {prof["score"].min():.4f})')

# clustering_tum.py warns when the optimum is the LARGEST resolution tested. The symmetric case
# matters just as much and is not covered there: an optimum at the smallest grid point means the
# sweep never looked in the direction the criterion was still improving in.
if best['resolution'] == prof['resolution'].min():
    print(f'\n[!] the optimum is the SMALLEST resolution tested ({best["resolution"]:g}): the grid '
          'was never extended below it, so this is a boundary solution, not an interior optimum.')
if best['resolution'] == prof['resolution'].max():
    print(f'\n[!] the optimum is the LARGEST resolution tested ({best["resolution"]:g}): extend the '
          'grid before trusting it.')
if spread < 0.05:
    print(f'\n[!] the profile spans {spread:.4f} NMI across the whole grid '
          f'({100 * spread / prof["score"].max():.1f}% of the maximum): the criterion barely '
          'discriminates between resolutions, so the selected one should be REPORTED as a '
          'near-arbitrary choice rather than quoted as an optimum.')

prof.round(4)

# ================= cella 40 =================
print('>>> cella 40', flush=True)
# The sweep on the UMAP: the low end plus the annotation as the reference panel. Only the low
# resolutions are shown when the profile peaks there - the higher grids split the same structure
# further and further without a label to justify the extra clusters.
def plot_resolution_grid(columns, filename, ncols=3):
    nrows = int(np.ceil(len(columns) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5.5 * nrows))
    for col, ax in zip(columns, np.atleast_1d(axes).flat):
        sc.pl.umap(adata_plot, color=col, ax=ax, show=False, title=col,
                   legend_loc='none' if col == LABEL_FOR_PLOTS else 'on data',
                   legend_fontsize=7)
    for ax in np.atleast_1d(axes).flat[len(columns):]:
        ax.axis('off')
    fig.savefig(os.path.join(FIG_DIR, filename), dpi=300, bbox_inches='tight')
    plt.show()

# The sweep columns, sorted by resolution: the bare CLUSTER_KEY (the selected clustering, no
# numeric suffix) is left out, only the per-resolution ones are swept.
cluster_cols = sorted(
    (c for c in adata.obs.columns
     if c.startswith(f'{CLUSTER_KEY}_') and c.rsplit('_', 1)[1].replace('.', '', 1).isdigit()),
    key=lambda c: float(c.rsplit('_', 1)[1]),
)
print(f'{len(cluster_cols)} resolutions in the object: '
      f'{cluster_cols[0].rsplit("_", 1)[1]} .. {cluster_cols[-1].rsplit("_", 1)[1]}')

low = [c for c in cluster_cols if float(c.rsplit('_', 1)[1]) <= 0.5]
plot_resolution_grid(low + [LABEL_FOR_PLOTS], f'umap_leiden_resolutions_{SUFFIX}_unintegrated.png')

# ================= cella 41 =================
print('>>> cella 41', flush=True)
# Clusters against the annotation: how much of the label structure the selected clustering
# recovers, and where it does not. Row-normalized, so each cluster's row sums to 100%.
# Under `tum` a cluster that spreads evenly across the pre-CNV labels is the EXPECTED outcome, not
# a bad clustering - the labels are borrowed and the NMI behind them is ~0.34.
ct_cluster = pd.crosstab(adata.obs[CLUSTER_KEY], adata.obs[LABEL_FOR_PLOTS])
ct_cluster = ct_cluster[[c for c in label_order if c in ct_cluster.columns]]
ct_pct = ct_cluster.div(ct_cluster.sum(axis=1).replace(0, np.nan), axis=0) * 100

fig, ax = plt.subplots(figsize=(1.5 + 0.6 * ct_cluster.shape[1], 1.5 + 0.34 * ct_cluster.shape[0]))
sns.heatmap(ct_pct, cmap='Reds' if C.cell_set() == 'tum' else 'Blues', vmin=0, vmax=100,
            annot=True, fmt='.0f', annot_kws={'fontsize': 7},
            cbar_kws={'label': '% of the cluster'}, ax=ax)
ax.set_xlabel(f'{LABEL_FOR_PLOTS} ({LABEL_TITLE})')
ax.set_ylabel(f'{CLUSTER_KEY} (res {best["resolution"]:g})')
# seaborn rotates numeric-looking tick labels; cluster ids read better flat.
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
ax.set_title(f'Leiden clusters vs {LABEL_TITLE} ({SET_TITLE}), row-normalized', fontsize=10)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, f'heatmap_clusters_vs_label_{SUFFIX}_unintegrated.png'),
            dpi=300, bbox_inches='tight')
plt.show()

out = os.path.join(TABLE_DIR, f'clusters_vs_label_{SUFFIX}.csv')
ct_cluster.to_csv(out)
print('wrote', out)
ct_cluster

# ================= cella 42 =================
print('>>> cella 42', flush=True)
# Cluster composition by cohort: the counterpart of the panel above, and the one that says
# whether a leiden cluster is a biological group or one patient. On the malignant subset a cluster
# dominated by a single cohort is the expected case, not an artefact - one tumour is one clone.
ct_cohort = pd.crosstab(adata.obs[CLUSTER_KEY], adata.obs[BATCH_KEY])
top_share = (ct_cohort.max(axis=1) / ct_cohort.sum(axis=1) * 100).round(1)
dominant = ct_cohort.idxmax(axis=1)

summary_clusters = pd.DataFrame({
    'cells': ct_cohort.sum(axis=1),
    'cohorts': (ct_cohort > 0).sum(axis=1),
    'dominant cohort': dominant,
    '% from it': top_share,
})

fig, ax = plt.subplots(figsize=(9, 4.5))
ax.bar(summary_clusters.index.astype(str), summary_clusters['% from it'], color=QC_COLOR)
ax.axhline(50, color='k', ls='--', lw=1, label='50% from a single cohort')
ax.set_xlabel(CLUSTER_KEY)
ax.set_ylabel('% of the cluster from its dominant cohort')
ax.set_title(f'How patient-specific each cluster is ({SET_TITLE})')
ax.legend(frameon=False, fontsize=8)
ax.grid(axis='y', alpha=0.3)
ax.set_axisbelow(True)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, f'barplot_cluster_cohort_dominance_{SUFFIX}_unintegrated.png'),
            dpi=300, bbox_inches='tight')
plt.show()

summary_clusters

# ================= cella 44 =================
print('>>> cella 44', flush=True)
summary = pd.Series({
    'cell set': C.cell_set(),
    'cells': adata.n_obs,
    'genes': adata.n_vars,
    'cohorts': adata.obs[BATCH_KEY].nunique(),
    'post-CNV labels (cell_type)': adata.obs[LABEL_KEY].nunique(),
    'pre-CNV labels (cell_type_01_4)': adata.obs[PRIOR_KEY].nunique(),
    'HVGs': int(adata.var['highly_variable'].sum()),
    'PCs': adata.obsm['X_pca'].shape[1],
    'leiden NMI target': nmi_target,
    'leiden resolution (selected)': float(best['resolution']),
    'leiden NMI (selected)': round(float(best['score']), 4),
    'leiden clusters (selected)': adata.obs[CLUSTER_KEY].nunique(),
    'cycling cells (S + G2M)': f'{phase_counts[["S", "G2M"]].sum() / adata.n_obs:.1%}',
    'median genes per cell': int(adata.obs['n_genes_by_counts'].median()),
    'median counts per cell': int(adata.obs['total_counts'].median()),
    'median pct_counts_mt': round(float(adata.obs['pct_counts_mt'].median()), 2),
    'median cnv_score': round(float(adata.obs['cnv_score'].median()), 4),
}, name=f'{SET_TITLE} subset')

out = os.path.join(TABLE_DIR, f'summary_{SUFFIX}.csv')
summary.to_frame().to_csv(out)
print('wrote', out)
summary.to_frame()
