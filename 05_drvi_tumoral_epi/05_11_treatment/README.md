# 05_11_treatment — the named programmes across the three timepoints

Every step up to 05_10 reads the malignant latent space as a whole. This one asks the
biological question of the trial: **what did the treatment do to the tumour cells?**
Shiao et al. biopsied each patient three times, before treatment (`BASE`), after one cycle
of pembrolizumab (`PD1`) and after pembrolizumab plus radiotherapy (`RTPD1`), and grouped
the patients by response (`NR`, `R1` early responders, `R2` responders only once RT was added).

Two questions, two scripts:

| script | question | unit |
|---|---|---|
| `treatment_programmes_tum.py` | does a **named programme** (a bar of `programmes_named`) go up or down? | the biopsy, paired within the patient |
| `milo_treatment_tum.py` (+ `milo_edger.R`) | does the **composition** of malignant states change, wherever it is in the space? | Milo neighbourhoods, counted per biopsy |

`treatment_common.py` holds what the two scripts must agree on: the sample definition, the
cell floor, the readouts and the hypotheses.

## Why one joint DRVI run, and not one per timepoint

`drvi_tum_64_nomt` was trained with `batch_key='cohort'`: it corrects the **patient**. The
timepoint varies inside the patient and was never given to the model, so the treatment
effect is still in the latent space, and BASE, PD1 and RTPD1 cells share one coordinate
system. Three separate runs would produce three sets of axes with no correspondence between
them. Each would also be trained on a different number of cells, so a programme missing at
RTPD1 could mean that the treatment removed it or that the smaller run had too few cells to
learn it, and nothing would tell the two apart.

**Do not** retrain with `treatment` as a covariate: that would remove exactly the signal
this step measures.

## The unit is the biopsy

A test over 42,000 cells calls everything significant, and it is pseudoreplication: the
cells of one biopsy are not independent. The sample is `cohort × treatment` (one biopsy, its
CD45+ and CD45− libraries pooled, since they come from the same piece of tumour). A biopsy
needs `--min-cells` malignant cells (default 20) to be summarised. Every contrast is paired
within the patient, so each patient serves as their own control.

## `treatment_programmes_tum.py`

- **Readouts.** Each named direction (`DR 18-`) is oriented so that higher means more of
  the programme and z-scored over all cells, so deltas are in sd. A **programme** readout is
  the mean of its directions: one hypothesis, one readout.
- **Test.** For each of `PD1_vs_BASE`, `RTPD1_vs_BASE` and `RTPD1_vs_PD1`, a Wilcoxon
  signed-rank test on the per-patient difference of the biopsy means. Fewer than 6 pairs →
  not tested, because with 5 pairs the test cannot reject.
- **Two BH families.** The hypotheses in `HYPOTHESES` were fixed before looking at the data.
  They get their own correction:

  | programme | contrast | expected |
  |---|---|---|
  | HALLMARK_IFNG, MP17_INTERFERON_MHC_II_1 | PD1 vs BASE | up (IFN-γ from re-activated T cells) |
  | HALLMARK_IFNA | RTPD1 vs PD1 | up (RT-induced type-I IFN) |
  | MP19_EPITHELIAL_SENESCENCE | RTPD1 vs PD1 | up (RT-induced senescence) |
  | MP1 / MP2 cell cycle | RTPD1 vs PD1 | down (RT-induced arrest) |

  Everything else is one exploratory family, and a hit there is a lead.
- **Gene-level check** (a priori programmes only). The same test on 05_6's raw signature
  score, with no latent space in between. A DR that moves while its genes do not is a
  property of the embedding.
- **Response** is described rather than tested as an interaction, because there are 4–8
  patients per group. The output gives the median delta per group with its n, plus an
  exploratory Kruskal–Wallis. The figures draw every patient.

## `milo_treatment_tum.py`

The pipeline follows miloR. It builds a k-NN graph on the DRVI latent, picks neighbourhood
indices by refined sampling, counts cells per neighbourhood per biopsy, fits edgeR
quasi-likelihood with `~ cohort + treatment`, and applies the spatial FDR (weighted BH). A
second model, `~ cohort + response:timepoint`, is written as exploratory.

Each neighbourhood is then read by the programmes: the mean of each programme readout over
its cells, plus its `top_programme`. The association table gives the Spearman between a
contrast's logFC and each programme, which answers "are the neighbourhoods that grow after
RT the senescent ones?". There is no p-value, because the neighbourhoods overlap.

Only malignant cells are in the graph, so a logFC is a change in **share of the biopsy's
tumour cells**. A tumour that shrinks uniformly moves nothing here; the shrinkage itself is
in `treatment_cell_counts`.

`--k` defaults to 100. Milo's rule of thumb is a mean neighbourhood size of about 5 × the
number of biopsies, and with a few hundred malignant cells per biopsy a small k leaves
counts of 0–1 per biopsy. The run prints where it lands against the rule.

## Inputs

| file | written by | where |
|---|---|---|
| `embed_drvi_tum_64_nomt.h5ad` (latent, `obs` with `cohort`/`treatment`/`response`, `X_umap`) | 05_3 | `$DATA_DIR/05_tum/` |
| `programme_dimensions_consensus_drvi_tum_64_nomt.csv` | 05_8 `dr_naming_tum.py` | `tables_pruned_rho030/consensus/drvi_tum_64_nomt/` (in git) |
| `signature_scores_scie_drvi_tum_64_nomt.csv`, `signature_scores_gavish_tnbc_drvi_tum_64_nomt.csv` — optional, gene-level check only | 05_6 | `$DATA_DIR/05_tum/` |

Environment: `benchmark-py-r`, which provides anndata, sklearn, scipy and statsmodels, and
whose R has `edgeR` and `limma`. `milo_edger.R` is run with the `Rscript` next to the active
Python, unless `--rscript` or `$RSCRIPT` says otherwise.

## Run

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
cd 05_drvi_tumoral_epi/05_11_treatment
N_LATENT=64 HVG_SET=nomt python3 treatment_programmes_tum.py
N_LATENT=64 HVG_SET=nomt python3 milo_treatment_tum.py            # k=100, prop=0.05
N_LATENT=64 HVG_SET=nomt python3 milo_treatment_tum.py --k 50     # sensitivity to k
```

The programme table defaults to `PRUNE_VANISHED=1 OUT_TAG=_pruned_rho030`, the run the
thesis reports. Export either variable to read another naming table; the outputs follow
the same tag.

## Outputs

`tables_pruned_rho030/consensus/drvi_tum_64_nomt/`:

- `treatment_cell_counts_*`: malignant cells per patient × timepoint
- `treatment_pseudobulk_*`: biopsy means and the fraction above 1 sd, per readout
- `treatment_paired_tests_*`: one row per readout × contrast (family, FDR, per-response medians, gene-level columns)
- `milo_nhoods_k<k>_*`, `milo_summary_k<k>_*`, `milo_programme_association_k<k>_*`

`figures_pruned_rho030/05_11_treatment/consensus/drvi_tum_64_nomt/`:

- `sample_cell_counts`: the denominator of everything else
- `a_priori_trajectories`: one line per patient, response medians in bold, the test underneath
- `programme_trajectories`: the same for every named programme
- `paired_delta_programmes`, `paired_delta_directions`: median paired change, a priori cells outlined
- `milo_nhood_umap`, `milo_beeswarm`, `milo_programme_association`

The edgeR working files (counts, per-contrast fits) go to
`$DATA_DIR/05_tum/milo_<run><tag>_k<k>_p<prop>_min<n>_s<seed>/`.

## Limits

- **One biopsy per timepoint.** The timepoint cannot be separated from that biopsy's
  library or from where the needle went.
- **Survivor bias.** Responders lose tumour, so a late biopsy's malignant cells are the
  ones that survived. A shift can be selection as much as reprogramming, and the two
  scripts do not distinguish them. The cell count table shows where this applies.
- **Power.** About 19 patients, fewer with both biopsies of a contrast, and 4–8 per
  response group. A null result in a response group is weak evidence of no effect.
- **The programmes were named on all timepoints together.** A state that exists only after
  treatment, and was too rare to earn a direction, is invisible to the first script. Milo
  is the partial answer to that.
