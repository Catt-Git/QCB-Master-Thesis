# 05_10_pipeline_recall — the standard pipeline against the factor pipeline

[05_9](../05_9_embedding_control/README.md) asked whether a programme is **present** in a
coordinate system and answered **yes for Harmony** — on most readouts more cleanly than for
DRVI. That answer is real, it is in that step's own table, and it is the reason this step
exists: presence is not the question the phase turned on.

The question this step asks is the one the phase actually rests on:

> run a pipeline on these cells, naming no axis in advance. Does it **arrive** at the
> programme?

So it does not compare two spaces. It compares two ways of turning cells into a gene list,
and then runs the identical ORA on both.

```
standard   space → k-NN graph → Leiden → one-vs-rest DE → top N genes → ORA
factor     DRVI decoder → top N genes per dimension-direction → ORA
```

Everything downstream of *"top N genes"* is one code path. The gene sets, the background, the
list depth, the hypergeometric test and the Benjamini–Hochberg treatment are the same object
in both arms, so **the only thing that differs is where the gene list came from**. That is
the whole design, and it is what lets the difference be attributed to something.

## The arms, and why there are four

| arm | what it is | what it answers |
|---|---|---|
| `harmony_leiden` | 05_3b's corrected components, clustered | the pipeline anyone would actually run |
| `drvi_leiden` | 05_3's latent space, clustered **the same way** | separates the *space* from the *unit* |
| `pca_leiden` | the uncorrected PCA, clustered | what the batch correction bought for discovery |
| `drvi_axes` | 05_3's decoder, read per dimension-direction | the phase's own route, i.e. 05_7 |

`drvi_leiden` is the arm that makes this an experiment rather than a demonstration. If it
lands with `harmony_leiden` and both land below `drvi_axes`, the loss belongs to **clustering**
and not to Harmony — which is the honest reading and the stronger claim, because it does not
require Harmony to be the worse method at anything it promises. Phase 02 is where the two were
benchmarked on what they both claim; nothing here revisits that.

## The two outcomes, and why the obvious one is useless

The first thing anyone measures is *"did some gene list of this arm enrich for this signature
at FDR < 0.05"*. On these cells that number **saturates** — it is 17/17 for the metaprogram
collection in every real arm — because a 200-gene marker list drawn from a 2,000-gene
background overlaps something in a 22-set catalogue almost always. The saturation is not a
nuisance to be tuned away; it is the first result of the step, and it is reported as
`recall_any` precisely so that it cannot be quoted as agreement between the two pipelines.

The outcome that discriminates is the one that matches what **naming a dimension means** in
[05_8](../05_8_convergence/README.md):

> a unit of analysis — a cluster or a dimension-direction — is **identified by** the signature
> it enriches for most strongly, and a programme counts as **named** by an arm when some unit
> of that arm has it as its best hit.

`recall_named` is the fraction of the testable catalogue that ends up as some unit's identity,
and it is the headline. The gap between the two numbers is the whole phenomenon: an arm whose
`recall_any` is 1.00 and whose `recall_named` is 0.59 did not find seventeen programmes — it
found ten coarse identities, each significant for a dozen overlapping lists.
`mean_hits_per_list` is the same statement from the other side.

## The target list, and why it is not "the programmes DRVI named"

Recall needs a denominator neither pipeline chose. Using 05_8's 19 named programmes would be
circular: they are the output of one of the arms. The denominator here is a property of the
**catalogue** and the **object** — every signature with at least `MIN_SIGNATURE_GENES` genes
inside the ORA background, i.e. every set either arm could in principle enrich for. It is
computed before any arm runs, it is identical for all of them, and it is printed.

It is generous in one direction and harsh in another, and both are stated rather than
corrected: generous because a signature no cell in this object expresses is in the denominator
and neither arm can find it, harsh because such a set counts against both equally.
`--presence-filter` narrows it to the signatures 05_9 found present in the Harmony space at
≥ 3 sd over the matched-random level; the headline number is the unfiltered one.

## The three things that make it a measurement and not a demonstration

**1. The DE significance filter, and the run that forced it.** The first version handed the
ORA the top 200 genes of every cluster whether or not the cluster had a marker. Under it a
**size-matched random partition scored 10/10 on `scie`, exactly like both real arms** — 200
arbitrary genes of a 2,000-gene background overlap a 100-gene signature about as often as 200
real ones do. That is the instrument failing to discriminate, not a property of clustering.
A gene now enters a cluster's list only if the one-vs-rest test calls it (`--de-fdr`, default
0.05, and `--de-min-lfc`, default 0.0 — up-regulated is enough, Seurat's own default is 0.25);
`--n-top-genes` then truncates what is left, so a cluster with a thousand markers is read at
the same depth as a decoder direction. A cluster with **no** marker contributes an empty list,
tests nothing, and is counted rather than compensated for.

**2. The null partition.** The cluster labels of the reference resolution shuffled among the
cells: the number of groups and every group size preserved, membership destroyed, through the
identical DE and ORA. It is what caught the failure above and it stays in every run.

**3. The budget objection, answered empirically rather than argued away.** The factor arm
brings 2 × `n_latent` gene lists and Leiden at resolution 0.2 brings six. More lists is more
chances to hit *and* a bigger BH denominator. The scan runs to resolution 3.0, which puts the
cluster count in the same range as the direction count, and every row of
`pipeline_recall_*.csv` carries `n_gene_lists` and `n_tests` — so recall can be read against
the budget that bought it. The right-hand panel of `recall_by_resolution` plots exactly that.
A cluster arm still below the factor arm at matched list count has not been starved of tests.

## The diagnostic that explains the gap

`eta2_cluster` is the share of a signature's **per-cell score** that lies *between* clusters
rather than across them — one-way η², computed on 05_6's cached within-cohort z scores, the
same values 05_9 measured presence on, re-scored nowhere here.

A programme with a high η² is a programme the partition has a group for. A programme with a
low one is a **gradient across** the groups, and a gradient is what a one-vs-rest DE has
nothing to test. If the programmes the cluster arm cannot name are the low-η² ones, the miss
is *explained* rather than merely counted. `eta2_vs_recovery_<arm>.png` is that plot, and a
miss on the right-hand side of it is **not** explained by this figure and must not be claimed
to be.

`cohort_dominance`, in `cluster_composition`, is the other well-known failure and is reported
rather than filtered: a cluster that is 90 % one patient makes its DE a patient contrast, and
any programme it enriches for is that patient's. Filtering it would be deciding on the
pipeline's behalf which of its clusters were allowed to count.

## Run

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
cd 05_drvi_tumoral_epi/05_10_pipeline_recall

N_LATENT=64 conda run -n benchmark-py-r --no-capture-output bash ./pipeline_recall_all.sh
N_LATENT=64 python3 pipeline_recall_tum.py --collection gavish
N_LATENT=64 python3 pipeline_recall_tum.py --collection scie --arms harmony_leiden drvi_axes
```

`N_LATENT=64` is not optional in practice: 05_3b only ever ran at 64, so at any other value
the Harmony embedding this step needs is not on disk. `CELL_SET`, `HVG_SET` and
`PRUNE_VANISHED` behave exactly as in 05_4 – 05_9; under `PRUNE_VANISHED=1` the factor arm
loses its vanished directions and everything lands in `tables_pruned/` and `figures_pruned/`.

**`gavish` is the collection to read first.** 22 metaprograms is the only one of the three wide
enough for a recall number to have room to move; `scie` has ten heavily overlapping lists and
`emt` four inside the HVG background, and both are reported for completeness rather than as
evidence.

### The cache

A Leiden partition and its DE depend on the **space**, never on the collection, so they are
computed once per arm and reused by all three collections:

```
$DATA_DIR/05_tum/pipeline_lists_<arm>_<run_id>.json        the gene lists and the DE stats
$DATA_DIR/05_tum/pipeline_partitions_<arm>_<run_id>.csv.gz the labels, one column per resolution
```

Keyed on the parameters that can change them — `k`, the list depth and the two DE filters —
and dropped **whole** when any of them moves, so a stale list cannot survive a flag change.
The partitions are realigned on cell **name**, and a cache that does not cover every cell of
the object is discarded rather than lined up on a guess. `--overwrite` recomputes.

## What comes out

`../tables/<collection>/<run_id>/`:

| file | |
|---|---|
| `pipeline_recall_*.csv` | **the headline.** One row per (arm, resolution): `recall_named`, `recall_any`, `n_gene_lists`, `n_tests`, `mean_hits_per_list`, the cluster-size and cohort-dominance summaries |
| `programme_recovery_*.csv` | one row per signature, one block of columns per arm: named, best FDR, which list, how many lists, and `eta2_cluster` |
| `gene_list_hits_*.csv` | one row per gene list: its size, its DE count, how many signatures it is significant for |
| `cluster_composition_*.csv` | one row per cluster at the reference resolution: size, cohorts, dominance, markers, which signatures it hit |
| `programme_cluster_eta2_*.csv` | η²(cluster) for every (arm, resolution, signature) |
| `ora_all_pairs_*.csv` | every pair tested, before any filter |

Figures in `../figures/05_10_pipeline_recall/<collection>/<run_id>/`:

| figure | |
|---|---|
| `recall_by_resolution` | the headline, twice: against resolution, and against the list budget |
| `recovery_matrix` | signature × arm. `*` = some unit is identified by it, `.` = significant but never any unit's best hit |
| `eta2_vs_recovery_<arm>` | the explanation: what the cluster arm could not name, against how clustered it is |
| `hits_per_gene_list` | the entanglement read: how many programmes one unit carries |

The run id in every name is the **DRVI** one, as in 05_9, and that is not a mistake: it names
the *cells* and the *gene sets*, which are the same in every arm. The `arm` column says which
pipeline each row is from.

## What the first run said

`drvi_tum_64_nomt`, 42,096 malignant cells, 19 cohorts, the 2,000-HVG panel without the MT-
genes. Eight resolutions, four arms, three collections. **Read `gavish` and treat the other
two as what they are.**

### gavish, 17 testable metaprograms of 22

| arm | best single setting | union over the whole scan | `recall_any` | signatures per unit |
|---|---|---|---|---|
| DRVI decoder, 128 directions | **17/17** | 17/17 | 17/17 | 3.80 |
| DRVI decoder, **75 directions** (budget matched) | **17/17** | 17/17 | 17/17 | 4.01 |
| DRVI + Leiden + DE | 13/17 (res 3.0, 75 clusters) | 14/17 | 17/17 | 5.20 |
| Harmony + Leiden + DE | 12/17 (res 1.5, 33 clusters) | **12/17** | 17/17 | 4.97 |
| uncorrected PCA + Leiden + DE | 12/17 (res 3.0, 65 clusters) | 13/17 | 17/17 | 5.35 |
| size-matched random partition | 0/17 | 0/17 | 0/17 | 0.00 |

**The naive metric is 1.00 for every real arm and 0.00 for the null.** `recall_any` saturates
exactly as described above and carries no information beyond "the instrument discriminates
against noise". Everything below is `recall_named`.

**The budget objection does not survive.** Cut to the 75 directions of DRVI's own
reconstruction-effect ranking — the same list count as the widest Leiden partition in the run,
and a smaller BH denominator than three of the four cluster settings it is compared against —
the decoder still names **17/17**. The clustering arms at that width name 12 and 13.

**The loss is the unit of analysis, not Harmony.** `drvi_leiden` runs the identical clustering
on DRVI's own latent space and lands at 13–14/17, beside `harmony_leiden` at 12 and the
uncorrected PCA at 12–13 — all far from the 17 the same space reaches when it is read through
the decoder instead of through clusters. Nothing here is charged to the integration method.

**Five metaprograms are never named by Harmony + Leiden at any of the eight resolutions**:
`MP13_EMT_2`, `MP15_EMT_4`, `MP18_INTERFERON_MHC_II_2`, `MP22_SECRETED_1`,
`MP3_CELL_CYCLE_HMG_RICH`. Two of them — `MP13_EMT_2` and `MP22_SECRETED_1` — are among the
**19 programmes 05_8 actually reports** ([`programme_dimensions_consensus`](../tables_pruned/consensus/drvi_tum_64_nomt/programme_dimensions_consensus_drvi_tum_64_nomt.csv):
DR 19− at tier A, DR 21+ and DR 29− at tier A). On these cells they are invisible to the
standard pipeline at every granularity it was given.

**The misses are the gradients, with one instructive exception.** At the reference resolution
the median η²(cluster) is **0.218** for what Harmony + Leiden names and **0.151** for what it
does not, and the five lowest-η² metaprograms in the collection — `MP13_EMT_2` (0.105),
`MP22_SECRETED_1` (0.127), `MP7_STRESS_IN_VITRO` (0.139), `MP39_METAL_RESPONSE` (0.150),
`MP15_EMT_4` (0.185) — are all missed. The exception is on the other side of the plot and is
the second failure rather than a counterexample: `MP3_CELL_CYCLE_HMG_RICH` has η² = **0.366**,
higher than almost anything the arm did name, and is still missed — because `MP1_CELL_CYCLE_G2_M`
and `MP2_CELL_CYCLE_G1_S` take the same clusters and win the best-hit on every one of them.
A cluster carries **4.97 signatures on average and up to 12**; the third cell-cycle programme
has nowhere to be.

**Harmony did its job, and the cluster composition says so.** At resolution 1.0, 4 of 23
Harmony clusters are more than 80 % one patient, against **29 of 34** on the uncorrected PCA
and 18 of 36 on the DRVI space. The batch correction is working, the uncorrected space
clusters patients, and the discovery gap above is not a batch artefact — the same point 05_3b
made on η²(cohort) per dimension, now at the level of the clusters themselves.

### scie and emt: the limits of the claim, stated

On `scie` the gap **disappears**: every real arm reaches 9/10 over the scan, and the
budget-matched decoder is the *worst* of them at 7/10. Ten heavily overlapping lists with a
median of 247 genes leave a 17-name vocabulary nothing to distinguish — the collection is a
readout, not a catalogue, and it was never built to be one. `emt` has **four** lists inside the
HVG background and `drvi_leiden` (3/4) beats `drvi_axes` (2/4) on it; with four sets that is
noise, and it is reported rather than dropped.

So the result is about a **wide vocabulary**. It says that on a catalogue with enough distinct
programmes to have room, the clustering pipeline produces about two thirds of them as
identities and the decoder produces all of them at the same number of tests. It does **not**
say the decoder wins on any collection, and `scie` is the counterexample in this step's own
output.

## What this step does not claim

It does not rank Harmony and DRVI as integration methods — that is phase 02, on what they both
promise — and it does not say the standard pipeline is wrong. A cluster is the right unit for
a question about discrete populations. **This object has none**: it is one compartment of one
lineage, and [05's own README](../README.md) says in advance that partial EMT in carcinoma is
a continuum and that three self-naming clusters are not what to expect. What this step
measures is what that costs, in programmes, on these cells.
