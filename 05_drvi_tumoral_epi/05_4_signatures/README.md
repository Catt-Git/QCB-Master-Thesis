# 05_4_signatures

The prior knowledge, ingested and characterised, before either route touches it. Counterpart of
[04_3_signatures](../../04_drvi_epithelial/04_3_signatures/), on the malignant object.

Reads the collaborator's plain-text lists from `$DATA_DIR/signatures/` — shared with 04, they
are a property of the prior knowledge and not of a compartment — and writes one `.gmt` per
collection whose description field carries the provenance string, so the same file feeds Route A
and Route B and doubles as the Appendix table.

| File | What it does |
|---|---|
| `build_signatures_tum.py` | lists → `.gmt`, coverage, Jaccard, the two figures |
| `emt_vs_gavish_tum.py` | the one cross-collection check: do the collaborator's EMT lists and Gavish's EMT metaprograms name the same genes? |
| `signature_composition_tum.py` | what is actually inside each score: which genes carry it, and are they measurable. Runs after 05_6, to join the real scores |

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
python3 build_signatures_tum.py                      # scie, the default
python3 build_signatures_tum.py --collection emt     # the EMT lists
python3 build_signatures_tum.py --collection gavish  # the 22 TNBC-relevant metaprograms
python3 build_signatures_tum.py --collection gavish --all-metaprograms  # all 41 instead
python3 build_signatures_tum.py --allow-low-coverage # report, do not stop

# the two runs this phase reports, both on the default (MT-free) panel:
N_LATENT=64 python3 build_signatures_tum.py --collection gavish              # drvi_tum_64_nomt
CELL_SET=epi N_LATENT=64 python3 build_signatures_tum.py --collection gavish # drvi_epicnv_64_nomt

N_LATENT=64 python3 emt_vs_gavish_tum.py                   # the EMT lists against MP12 - MP15
N_LATENT=64 python3 emt_vs_gavish_tum.py --all-metaprograms  # MP16 (glioma mesenchymal) too
```

`CELL_SET`, `N_LATENT` and `HVG_SET` select **which 05_3 run** this reads — `drvi_tum_32_nomt`
by default, `N_LATENT=64` for `drvi_tum_64_nomt`, `CELL_SET=epi N_LATENT=64` for
`drvi_epicnv_64_nomt`. `HVG_SET` no longer has to be spelled out: `nomt` is the default panel
now (`HVG_SET=withmt` is the explicit opt-out), and the `_nomt` tag stays in the run id because
the bare `drvi_tum_64` has already meant the with-MT panel. See
[the phase README](../README.md#which-run-05_4---05_8-read); the run id is in the name of
everything written.

## Why this is re-run and not inherited from 04

The lists are the same; **the universes they are mapped onto are not**, and both matter:

- Route A scores on all genes of the object. 05_2's `min_cells = 3` filter ran on 42,096
  malignant cells and left **25,133 genes**, against 04's 26,371 on 74,441 epithelial ones — a
  gene detected in three epithelial cells need not be detected in three aneuploid ones.
- Route B can only test the part of a signature that survived HVG selection, and 05_2
  **re-selected the HVGs inside the tumour**. Which part of a list is testable is therefore a
  different set of genes here.

The second is not a rounding difference. Against 04, inside the 2,000-HVG background:

> **These four rows are from the superseded run** — the original 05_1 CNV call and the
> with-MT gene panel — and are kept because the *direction* is the point, not the digits.
> `scie` and `emt` have not been re-run on the `newcnv` objects; their v1 tables are in
> [`../tables_v1/`](../tables_v1/). Re-measure before quoting these in the write-up.

| signature | 04 | 05 | |
|---|---|---|---|
| `ESC_WONG` | 48 | 64 | +16 |
| `BENPORATH_ES1` | 64 | 77 | +13 |
| `HALLMARK_IFNA` | 40 | 48 | +8 |
| `LIM_STEM` | 160 | 143 | **−17** |

The embryonic-stemness and interferon lists gained representation among the genes that vary
inside the tumour; `LIM_STEM`, an adult mammary stem list, lost it. Route B's power per
signature moved accordingly, in both directions, which is exactly why the `n_in_hvg_background`
column is reported and why 05_7 must not be read against 04's version of this table.

Mapped coverage on the all-genes object barely moved (`KEGG_APM` 67 → 62 is the largest change),
so Route A's inputs are comparable between the phases.

## What comes out

`../tables/<collection>/<run_id>/`:

```
signatures_<collection>_<compartment>.gmt  the collection as actually used (mapped genes)
coverage_<collection>_<run_id>.csv       per list: mapped fraction, HVG count
jaccard_<collection>_<run_id>.csv        pairwise overlap
shared_genes_<collection>_<run_id>.csv   the same as counts
```

Figures in `../figures/05_4_signatures/<collection>/<run_id>/`.

The `.gmt` is the one output named after the **object** and not after the run, because that is
what it depends on: it holds the signature genes that exist in the object, and `shiao_tum.h5ad`
(25,133 genes) and `shiao_epicnv.h5ad` (26,023) do not have the same gene axis. `N_LATENT` and
`HVG_SET` are not in the name and must not be — neither changes which genes the object has — so
`drvi_tum_32` and `drvi_tum_64` correctly share one `.gmt` while `CELL_SET=epi` gets its own.
04 has a single object and so a single `signatures_<collection>.gmt`; this is one more of the
things that differ here, and it differs because this phase has a second cell set.

Both tables have to be read before any result of this stage is believed:

- **coverage** — a symbol that does not map is **NOT MEASURED**, which is not the same as not
  expressed. These lists date from 2007–2012 and carry deprecated symbols. Below
  `MIN_MAPPED_FRACTION = 0.60` the script **stops** rather than scoring a list that is no longer
  the list it is named after (`--allow-low-coverage` to override).
- **jaccard** — the immune four are largely nested and the embryonic stemness lists overlap
  heavily, so they are not independent tests and 05_7's FDR must not be presented as if they
  were. On the EMT collection the block structure is extreme: `EMT_B_MESENCHYMAL` and
  `EMT_C_MESENCHYMAL` share a Jaccard of 0.76.

Five EMT lists sit under the 10-gene floor inside the HVG background and are effectively
untestable on Route B (`EMT_A_EPITHELIAL` has 2 genes there). They are still scored in full by
Route A, which reads the all-genes object — the asymmetry is reported, not silently resolved.

`gavish_tnbc` clears the mapping floor on both objects with room to spare; nothing had to be
forced. The lowest is `MP7_STRESS_IN_VITRO` at 0.880, and it is the only list under 0.90 —
expected, since the metaprograms were published in 2023 against a modern reference where the
SCIE lists date from 2007-2012. `MP1_CELL_CYCLE_G2_M` maps **50/50** on both objects, the only
list in the collection that loses nothing: it keeps the paper's gene symbols and this dataset is
on the reference those were written against.

What the coverage table shows instead is a **Route B** limit, and it is the same one on both
cell sets: **five of the 22** have fewer than 10 genes inside the 2,000-HVG background.

| metaprogram | genes in HVG background, `tum` | `epicnv` |
|---|---|---|
| `MP4_CHROMATIN` | 3 | 3 |
| `MP8_PROTEASOMAL_DEGRADATION` | 6 | 5 |
| `MP9_UNFOLDED_PROTEIN_RESPONSE` | 4 | 2 |
| `MP20_MYC` | 5 | 7 |
| `MP21_RESPIRATION` | 4 | 4 |

MP1 is at the other end: 47 of its 50 genes are in the `tum` HVG panel and 45 in `epicnv`, the
best-covered list here.

**A low count here is not "untestable", and the run proves it.** `MP21_RESPIRATION` has 4 genes
in the panel and still comes out *convergent* on `DR 4+` (Route A ρ 0.469, Route B FDR 2.0e-03).
The ORA depends on the overlap with each dimension's top-200 decoder genes, not on the size of
the set in the background, so this table is a caution about power and never a reason to drop a
list. They are all scored in full by Route A, which reads the all-genes object; the convergence
table of 05_8 has to be read knowing which five those are. Note what they have in common — proteostasis, transcription, MYC,
respiration: these are housekeeping-adjacent programmes whose genes are expressed everywhere
and therefore *not* highly variable, so their absence from the HVG panel is a property of the
panel and not evidence that the programme is off.


## EMT twice: the collaborator's lists against Gavish's metaprograms

`emt_vs_gavish_tum.py`. EMT enters this phase from two unrelated directions — the `emt`
collection (a **hypothesis**: nine curated lists, three axes × three curations, 05_6 calls
cells in it) and the `emt` axis of `gavish_tnbc` (a **vocabulary**: MP12 – MP15, "EMT-I" to
"EMT-IV", four NMF metaprograms recurrent across ~1,000 tumours, calling no cells). 05_8 treats
an axis that both collections land on as an axis that does not depend on whose EMT list was
used — but that argument only holds if the two are not the same genes twice. This step measures
the genes, on the mapped universe, so it is on the same footing as the Jaccard matrix above.

```
jaccard_emt_vs_gavish_<run_id>.csv        the 13 × 13 matrix
shared_genes_emt_vs_gavish_<run_id>.csv   the same as counts
jaccard_bands_emt_vs_gavish_<run_id>.csv  the four reference bands, summarised
emt_vs_gavish_pairs_<run_id>.csv          the 36 cross pairs, with the shared genes spelled out
```

**They are not the same genes twice.** Across the 36 cross pairs the median Jaccard is
**0.025**, the maximum **0.098**, and 11 pairs share no gene at all. The number to read that
against is **not 1.0**: Gavish's own four EMT metaprograms reach a median of **0.060** against
*each other* (range 0.000 – 0.099), and MP14 and MP15 share nothing. The metaprograms are
clustered NMF programmes, not curated marker panels, so four mutually near-disjoint EMT
programmes is their result and not a defect of ours. The collaborator's lists are therefore
about as close to the metaprograms as the metaprograms are to themselves, while two curations
of the *same* list reach 0.32 – 0.76. The two vocabularies are independent measurements of
EMT, which is exactly what makes a dimension they both land on worth something.

**And yet the little overlap there is, is structured — one metaprogram per axis.** Every
collaborator axis has a single best-matching metaprogram, the same one for all three curations:

| axis | best MP | shared genes |
|---|---|---|
| mesenchymal | MP12_EMT_1 | ACTA2, FN1, MMP2, SNAI2, SPARC, VIM (6/6/6 on A/B/C) |
| hybrid | MP13_EMT_2 | ITGAV, LAMC2, PDPN, TNC on A and B; C adds ITGA3, LAMB3, TGFBI |
| epithelial | MP14_EMT_3 | CD24, CLDN4, KRT8, KRT18, KRT19, LMNA on B; A has no KRT19, C no CLDN4 or LMNA |
| — | MP15_EMT_4 | nothing above 0.03; it matches no axis |

Six genes out of 25,133 is a Jaccard of 0.10 and a hypergeometric FDR of 1.7e-11 at the same
time, and those two facts are not in conflict — the `expected_by_chance` and `hypergeom_fdr`
columns of the pairs table are there so the small Jaccard is not read as "no relation".

So MP12 is the mesenchymal end, MP14 the epithelial end, MP13 the hybrid/ECM-remodelling
middle, and **MP15 is a fourth EMT programme this project has no list for**. That is the one
actionable consequence: a latent dimension whose top metaprogram is MP15 cannot be checked
against the collaborator's lists at all, and must not be reported as agreeing with them.
