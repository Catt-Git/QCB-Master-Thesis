# EMT signature files (`EMT_[ABC]_{EPITHELIAL,HYBRID,MESENCHYMAL}.txt`)

Source: the collaborator's EMT/hysteresis gene lists (list A = experimentally validated,
list B = A extended with non-validated genes, list C = Tomás' list, originally mouse
symbols). Same format as every other file here: one HGNC symbol per line, no header.

The originals mix protein/CD names, complexes and heterodimers with gene symbols, so the
files below are the **transcribed** lists. Every substitution made is recorded here; nothing
was added that is not in the originals, and nothing was dropped except where noted.

## Symbol normalisation (all three lists)

| original | file | note |
|---|---|---|
| ZO-1 | `TJP1` | already present in the same list -> deduplicated |
| ZO-3 | `TJP3` | |
| OCLN-1 / CLDN-3 / CLDN-4 / SDC-1 / NID-1 / FAT-1 | `OCLN` `CLDN3` `CLDN4` `SDC1` `NID1` `FAT1` | hyphen dropped |
| MMP-14 / GSK3β / LEF-1 / CSF-1 | `MMP14` `GSK3B` `LEF1` `CSF1` | |
| CD106 / CD104 / CD51 / CD61 | `VCAM1` `ITGB4` `ITGAV` `ITGB3` | CD name -> gene |
| P-CAD | `CDH3` | |
| NRF2 | `NFE2L2` | |
| FSP1 | `S100A4` | |
| CD44+ALDH1+ | `CD44` + `ALDH1A1` | a flow-cytometry gate, split into its two genes |
| mTORC1, mTORC2 | `MTOR` | complexes, not genes; one entry |
| ITGA5:ITGB1, ITGAV:ITGB3, ITGA5:ITGB6 | `ITGA5` `ITGB1` `ITGAV` `ITGB3` `ITGB6` | heterodimers split |
| ACTA1 (aSMA) | `ACTA2` | **see below** |
| NFATc | `NFATC1` | ambiguous in the original (NFATC1/2/3/4); NFATC1 chosen by the collaborator |

## Decisions taken

**List B is the primary triad**, because it is the >=2-of-3 consensus of the three lists and
not because it is the longest - see "The three lists are nested" below. A is its experimentally
validated CORE (a strict subset), C the one semi-independent curation.

1. **ACTA1 vs ACTA2.** Lists A and B write `ACTA1 (aSMA)`, but alpha-SMA is **ACTA2**
   (ACTA1 is skeletal-muscle alpha-actin). List C writes `Acta2`. The files carry `ACTA2`,
   confirmed by the collaborator.
2. **`SNAI2` is in both the hybrid and the mesenchymal list** (B and C, as in the originals).
   It is kept in both: the two sets are not disjoint by construction, which is why the
   hybrid/mesenchymal overlap has to be reported before the two are read as independent.

## List C: mouse -> human

Uppercased, plus `Cd24a` -> `CD24`, `S100a4` -> `S100A4`, `Nfe2l2` -> `NFE2L2`,
`Mir21a` -> `MIR21`. `Krt8`/`Krt18` and `Itgb4` are listed twice in the original and were
deduplicated. `Atl1` is kept verbatim (ATL1, atlastin-1) although it looks like a typo in
the source list.

## Coverage on `shiao_epi.h5ad` (26,371 genes; HVG = the 2,000 DRVI features)

| file | n | mapped | in HVG background |
|---|---|---|---|
| EMT_A_EPITHELIAL | 15 | 15 (1.00) | 2 |
| EMT_A_HYBRID | 12 | 12 (1.00) | 6 |
| EMT_A_MESENCHYMAL | 19 | 19 (1.00) | 10 |
| EMT_B_EPITHELIAL | 28 | 28 (1.00) | 6 |
| EMT_B_HYBRID | 20 | 20 (1.00) | 6 |
| EMT_B_MESENCHYMAL | 33 | 33 (1.00) | 15 |
| EMT_C_EPITHELIAL | 32 | 29 (0.91) | 5 |
| EMT_C_HYBRID | 34 | 33 (0.97) | 12 |
| EMT_C_MESENCHYMAL | 43 | 41 (0.95) | 17 |

Unmapped, all in list C: `CLDN13` (no human ortholog) and the miRNAs `MIR200A`, `MIR34A`,
`MIR151A`, `MIR10B`, `MIR21` (not in the gene-expression object). They are kept in the files
so that 04_3's coverage table reports them rather than hiding them. `HOTAIR` and `MYOSLID`
(lncRNAs) do map.

Every list clears 04_3's 0.60 mapping floor and the 10-mapped-gene minimum. Two files fall
below 10 genes inside the 2,000-HVG ORA background (`EMT_A_EPITHELIAL`: 2, `EMT_C_EPITHELIAL`: 5),
i.e. they are effectively untestable in the factor-first route even though the cell-first
route scores them on all genes.

The primary triad inside the HVG background, which is the tested form in the factor-first route:

| file | in HVG | which |
|---|---|---|
| EMT_B_EPITHELIAL | 6/28 | KRT19, SDC1, GATA3, KLF4, CLDN3, CLDN4 |
| EMT_B_HYBRID | 6/20 | VCAM1, TNC, CDH3, CD44, PDPN, ALDH1A1 |
| EMT_B_MESENCHYMAL | 15/33 | ZEB1, ZEB2, VIM, FN1, S100A4, SPARC, ACTA2, MMP2, MMP3, MMP9, CTNNB1, LEF1, ADAM12, ITGA5, ITGB6 |

This is a limit of the factor-first route only, NOT a reason to retrain DRVI on all genes:
the cell-first route scores every list on the full 26,371-gene object.

## The three lists are nested, and B is their consensus

Measured on the transcribed files, per axis:

| axis | \|A\| | \|B\| | \|C\| | A \\ B | A n B n C | only in B | only in C | >=2 of 3 | union |
|---|---|---|---|---|---|---|---|---|---|
| epithelial  | 15 | 28 | 32 | 0 | 13 | 2 | 8  | 26 | 36 |
| hybrid      | 12 | 20 | 34 | 0 | 11 | 2 | 17 | 18 | 37 |
| mesenchymal | 19 | 33 | 43 | 0 | 19 | 1 | 11 | 32 | 44 |

Three consequences, and they change how the lists have to be described:

1. **A is a strict subset of B on every axis** (`A \ B` is empty), which is what "B = A extended
   with non-validated genes" means at the level of the files. A is also almost entirely inside C
   (13/15, 11/12, and 19/19 on the mesenchymal axis, where `A subset C` holds exactly).
2. **B carries only 2 / 2 / 1 genes that no other list has.** Removing them turns B into the
   >=2-of-3 majority-vote set (26 / 18 / 32 genes). B is therefore, numerically, the consensus of
   A, B and C - which is the reason to make it primary. "It has the higher n" is not.
3. **A is not a replicate of B and must not be called one.** A subset cannot replicate its own
   superset. The A<->B comparison is a *sensitivity analysis on the non-validated genes*: exactly
   one thing changes between them, the validation status of the extra genes, which is what makes
   it interpretable. B<->C changes curator, species of origin, size and validation status all at
   once, so a difference there is attributable to nothing in particular. C is the only comparison
   that is a replicate in the ordinary sense, and even it inherits A.

## Nesting at the gene level is NOT redundancy at the readout level

The obvious conclusion from the table - drop A, it is contained in the others - is wrong on this
data. Scores are means over the set, so adding 13-14 genes does not extend a score, it MOVES it.
Jaccard of the cell sets each version calls (`quadrant_stability_emt_drvi_epi_64.csv`):

| | list B | list A | list C |
|---|---|---|---|
| **list B** | 1.00 | **0.32** | 0.39 |
| **list A** | 0.32 | 1.00 | 0.24 |
| **list C** | 0.39 | 0.24 | 1.00 |

A is a strict subset of B at the gene level and still shares only 32% of its called cells with it -
*less* than B shares with C. The nested list is the most divergent readout of the three, not the
redundant one.

It is also frequently the strongest. Over the 128 dimension-directions of
`convergence_emt_drvi_epi_64.csv`, the readout with the largest Route A correlation is
`EMT_A_MESENCHYMAL` in **33** rows against **8** for `EMT_B_MESENCHYMAL`, which contains it
entirely (19/19 genes).

**But read that with `signature_concentration`, not on its own.** The mesenchymal lists are
VIM-dominated in inverse proportion to their size:

| list | n mapped | effective n | VIM share of variance | rho(VIM, own score) |
|---|---|---|---|---|
| EMT_A_MESENCHYMAL | 19 | 3.75 | 0.380 | 0.746 |
| EMT_B_MESENCHYMAL | 33 | 6.92 | 0.263 | 0.662 |
| EMT_C_MESENCHYMAL | 41 | 9.50 | 0.210 | 0.598 |

So A does not win because the non-validated genes dilute a real signal; it wins because it is the
purest VIM readout of the three, and VIM is precisely the gene the
`mesenchymal_may_be_ambient_or_doublet` flag exists for. The honest statement is that on the
mesenchymal axis all three lists are largely measuring VIM, and the ranking between them tracks how
much else is averaged in. That is an argument for reporting the flag, not for changing the primary.

## Pairwise Jaccard (mapped genes) - these are NOT independent readouts

A/B/C of the same axis overlap heavily by construction: B_MESENCHYMAL vs C_MESENCHYMAL 0.76,
B_EPITHELIAL vs C_EPITHELIAL 0.73, A_HYBRID vs B_HYBRID 0.63. Only one of A/B/C is primary; the
other two are a sensitivity analysis (A) and a replicate (C), and neither is an independent test.
