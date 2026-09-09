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
| `signature_composition_tum.py` | what is actually inside each score: which genes carry it, and are they measurable. Runs after 05_6, to join the real scores |

```bash
export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
python3 build_signatures_tum.py                      # scie, the default
python3 build_signatures_tum.py --collection emt     # the EMT lists
python3 build_signatures_tum.py --collection gavish  # the 27 TNBC-relevant metaprograms
python3 build_signatures_tum.py --collection gavish --all-metaprograms  # all 40 instead
python3 build_signatures_tum.py --allow-low-coverage # report, do not stop

N_LATENT=64 python3 build_signatures_tum.py --collection gavish   # the drvi_tum_64 run
CELL_SET=epi HVG_SET=nomt N_LATENT=64 \
    python3 build_signatures_tum.py --collection gavish           # drvi_epicnv_64_nomt
```

`CELL_SET`, `N_LATENT` and `HVG_SET` select **which 05_3 run** this reads — `drvi_tum_32` by
default, `N_LATENT=64` for `drvi_tum_64`, `CELL_SET=epi HVG_SET=nomt N_LATENT=64` for
`drvi_epicnv_64_nomt`. See [the phase README](../README.md#which-run-05_4---05_8-read); the
run id is in the name of everything written.

## Why this is re-run and not inherited from 04

The lists are the same; **the universes they are mapped onto are not**, and both matter:

- Route A scores on all genes of the object. 05_2's `min_cells = 3` filter ran on 36,192
  malignant cells and left **24,779 genes**, against 04's 26,371 on 74,441 epithelial ones — a
  gene detected in three epithelial cells need not be detected in three aneuploid ones.
- Route B can only test the part of a signature that survived HVG selection, and 05_2
  **re-selected the HVGs inside the tumour**. Which part of a list is testable is therefore a
  different set of genes here.

The second is not a rounding difference. Against 04, inside the 2,000-HVG background:

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
(24,779 genes) and `shiao_epicnv.h5ad` (26,379) do not have the same gene axis. `N_LATENT` and
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

All three collections pass the mapping floor on this object; nothing had to be forced. On
`gavish` (all 40, which is what these numbers were measured on) the lowest is `MP24_CILIA` at 0.796 and nothing else is under 0.83, which is expected:
the metaprograms were published in 2023 against a modern reference, where the SCIE lists date
from 2007-2012. What the coverage table shows there instead is a **Route B** limit — nine of the
40 have fewer than 10 genes inside the 2,000-HVG background and are effectively untestable by
the ORA of 05_7, `MP4_CHROMATIN`, `MP20_MYC` and `MP21_RESPIRATION` among them. They are still
scored per cell by Route A; a dimension can only fail to match them on one route, and the
convergence table has to be read knowing which nine those are.
