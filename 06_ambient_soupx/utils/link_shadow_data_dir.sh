#!/usr/bin/env bash
#
# Make $DATA_DIR/06_amb usable as a DATA_DIR in its own right.
#
# This is the whole mechanism by which phase 05 can be re-run on the ambient-corrected data
# **without editing a single line of phase 05**. Every script of 05 reads its inputs as
# `$DATA_DIR/<name>` - `shiao.h5ad`, `Cells_Adult_Breast.pkl`,
# `regev_lab_cell_cycle_genes.txt`, `signatures/` - and writes its outputs under
# `$DATA_DIR/05_tum/`. Point DATA_DIR at 06_amb/ and all of that moves with it: 05 reads the
# corrected object (06_4 named it `shiao.h5ad` for exactly this reason) and writes into
# `06_amb/05_tum/`, leaving `datasets/05_tum/` untouched.
#
# The auxiliary files are symlinked rather than copied: they are inputs nobody writes to,
# and one of them (`signatures/`) is a directory phase 05 only reads.
#
# What this does NOT solve, and there is no way to solve it from here: phase 05 writes its
# FIGURES and TABLES into the repo, at `05_drvi_tumoral_epi/figures/` and
# `05_drvi_tumoral_epi/tables/`, and those paths are derived from the location of the
# scripts, not from DATA_DIR. A re-run would overwrite them. If the existing 05 figures
# matter, copy them somewhere first:
#     cp -r 05_drvi_tumoral_epi/figures 05_drvi_tumoral_epi/figures_uncorrected
#     cp -r 05_drvi_tumoral_epi/tables  05_drvi_tumoral_epi/tables_uncorrected
# `git status` will show anything that changed either way.
#
# Usage:
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets
#   ./link_shadow_data_dir.sh
#   ./link_shadow_data_dir.sh --check     # report what is there, link nothing
#
# Then, from the repo root:
#   export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets/06_amb
#   cd 05_drvi_tumoral_epi/05_1_infercnv && ./infercnv_all.sh --threads 12
#   ...and the rest of phase 05, verbatim.

set -euo pipefail

CHECK=0
[ "${1:-}" = "--check" ] && CHECK=1

: "${DATA_DIR:?set DATA_DIR to the directory holding the datasets (the datasets/ of this repo)}"
DATA_DIR="$(cd "$DATA_DIR" && pwd)"
AMB_DIR="$DATA_DIR/06_amb"

[ -d "$AMB_DIR" ] || { echo "no $AMB_DIR; run 06_2 and 06_4 first" >&2; exit 1; }

echo "source DATA_DIR : $DATA_DIR"
echo "shadow DATA_DIR : $AMB_DIR"
echo

# The corrected object itself, which 06_4 wrote. Without it the shadow directory is not a
# DATA_DIR at all, so this is checked and never linked.
if [ -s "$AMB_DIR/shiao.h5ad" ]; then
  echo "[ok]   shiao.h5ad  ($(du -h "$AMB_DIR/shiao.h5ad" | cut -f1), written by 06_4)"
else
  echo "[MISS] shiao.h5ad: run 06_4_recelltypist/recelltypist_soupx.py" >&2
fi

# Read-only inputs phase 05 expects next to it.
for item in Cells_Adult_Breast.pkl regev_lab_cell_cycle_genes.txt signatures; do
  src="$DATA_DIR/$item"
  dst="$AMB_DIR/$item"
  if [ ! -e "$src" ]; then
    echo "[MISS] $item is not in $DATA_DIR" >&2
    continue
  fi
  if [ -e "$dst" ] || [ -L "$dst" ]; then
    echo "[have] $item"
    continue
  fi
  if [ "$CHECK" -eq 1 ]; then
    echo "[link] $item  (--check: not created)"
  else
    ln -s "$src" "$dst"
    echo "[link] $item -> $src"
  fi
done

# inferCNV's gene ordering file: 05_1 downloads it if absent, so linking it is a courtesy
# and not a requirement.
GENE_ORDER="$DATA_DIR/05_tum/gene_order_hg38_gencode_v27.txt"
if [ -f "$GENE_ORDER" ]; then
  mkdir -p "$AMB_DIR/05_tum"
  dst="$AMB_DIR/05_tum/gene_order_hg38_gencode_v27.txt"
  if [ -e "$dst" ] || [ -L "$dst" ]; then
    echo "[have] 05_tum/gene_order_hg38_gencode_v27.txt"
  elif [ "$CHECK" -eq 1 ]; then
    echo "[link] 05_tum/gene_order_hg38_gencode_v27.txt  (--check: not created)"
  else
    ln -s "$GENE_ORDER" "$dst"
    echo "[link] 05_tum/gene_order_hg38_gencode_v27.txt (saves 05_1 the download)"
  fi
fi

echo
echo "next:"
echo "  export DATA_DIR=$AMB_DIR"
echo "  cd 05_drvi_tumoral_epi/05_1_infercnv && ./infercnv_all.sh"
echo "Phase 05's own figures/ and tables/ are NOT redirected; see the note in this script."
