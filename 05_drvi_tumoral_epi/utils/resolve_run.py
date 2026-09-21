#!/usr/bin/env python3
"""Print `<run_id> <compartment> <out_tag>` for the current environment, one line.

Exists so that `signature_interpretation_all.sh` can build the paths it checks for without a
second, bash-side copy of two rules that already live in Python: `cell_set` owns the run id
(CELL_SET x HVG_SET x N_LATENT), `signature_common` owns the output tag ($PRUNE_VANISHED).
A driver that re-derives either in bash is how one ends up reporting [have] on files the
steps do not write.

The tag is printed as `-` when it is empty, so the line always has three fields.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "05_2_subsetting"))

import cell_set as CS              # noqa: E402
import signature_common as C       # noqa: E402

print(C.RUN_ID, CS.compartment(), C.OUT_TAG or "-")
