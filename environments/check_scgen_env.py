"""Validate the scgen-py environment before it is trusted with a real run.

scGen is the only method that lives outside the main environment, and its stack
is old enough that a broken install does not always announce itself at import
time. This script fails loudly and early instead.

It checks four things, in increasing order of cost:

1. numpy is a real package, not an empty namespace directory. A partially removed
   numpy (subpackages present, ``__init__.py`` missing) still satisfies
   ``import numpy``; the failure then surfaces much later and somewhere else.
2. numpy and pandas agree on the C ABI. pandas built against numpy 1.x raises
   "numpy.dtype size changed" the moment it is imported next to numpy 2.x.
3. anndata can read an .h5ad written by the main environment (anndata 0.13.x).
   This env is two minor versions behind and has to read those files.
4. The exact scGen call path used by ``scib.integration.scgen`` runs end to end
   on a synthetic dataset: setup_anndata -> SCGEN -> train -> batch_removal.

Usage:
    conda activate scgen-py
    export DATA_DIR=~/Desktop/QCB-Master-Thesis/datasets   # optional, for check 3
    python environments/check_scgen_env.py

Exit code 0 if every check passes, 1 otherwise.
"""

import os
import sys
import traceback

# Tried in order: the file scGen will actually be given, then the small fixture.
CANDIDATE_INPUTS = ["shiao_hvg_2k.h5ad", "smoke_fixture.h5ad"]

ENCODING_REMEDY = """
        The file uses encodings introduced after anndata 0.10, so this env cannot
        read it. Writing from the main environment (pandas 3.x, anndata 0.13)
        produces three of them, and all three have to be avoided:

          nullable-string-array   pandas 3 stores the index and string columns as
                                  StringDtype -> cast them back to object dtype
          nullable-boolean        pandas BooleanDtype columns -> cast to numpy
                                  bool (only lossless when the column has no NA)
          null                    None values inside .uns, e.g. uns['log1p']
                                  ['base'] set by scanpy -> drop the key

        The fix belongs to the writer, not here. 02_1_prepare sanitises the
        objects it writes; if you hit this on a file produced elsewhere, apply the
        same three rules before write_h5ad()."""


def fail(message):
    print(f"  FAIL  {message}", flush=True)
    return False


def check_numpy_is_complete():
    """Catch the empty-namespace-package failure mode described above."""
    import numpy as np

    if getattr(np, "__file__", None) is None:
        return fail(
            "numpy imports but has no __file__: it is an empty namespace package, "
            f"not an installation. Directory: {list(getattr(np, '__path__', []))}\n"
            "        Rebuild the environment from scratch; this is not repairable "
            "in place."
        )
    if not hasattr(np, "ndarray"):
        return fail("numpy has no ndarray attribute: the installation is incomplete")

    print(f"  ok    numpy {np.__version__} ({np.__file__})", flush=True)
    return True


def check_numpy_pandas_abi():
    """pandas compiled against a different numpy ABI fails on first real use."""
    import numpy as np
    import pandas as pd

    try:
        frame = pd.DataFrame(np.zeros((3, 2), dtype=np.float32), columns=["a", "b"])
        _ = frame.to_numpy().sum()
    except Exception as error:  # noqa: BLE001 - we want the message, whatever it is
        return fail(f"numpy/pandas ABI mismatch: {type(error).__name__}: {error}")

    print(f"  ok    pandas {pd.__version__} interoperates with numpy", flush=True)
    return True


def check_reads_main_env_h5ad():
    """The h5ad files this env consumes are written by a newer anndata."""
    import anndata

    data_dir = os.environ.get("DATA_DIR")
    if not data_dir:
        print("  skip  DATA_DIR unset, cannot test .h5ad compatibility", flush=True)
        return True

    path = next(
        (p for p in (os.path.join(data_dir, n) for n in CANDIDATE_INPUTS)
         if os.path.exists(p)),
        None,
    )
    if path is None:
        print(
            f"  skip  none of {CANDIDATE_INPUTS} found under DATA_DIR",
            flush=True,
        )
        return True

    name = os.path.basename(path)
    try:
        adata = anndata.read_h5ad(path)
    except Exception as error:  # noqa: BLE001
        message = (
            f"anndata {anndata.__version__} cannot read {name}, written by the "
            f"main environment: {type(error).__name__}: {error}"
        )
        if "nullable-string-array" in str(error):
            message += ENCODING_REMEDY
        return fail(message)

    print(
        f"  ok    anndata {anndata.__version__} read {name} "
        f"({adata.n_obs} x {adata.n_vars})",
        flush=True,
    )
    return True


def check_scgen_runs():
    """Exercise the same call sequence as scib.integration.scgen."""
    import anndata
    import numpy as np
    from scgen import SCGEN

    rng = np.random.default_rng(0)
    n_cells, n_genes = 200, 50

    # Two batches with a deliberate offset, so batch_removal has something to do.
    counts = rng.poisson(2.0, size=(n_cells, n_genes)).astype(np.float32)
    counts[: n_cells // 2] += 1.0
    adata = anndata.AnnData(np.log1p(counts))
    adata.obs["batch"] = ["b0"] * (n_cells // 2) + ["b1"] * (n_cells // 2)
    adata.obs["cell_type"] = (["t0", "t1"] * (n_cells // 2))[:n_cells]
    adata.obs["batch"] = adata.obs["batch"].astype("category")
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

    try:
        SCGEN.setup_anndata(adata, batch_key="batch", labels_key="cell_type")
        model = SCGEN(adata)
        model.train(max_epochs=2, batch_size=32, early_stopping=False)
        corrected = model.batch_removal()
    except Exception:  # noqa: BLE001
        print("  FAIL  the scGen call path raised:", flush=True)
        traceback.print_exc()
        return False

    if corrected.shape != adata.shape:
        return fail(
            f"batch_removal changed the shape: {adata.shape} -> {corrected.shape}"
        )
    if not np.all(np.isfinite(corrected.X)):
        return fail("batch_removal produced non-finite values")

    import torch

    device = "GPU" if torch.cuda.is_available() else "CPU"
    print(f"  ok    scGen trained and corrected {corrected.shape} on {device}", flush=True)
    return True


def main():
    import importlib.metadata as md

    print("=" * 70, flush=True)
    print("scgen-py ENVIRONMENT CHECK", flush=True)
    print("=" * 70, flush=True)
    print(f"  python      {sys.version.split()[0]}  ({sys.executable})", flush=True)
    for package in ["numpy", "pandas", "anndata", "scanpy", "torch",
                    "pytorch-lightning", "jax", "scvi-tools", "scgen"]:
        try:
            print(f"  {package:<18} {md.version(package)}", flush=True)
        except md.PackageNotFoundError:
            print(f"  {package:<18} NOT INSTALLED", flush=True)
    print(flush=True)

    checks = [
        ("numpy installation", check_numpy_is_complete),
        ("numpy/pandas ABI", check_numpy_pandas_abi),
        ("cross-version .h5ad", check_reads_main_env_h5ad),
        ("scGen call path", check_scgen_runs),
    ]

    failed = []
    for name, check in checks:
        print(f"[{name}]", flush=True)
        try:
            passed = check()
        except Exception:  # noqa: BLE001
            print("  FAIL  the check itself raised:", flush=True)
            traceback.print_exc()
            passed = False
        if not passed:
            failed.append(name)
        print(flush=True)

    print("=" * 70, flush=True)
    if failed:
        print(f"FAILED: {', '.join(failed)}", flush=True)
        print(
            "Do not run scGen with this environment. Rebuild it from scratch:\n"
            "    conda env remove -n scgen-py\n"
            "    conda env create -f environments/scgen-py.yml",
            flush=True,
        )
        return 1

    print("All checks passed: scgen-py is usable.", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
