"""
02 utils - compatibility of scib 1.1.7 with numpy 2.x / pandas 3.x

scib 1.1.7 dates from 2023 and uses a handful of APIs that numpy and pandas have
since removed. In the current environment (numpy 2.4, pandas 3.0) this surfaces
as an AttributeError halfway through the metric computation, not at import time:
the job runs for hours and then dies on the second-to-last metric.

Importing this module BEFORE scib restores the missing APIs.

    import scib_compat  # noqa: F401  (must precede the scib import)
    import scib

What is restored, and where it is needed:

- np.in1d          removed in numpy 2.4 in favour of np.isin.
                   Used by scib.preprocessing (hvg_batch) and scib.metrics.kbet.
- pd.value_counts  removed in pandas 3.0 in favour of the .value_counts() method.
                   Used by scib.metrics.graph_connectivity, .kbet and .utils.

It also fixes anndata2ri.activate/deactivate, which in recent versions no longer
register the converter the way scib expects.

If scib is ever updated, or the environment repinned to numpy<2.4 and pandas<3,
this module becomes a harmless no-op: every patch is applied only if the API is
actually missing.
"""

import os
import shutil

import numpy as np
import pandas as pd

# numpy 

if not hasattr(np, "in1d"):
    np.in1d = np.isin

# pandas

if not hasattr(pd, "value_counts"):

    def _value_counts(values, *args, **kwargs):
        """Equivalent of the module-level function removed in pandas 3.0."""
        series = values if isinstance(values, pd.Series) else pd.Series(values)
        return series.value_counts(*args, **kwargs)

    pd.value_counts = _value_counts

# rpy2 / anndata2ri


def require_r():
    """Stop immediately, with a readable message, if R cannot be reached.

    rpy2 locates R through R_HOME or through the executable on PATH. Invoking
    .../envs/<env>/bin/python directly, without activating the conda environment,
    leaves R unreachable and the anndata2ri import fails with an unhelpful
    traceback.
    """
    if not os.environ.get("R_HOME") and shutil.which("R") is None:
        raise SystemExit(
            "R is not reachable: rpy2 cannot start, and without rpy2 neither kBET "
            "nor the other metrics that go through R can be computed.\n"
            "Activate the environment before launching the script:\n"
            "    conda activate benchmark-py-r   # locally\n"
            "    conda activate catalano_env     # on the cluster"
        )


require_r()

import anndata2ri  # noqa: E402
from rpy2.robjects import conversion, default_converter  # noqa: E402


def _activate():
    conversion.set_conversion(default_converter + anndata2ri.converter)


def _deactivate():
    conversion.set_conversion(default_converter)


anndata2ri.activate = _activate
anndata2ri.deactivate = _deactivate

# pandas: positional access on Series 
#
# scib uses the `series[0]` idiom in several places to read the first element
# (for example `adata.obs[batch_key][0]` in hvg_overlap and precompute_hvg_batch).
# Up to pandas 2 this fell back to positional access on a Series with a
# non-integer index, with a FutureWarning; in pandas 3 the fallback was removed
# and the same line raises KeyError.
#
# The occurrences are scattered across several functions, so rather than rewrite
# them one by one we restore the fallback: only when the key is an integer, the
# label lookup has failed, and the index is not integer-typed. In every other
# case pandas 3 behaviour is left intact.

_series_getitem = pd.Series.__getitem__


def _series_getitem_positional_fallback(self, key):
    try:
        return _series_getitem(self, key)
    except KeyError:
        if isinstance(key, (int, np.integer)) and not pd.api.types.is_integer_dtype(
            self.index
        ):
            return self.iloc[key]
        raise


pd.Series.__getitem__ = _series_getitem_positional_fallback
