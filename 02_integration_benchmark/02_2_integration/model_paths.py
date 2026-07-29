"""02_2 integration: where the trained methods keep their checkpoints.

scVI, scANVI and scGen are the three methods of the phase that fit a model, and
that model is worth keeping: training is the long half of those runs and what
follows it (batch_removal for scGen, the latent pass for the others) is the
fragile half, so a crash afterwards should not cost the training again.

The checkpoint does NOT go next to the integrated object. datasets/02_integration/
holds exactly the objects the metrics read, and a model file sitting among them
would be neither an input nor an output of the benchmark. Each trained method gets
a flat directory of its own instead, matching the datasets/02_drvi/ the DRVI
notebook already writes:

    datasets/02_scvi/scvi_unscaled_model.pt
    datasets/02_scanvi/scanvi_unscaled_scvi_model.pt
    datasets/02_scanvi/scanvi_unscaled_scanvi_model.pt
    datasets/02_scgen/scgen_unscaled_model.pt
    datasets/02_scgen/scgen_scaled_model.pt

One file per checkpoint, no per-run subdirectory: scvi-tools' save() writes a
directory rather than a file, but it takes a `prefix` that is prepended to the
name inside it, so the run id goes in the file name and the directory stays flat.
That is why the callers deal in a (directory, prefix) pair. scANVI needs two
checkpoints, one per training stage, hence the two extra prefixes above.

run_all.sh and submit_integration.slurm pass the directory explicitly with
--model-dir, since they own the layout (as they do for 02_embeddings/); the prefix
follows from the output name, so it never has to be passed. This module is the
fallback for running a script by hand, and the single definition of the
convention, so the callers cannot drift apart.

Deleting a checkpoint is how a run is made to train again.
"""

from __future__ import annotations

import os


def run_id_of(output_path):
    """Grid run id of an integration output: <run_id>.h5ad -> <run_id>."""
    return os.path.splitext(os.path.basename(os.path.abspath(output_path)))[0]


def default_model_dir(output_path, method):
    """Checkpoint directory for `method`, derived from its grid output path.

    The grid writes every integration to <DATA_DIR>/02_integration/<run_id>.h5ad,
    so the data root is the output's grandparent. A path that does not follow that
    layout still lands somewhere sensible: the method directory is created beside
    whatever directory the output is in.
    """
    data_dir = os.path.dirname(os.path.dirname(os.path.abspath(output_path)))
    return os.path.join(data_dir, f"02_{method}")


def model_prefix(output_path):
    """File-name prefix for this run's checkpoint(s) inside the model directory."""
    return f"{run_id_of(output_path)}_"


def saved_model(model_dir, prefix):
    """Path scvi-tools' save(prefix=...) writes, or None if there is no directory."""
    if model_dir is None:
        return None
    return os.path.join(model_dir, f"{prefix}model.pt")


def has_saved_model(model_dir, prefix):
    """True if that checkpoint is already on disk."""
    path = saved_model(model_dir, prefix)
    return path is not None and os.path.isfile(path)
