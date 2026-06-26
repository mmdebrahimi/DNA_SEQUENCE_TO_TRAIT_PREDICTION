"""Functional-unit alphabet for the non-neural functional-alphabet probe.

The hypothesis (from /hypothesise + /probe): a FUNCTION-LEVEL token alphabet (codon/allele/mechanism)
might carry resistance signal the base-level k-mer alphabet cannot — k-mer top-N vocab picks the core
genome (lineage), so the probe asks whether function-level tokens beat k-mer on the de-confounded
WITHIN-LINEAGE metric.

This module is the functional tokenizer. It REUSES `point_baseline.strain_point_features` (already a
QRDR-allele + plasmid-mechanism tokenizer for cipro) and rolls up coarse mechanism-class tokens. It
DELIBERATELY EXCLUDES raw lineage-identity tokens (`mlst:<ST>`): under leave-one-MLST-out CV a held-out
lineage's ST never appears in training, so an ST-identity token is both useless AND would be "lineage in
a thin disguise" (the /probe Open-Q2 honesty point). The alphabet here is mechanism-only.

Pure + offline: reads the committed AMRFinder run caches (`data/amrfinder_runs/<accession>/`); no Docker,
no GPU, no network.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from dna_decode.eval.point_baseline import strain_point_features


def strain_functional_tokens(row: pd.Series, runs_root, drug: str) -> set[str] | None:
    """Mechanism-only functional tokens for one strain, or None if its AMRFinder cache is absent.

    Tokens = the point_baseline allele/plasmid tokens (e.g. `gyrA_S83L`, `plasmid:qnrS`) PLUS coarse
    mechanism-class rollups (`mech:qrdr`, `mech:plasmid`). NO lineage-identity tokens.
    """
    base = strain_point_features(runs_root, row["assembly_accession"], drug)
    if base is None:
        return None
    toks: set[str] = set(base)
    if any(t.startswith("plasmid:") for t in base):
        toks.add("mech:plasmid")
    if any(not t.startswith("plasmid:") for t in base):
        toks.add("mech:qrdr")  # a non-plasmid point_baseline token is a QRDR target-site allele
    return toks


def build_feature_matrix(cohort: pd.DataFrame, runs_root, drug: str):
    """Presence/absence functional-token matrix over admitted strains.

    Returns (X, vocab, strain_ids, dropped):
      X         -> (n_admitted, n_vocab) float matrix (1.0 = token present)
      vocab     -> sorted token list (deterministic column order)
      strain_ids-> admitted strain_ids, row-aligned to X
      dropped   -> strain_ids whose AMRFinder cache was absent (token set None)
    """
    per_strain: dict[str, set[str]] = {}
    dropped: list[str] = []
    for _, row in cohort.iterrows():
        toks = strain_functional_tokens(row, runs_root, drug)
        if toks is None:
            dropped.append(str(row["strain_id"]))
        else:
            per_strain[str(row["strain_id"])] = toks
    strain_ids = [s for s in cohort["strain_id"].astype(str) if s in per_strain]
    vocab = sorted({t for toks in per_strain.values() for t in toks})
    idx = {t: j for j, t in enumerate(vocab)}
    X = np.zeros((len(strain_ids), len(vocab)), dtype=float)
    for i, s in enumerate(strain_ids):
        for t in per_strain[s]:
            X[i, idx[t]] = 1.0
    return X, vocab, strain_ids, dropped
