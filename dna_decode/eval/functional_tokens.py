"""Drug-general functional-unit alphabet for the functional-alphabet probe.

The hypothesis (from /hypothesise + /probe): a FUNCTION-LEVEL token alphabet (gene/allele/mechanism) might
carry resistance signal the base-level k-mer alphabet cannot. The probe asks whether function-level tokens
beat k-mer on the de-confounded WITHIN-LINEAGE metric.

This module is the functional tokenizer. It is DRUG-GENERAL: the alphabet for a drug is its curated
mechanism loci from `mic_tiers.loci_by_mechanism_for(drug)` (cipro -> QRDR alleles + plasmid qnr; tet ->
tetA/tetB efflux + tetM/tetO ribosomal-protection + tetX enzymatic; etc.). For each strain's AMRFinder run
it emits: a `gene:<locus>` token per acquired gene whose locus is in the drug's mechanism catalog
(main.tsv), the specific non-synonymous point allele on a target locus (mutations.tsv), and a coarse
`mech:<mechanism-class>` rollup. It DELIBERATELY EXCLUDES lineage-identity tokens (no `mlst:<ST>`) -- under
leave-one-MLST-out CV an ST-identity token is useless AND would be "lineage in a thin disguise".

Pure + offline: reads the committed AMRFinder caches (`data/amrfinder_runs/<accession>/`); no Docker/GPU/net.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pandas as pd

from dna_decode.data.mic_tiers import classify_gene_symbol
from dna_decode.eval.point_baseline import _is_synonymous, _locus_of


def _strip_parens(sym: str) -> str:
    """AMRFinder reports `tet(A)` / `aac(6')-Ib`; the curated catalog uses `tetA` / `aac6-Ib`. Strip the
    parentheses so the tolerant `classify_gene_symbol` matcher can resolve the mechanism."""
    return sym.replace("(", "").replace(")", "")


def _symbols(tsv: Path) -> list[str]:
    if not tsv.exists():
        return []
    with open(tsv, encoding="utf-8") as f:
        return [(row.get("Element symbol") or "").strip() for row in csv.DictReader(f, delimiter="\t")]


def strain_functional_tokens(row: pd.Series, runs_root, drug: str) -> set[str] | None:
    """Mechanism functional tokens for one strain (drug-general), or None if its AMRFinder cache is absent.

    Tokens: `gene:<locus>` (acquired gene on a mechanism locus, from main.tsv) + the specific
    non-synonymous point allele on a target locus (from mutations.tsv) + `mech:<mechanism-class>` rollups.
    No lineage-identity tokens.
    """
    d = Path(runs_root) / str(row["assembly_accession"])
    main_tsv, mut_tsv = d / "main.tsv", d / "mutations.tsv"
    if not main_tsv.exists() and not mut_tsv.exists():
        return None
    toks: set[str] = set()
    # acquired genes (main.tsv): presence per mechanism locus (tolerant match handles tet(A)->tetA etc.)
    for sym in _symbols(main_tsv):
        n = _strip_parens(sym)
        mech = classify_gene_symbol(drug, n)
        if mech:
            toks.add(f"gene:{_locus_of(n)}")
            toks.add(f"mech:{mech}")
    # point mutations (mutations.tsv): non-synonymous allele on a mechanism locus
    for sym in _symbols(mut_tsv):
        if not sym or _is_synonymous(sym):
            continue
        n = _strip_parens(sym)
        mech = classify_gene_symbol(drug, _locus_of(n))
        if mech:
            toks.add(n)                          # the specific allele, e.g. gyrA_S83L
            toks.add(f"mech:{mech}")
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
